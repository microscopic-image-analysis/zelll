use clap::{Parser, Subcommand};
use nuts_rs::{Chain, CpuMath, DiagGradNutsSettings, Settings};
use pdbtbx::{Atom, Model, PDB, StrictnessLevel, open, save};
use psssh::Angstrom;
use psssh::io::PointCloud;
use psssh::sdf::SmoothDistanceField;
use psssh::utils::approx_geodesic_dist;
// use rand::prelude::*;
use csv::Writer;
use std::path::PathBuf;
use zelll::{CellGrid, Particle};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Sample points on the surface of a single protein structure.
    Sample {
        /// Protein structure file to sample on.
        pdb: PathBuf,
        /// File path to save sampled surface to. defaults to input path + ".psssh.pdb".
        out: Option<PathBuf>,
        /// Neighborhood cutoff treshold used for sampling.
        #[arg(short, long, default_value_t = 10.0)]
        cutoff: f64,
        /// Number of samples to produces.
        #[arg(short = 'n', long = "samples", default_value_t = 2000)]
        n: usize,
        /// Number of samples to discard before sampling 'n' samples
        #[arg(short = 'b', long = "burn-in", default_value_t = 1000)]
        b: usize,
        /// Nistance to the protein structure at which the surface will be sampled.
        #[arg(short = 'l', long, default_value_t = 1.05)]
        surface_level: f64,
        /// Maximum tree depth for NUTS. The default value is robust enough for this application.
        /// Lower values are cheaper and may suffice if it's not required to cover the complete
        /// surface with the sampled points or the sample size is large enough.
        #[arg(short = 'd', long, default_value_t = 7)]
        nuts_depth: u64,
    },
    /// Analyze a sampled surface via the smooth distance function of a protein structure.
    /// The SDF is used to generate normal vectors for each surface sample in order
    /// to approximate the geodesic distance on the surface.
    Stats {
        /// The protein structure from which the surface was sampled.
        structure: PathBuf,
        /// Surface sample file.
        surface: PathBuf,
        /// CSV file path to write stats to. defaults to structure path + ".csv".
        csv: Option<PathBuf>,
        /// Cutoff radius determining the local neighborhood in which approximate geodesic distances
        /// should be measured.
        #[arg(short, long, default_value_t = 5.0)]
        cutoff: f64,
    },
}

fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Commands::Sample {
            pdb,
            out,
            cutoff,
            n,
            b,
            surface_level,
            nuts_depth,
        } => {
            let out = out.clone().unwrap_or(pdb.with_extension("psssh.pdb"));

            let (data, _) = open(&pdb.to_str().expect("Expected a valid file path"))
                .expect("Expected a valid PDB file");
            let data = PointCloud::from_pdb_atoms(data.atoms());

            let sdf =
                SmoothDistanceField::new(&data, cutoff.abs()).with_surface_radius(*surface_level);

            let mut settings = DiagGradNutsSettings::default();
            settings.num_tune = 1000;
            settings.maxdepth = *nuts_depth;
            // sometimes, NUTS gets stuck with very small step sizes
            // this *occasionally* helps in those cases
            settings.adapt_options.dual_average_options.initial_step = 0.1;

            let chain = 0;
            let math = CpuMath::new(sdf);
            let mut rng = rand::rng();
            let mut sampler = settings.new_chain(chain, math, &mut rng);

            // let init = data
            //     .points
            //     .choose(&mut rng)
            //     .map_or([0.0; 3], |atom| atom.coords());
            // TODO: alternatively, let's just use the first atom, assuming it's at one of the ends
            // TODO: we're discarding the first `b` samples anyway
            // FIXME: this seems to help against NUTS ending up with extremely small step sizes
            let init = data.points.get(0).map_or([0.0; 3], |atom| atom.coords());

            sampler
                .set_position(init.as_slice())
                .expect("Unrecoverable error during init");
            let mut trace = vec![];

            // burn-in period
            for _ in 0..*b {
                let (_draw, _info) = sampler.draw().expect("Unrecoverable error during sampling");
            }

            for _ in 0..*n {
                let (draw, _info) = sampler.draw().expect("Unrecoverable error during sampling");
                trace.push(draw);
            }

            let mut out_pdb = PDB::new();
            let mut model_out = Model::new(0);

            let atoms = trace.iter().enumerate().filter_map(|(i, coords)| {
                let [x, y, z]: [f64; 3] = coords[..].try_into().ok()?;
                Atom::new(false, i, "H", x, y, z, 1.0, 0.0, "H", 0)
            });

            for (i, atom) in atoms.enumerate() {
                model_out.add_atom(atom, "X", (i as isize, None), ("PSH", None));
            }

            out_pdb.add_model(model_out);

            save(
                &out_pdb,
                &out.to_str().expect("Expected a valid file path"),
                StrictnessLevel::Loose,
            )
            .expect("Saving sampled structure to PDB file failed");
        }
        Commands::Stats {
            structure,
            surface,
            csv,
            cutoff,
        } => {
            let csv_path = csv.clone().unwrap_or(surface.with_extension("csv"));
            let (structure, _) = open(&structure.to_str().expect("Expected a valid file path"))
                .expect("Expected a valid PDB file");
            let (surface, _) = open(&surface.to_str().expect("Expected a valid file path"))
                .expect("Expected a valid PDB file");

            let structure = PointCloud::from_pdb_atoms(structure.atoms());
            let surface = PointCloud::from_pdb_atoms(surface.atoms());
            let sdf = SmoothDistanceField::new(&structure, cutoff.abs());

            let surface_grid = CellGrid::new(surface.points.iter().copied(), *cutoff);

            let mut wtr = Writer::from_path(csv_path).expect("Could not create Writer from path");
            wtr.write_record(&["i", "j", "sd_i", "sd_i", "approx. geodesic distance"])
                .expect("Could not write CSV header");

            surface_grid
                .particle_pairs()
                .filter_map(|((i, p), (j, q))| {
                    let pc: [Angstrom; 3] = p.coords();
                    let qc: [Angstrom; 3] = q.coords();

                    let (p_sd, p_normal) = sdf.evaluate(pc)?;
                    let (q_sd, q_normal) = sdf.evaluate(qc)?;

                    let agd = approx_geodesic_dist(p, q, p_normal, q_normal);

                    if agd <= *cutoff {
                        Some((i, j, p_sd, q_sd, agd))
                    } else {
                        None
                    }
                })
                .for_each(|rec| {
                    wtr.serialize(rec).expect("Could not write CSV record");
                });
        }
    }
}
