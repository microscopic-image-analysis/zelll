use clap::{Parser, Subcommand};
use nalgebra::{Matrix4, Point3, Vector3};
use nuts_rs::{Chain, CpuMath, DiagGradNutsSettings, Settings};
use pdbtbx::{Atom, Model, PDB, StrictnessLevel, open, save};
use psssh::io::PointCloud;
use psssh::sdf::SmoothDistanceField;
use std::path::PathBuf;
use std::time::Instant;
use zelll::Particle;

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
        /// Number of samples to produce.
        #[arg(short = 'n', long = "samples", default_value_t = 2000)]
        n: usize,
        /// Number of samples to discard before sampling 'n' samples.
        #[arg(short = 'b', long = "burn-in", default_value_t = 1000)]
        b: usize,
        /// Distance to the protein structure at which the surface will be sampled.
        #[arg(short = 'l', long, default_value_t = 1.05)]
        surface_level: f64,
        /// Force constant used for sampling on the protein surface.
        /// Smaller values might work better with smaller cutoff radii
        /// but might also require adjusting the surface level.
        #[arg(short = 'f', long, default_value_t = 10.0)]
        force_constant: f64,
        /// Maximum tree depth for NUTS. The default value is robust enough for this application.
        /// Lower values are cheaper and may suffice if it's not required to cover the complete
        /// surface with the sampled points or the sample size is large enough.
        #[arg(short = 'd', long, default_value_t = 7)]
        nuts_depth: u64,
    },

    /// Evaluate SDF of a PDB file for a specified grid of q=l³ query points
    /// in a gr
    Eval {
        /// Protein structure file to generate SDF from.
        pdb: PathBuf,
        /// Neighborhood cutoff treshold used for sampling.
        #[arg(short, long, default_value_t = 10.0)]
        cutoff: f64,
        /// Axis length for the query grid. Results in q=l³ query points
        #[arg(short = 'l', long = "axis-length", default_value_t = 256)]
        l: usize,
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
            force_constant,
            nuts_depth,
        } => {
            let out = out.clone().unwrap_or(pdb.with_extension("psssh.pdb"));

            let (data, _) = open(pdb.to_str().expect("Expected a valid file path"))
                .expect("Expected a valid PDB file");
            let data = PointCloud::from_pdb_atoms(data.atoms());

            let sdf = SmoothDistanceField::new(&data, cutoff.abs())
                .with_surface_radius(*surface_level)
                .with_k_force(*force_constant);

            let mut settings = DiagGradNutsSettings {
                num_tune: 1000,
                maxdepth: *nuts_depth,
                ..Default::default()
            };
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
            let init = data.points.first().map_or([0.0; 3], |atom| atom.coords());

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
                out.to_str().expect("Expected a valid file path"),
                StrictnessLevel::Loose,
            )
            .expect("Saving sampled structure to PDB file failed");
        }

        // We're not using criterion for this purely for convenience
        // i.e., we want to time `n` queries for PDB files of varying sizes
        // which is a bit annoying to set up in a criterion benchmark
        // Using hyperfine is also not ideal because we'd like to avoid measuring file IO
        Commands::Eval { pdb, cutoff, l } => {
            let (pdbf, _) = open(pdb.to_str().expect("Expected a valid file path"))
                .expect("Expected a valid PDB file");
            let data = PointCloud::from_pdb_atoms(pdbf.atoms());

            let sdf = SmoothDistanceField::new(&data, cutoff.abs());

            let inf = Vector3::from(sdf.grid().info().bounding_box().inf());
            let sup = Vector3::from(sdf.grid().info().bounding_box().sup());
            // bounding box volume
            let vol = (sup - inf).product();
            let transform = Matrix4::new_nonuniform_scaling(&(sup - inf)).append_translation(&inf);

            let unit_interval: Vec<_> = (0..*l).map(|i| i as f64 / (*l - 1) as f64).collect();
            let unit_grid = unit_interval.iter().flat_map(|x| {
                unit_interval
                    .iter()
                    .flat_map(|y| unit_interval.iter().map(|z| [*x, *y, *z]))
            });

            let grid: Vec<_> = unit_grid
                .map(|p| Point3::from(p))
                .map(|p| transform.transform_point(&p))
                .collect();

            let t = Instant::now();

            grid.into_iter().for_each(|p| {
                sdf.evaluate(p.into());
            });

            let t = t.elapsed().as_nanos();

            // println!("PDB, n, cutoff, queries, ns");
            println!(
                "{:?}, {}, {}, {}, {}, {}",
                pdbf.identifier
                    .clone()
                    .or(pdb.file_stem().map(|f| f.display().to_string())),
                pdbf.atom_count(),
                vol,
                cutoff,
                l.pow(3),
                t
            );
        }
    }
}
