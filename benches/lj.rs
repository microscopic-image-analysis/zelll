//! RUSTFLAGS="-C target-cpu=native" cargo bench --features quick_bench -- Lennard-Jones
use criterion::{
    AxisScale, BenchmarkId, Criterion, PlotConfiguration, SamplingMode, criterion_group,
    criterion_main,
};
use nalgebra::{Point, Point3, Vector3, distance_squared};
use rand::distributions::Standard;
use rand::prelude::*;
use zelll::{CellGrid, WrappedParticle};

type F32or64 = f64;

type PointCloud<const N: usize> = Vec<Point<F32or64, N>>;
/// Generate a uniformly random 3D point cloud of size `n` in a cuboid of edge lengths `vol` centered around `origin`.
fn generate_points_random(
    n: usize,
    vol: [F32or64; 3],
    origin: [F32or64; 3],
    seed: Option<u64>,
) -> PointCloud<3> {
    // with fixed seed for reproducability
    let mut rng = StdRng::seed_from_u64(seed.unwrap_or(3079380797442975911));

    std::iter::repeat_with(|| {
        Point3::<F32or64>::from(
            (Vector3::from_iterator((&mut rng).sample_iter(Standard))
                - Vector3::new(0.5, 0.5, 0.5)
                + Vector3::from(origin))
            .component_mul(&Vector3::from(vol)),
        )
    })
    .take(n)
    .collect()
}

// cf. https://docs.lammps.org/pair_lj.html
// and https://docs.lammps.org/units.html
// and https://github.com/lammps/lammps/blob/develop/src/pair_lj_cut.cpp
// for dimensionless `lj/cut`
// we use the squared euclidean distance as input here
// panics if dsq is zero
fn lj(dsq: F32or64) -> F32or64 {
    // (1/r)⁶
    let tmp = dsq.recip().powi(3);
    // 4*((1/r)¹² - (1/r)⁶) (with epsilon, sigma = 1)
    4.0 * tmp * (tmp - 1.0)
}

pub fn bench_lj(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("Lennard-Jones");
    group
        .sampling_mode(SamplingMode::Flat)
        .plot_config(plot_config.clone());

    #[cfg(feature = "quick_bench")]
    group.sample_size(10);

    for size in (2..=8).map(|exp| 10usize.pow(exp)) {
        let cutoff: F32or64 = 10.0;
        let conc = 10.0 / cutoff.powi(3); //i.e. 10mol per cutoff^3 volume units
        let a = 3.0 * cutoff;
        let b = 3.0 * cutoff;
        let c = (size as F32or64 / conc) / a / b;

        let pointcloud = generate_points_random(size, [a, b, c], [0.0, 0.0, 0.0], None);
        // pointcloud.sort_unstable_by(|p, q| p.z.partial_cmp(&q.z).unwrap());

        // using a smaller scope here so the allocated memory gets dropped again
        // so Massif doesn't keep it in its snapshots the whole time
        {
            let cutoff_squared = cutoff.powi(2);
            let cg = CellGrid::new(
                pointcloud
                    .iter()
                    .map(|p| p.coords)
                    .map(WrappedParticle::from),
                cutoff,
            );
            let potential_energy: F32or64 = cg
                .particle_pairs()
                .filter_map(|((_i, p), (_j, q))| {
                    let dsq = distance_squared(&(*p).into(), &(*q).into());
                    if dsq < cutoff_squared {
                        Some(dsq)
                    } else {
                        None
                    }
                })
                .map(|dsq| lj(dsq))
                .sum();
            dbg!(potential_energy / size as F32or64);
        }

        group.bench_with_input(
            BenchmarkId::new("sequential", size),
            &pointcloud,
            |b, pointcloud| {
                b.iter(|| {
                    let cutoff_squared = cutoff.powi(2);
                    let cg = CellGrid::new(
                        pointcloud
                            .iter()
                            .map(|p| p.coords)
                            .map(WrappedParticle::from),
                        cutoff,
                    );
                    let _potential_energy: F32or64 = cg
                        .particle_pairs()
                        .filter_map(|((_i, p), (_j, q))| {
                            let dsq = distance_squared(&(*p).into(), &(*q).into());
                            if dsq < cutoff_squared {
                                Some(dsq)
                            } else {
                                None
                            }
                        })
                        .map(|dsq| lj(dsq))
                        .sum();
                    // dbg!(_potential_energy / size as F32or64);
                })
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_lj);
criterion_main!(benches);
