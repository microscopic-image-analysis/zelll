//! RUSTFLAGS="-C target-cpu=native" cargo bench --features quick_bench,rayon -- Iteration
use criterion::{
    AxisScale, BenchmarkId, Criterion, PlotConfiguration, SamplingMode, criterion_group,
    criterion_main,
};
use nalgebra::{Point, Point3, Vector3, distance_squared};
use rand::distributions::Standard;
use rand::prelude::*;
#[cfg(feature = "rayon")]
use rayon::ThreadPoolBuilder;
#[cfg(feature = "rayon")]
use zelll::rayon::ParallelIterator;
use zelll::{CellGrid, Particle};

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

pub fn bench_iters(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("Iteration");
    group
        .sampling_mode(SamplingMode::Flat)
        .plot_config(plot_config.clone());

    #[cfg(feature = "quick_bench")]
    group.sample_size(10);

    for size in (1..=8).map(|exp| 10usize.pow(exp)) {
        let cutoff: F32or64 = 10.0;
        let conc = 10.0 / cutoff.powi(3); //i.e. 10mol per cutoff^3 volume units
        let a = 3.0 * cutoff;
        let b = 3.0 * cutoff;
        let c = (size as F32or64 / conc) / a / b;

        let pointcloud = generate_points_random(size, [a, b, c], [0.0, 0.0, 0.0], None);
        // pointcloud.sort_unstable_by(|p, q| p.z.partial_cmp(&q.z).unwrap());
        let cg = CellGrid::new(
            pointcloud
                .iter()
                .map(|p| p.coords)
                .map(Particle::from)
                .enumerate(),
            cutoff,
        );

        #[cfg(not(feature = "rayon"))]
        group.bench_with_input(BenchmarkId::new("sequential", size), &cg, |b, cg| {
            let cutoff_squared = cutoff.powi(2);
            b.iter(|| {
                cg.particle_pairs()
                    .filter(|&((_i, p), (_j, q))| {
                        distance_squared(&(*p).into(), &(*q).into()) <= cutoff_squared
                    })
                    .for_each(|_| {});
            })
        });

        #[cfg(feature = "rayon")]
        for nthreads in 1..=16 {
            let pool = ThreadPoolBuilder::new()
                .num_threads(nthreads)
                .build()
                .unwrap();

            group.bench_with_input(
                BenchmarkId::new(format!("{} threads", pool.current_num_threads()), size),
                &cg,
                |b, cg| {
                    let cutoff_squared = cutoff.powi(2);
                    b.iter(|| {
                        pool.install(|| {
                            cg.par_particle_pairs().for_each(|((_i, p), (_j, q))| {
                                if distance_squared(&(*p).into(), &(*q).into()) <= cutoff_squared {
                                } else {
                                }
                            })
                        });
                    })
                },
            );
        }
    }

    group.finish();
}

criterion_group!(benches, bench_iters);
criterion_main!(benches);
