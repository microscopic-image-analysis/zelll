use criterion::{
    black_box, criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion,
    PlotConfiguration, SamplingMode,
};
use nalgebra::{distance_squared, Norm, Point, Point3, UniformNorm, Vector3};
use rand::distributions::Standard;
use rand::prelude::*;
use zelll::cellgrid::CellGrid;
#[cfg(feature = "rayon")]
use zelll::cellgrid::ParallelIterator;

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

pub fn bench_cellgrid_concentration(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("CellGrid");
    group
        .sampling_mode(SamplingMode::Flat)
        .plot_config(plot_config.clone());

    group.bench_with_input(BenchmarkId::new("Gonnet2007", 1000), &1000, |b, size| {
        let pointcloud = generate_points_random(*size, [3.166; 3], [0.0, 0.0, 0.0], None);
        let cg = CellGrid::new(pointcloud.iter().map(|p| p.coords), 1.0);
        b.iter(|| {
            cg.point_pairs()
                .filter(|&((_i, p), (_j, q))| distance_squared(&p.into(), &q.into()) <= 1.0)
                .for_each(|_| {});
        })
    });

    #[cfg(feature = "rayon")]
    group.bench_with_input(
        BenchmarkId::new("Gonnet2007_par", 1000),
        &1000,
        |b, size| {
            let pointcloud = generate_points_random(*size, [3.166; 3], [0.0, 0.0, 0.0], None);
            let cg = CellGrid::new(pointcloud.iter().map(|p| p.coords), 1.0);
            b.iter(|| {
                cg.par_point_pairs()
                    .filter(|&((_i, p), (_j, q))| distance_squared(&p.into(), &q.into()) <= 1.0)
                    .for_each(|_| {});
            })
        },
    );

    for size in (2..=7).map(|exp| 10usize.pow(exp)) {
        let cutoff: F32or64 = 10.0;
        let conc = 10.0 / cutoff.powi(3); //i.e. 10mol per 10^3 volume units
        let a = 3.0 * cutoff;
        let b = 3.0 * cutoff;
        let c = (size as F32or64 / conc) / a / b;
        let vol_edges = (size as F32or64 / conc).cbrt();
        let pointcloud = generate_points_random(size, [a, b, c], [0.0, 0.0, 0.0], None);

        group.bench_with_input(
            BenchmarkId::new("::new()", size),
            &pointcloud,
            |b, pointcloud| {
                b.iter(|| {
                    CellGrid::new(pointcloud.iter().map(|p| p.coords), cutoff);
                })
            },
        );
    }

    for size in (2..=7).map(|exp| 10usize.pow(exp)) {
        let cutoff: F32or64 = 10.0;
        let conc = 10.0 / cutoff.powi(3); //i.e. 10mol per cutoff^3 volume units
        let a = 3.0 * cutoff;
        let b = 3.0 * cutoff;
        let c = (size as F32or64 / conc) / a / b;
        let vol_edges = (size as F32or64 / conc).cbrt();

        let pointcloud = generate_points_random(size, [a, b, c], [0.0, 0.0, 0.0], None);
        let cg = CellGrid::new(pointcloud.iter().map(|p| p.coords), cutoff);

        group.bench_with_input(BenchmarkId::new("::point_pairs()", size), &cg, |b, cg| {
            let cutoff_squared = cutoff.powi(2);
            b.iter(|| {
                cg.point_pairs()
                    .filter(|&((_i, p), (_j, q))| {
                        distance_squared(&p.into(), &q.into()) <= cutoff_squared
                    })
                    .for_each(|_| {});
            })
        });

        #[cfg(feature = "rayon")]
        group.bench_with_input(
            BenchmarkId::new("::par_point_pairs()", size),
            &cg,
            |b, cg| {
                let cutoff_squared = cutoff.powi(2);
                b.iter(|| {
                    cg.par_point_pairs().for_each(|((_i, p), (_j, q))| {
                        if distance_squared(&p.into(), &q.into()) <= cutoff_squared {
                        } else {
                        }
                    });
                })
            },
        );

        let pointcloud = generate_points_random(
            size,
            [vol_edges; 3],
            [0.0, 0.0, 0.0],
            Some(3079380797442975920),
        );

        // FIXME: this does not scale linearly because:
        // FIXME: 1. we generate new random point cloud (so FlatIndex does change massively)
        // FIXME: 2. cutoff remains unchanged
        // FIXME: Therefore it's increasingly unlikely that rebuild_mut() can spare some extra work
        // FIXME: In fact, for size=10_000_000 runtime is comparable to CellGrid::new()/rebuild()
        // FIXME: for size=1_000_000 this is not the case
        // FIXME: maybe some more investigation is needed.
        // FIXME: actually 3s warmup phase does differ from ::new()... why?
        // FIXME: (estimating target time ~130s vs 60s for 100 samples)
        group.bench_with_input(BenchmarkId::new("::rebuild_mut()", size), &cg, |b, cg| {
            let mut cg = cg.clone();
            b.iter(|| {
                cg.rebuild_mut(pointcloud.iter().map(|p| p.coords), None);
            })
        });
    }

    group.finish();
}

criterion_group!(benches, bench_cellgrid_concentration);
criterion_main!(benches);
