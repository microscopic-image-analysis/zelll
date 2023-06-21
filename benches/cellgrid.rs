use criterion::{
    black_box, criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion,
    PlotConfiguration, SamplingMode,
};
use nalgebra::{distance_squared, Point3, Vector3};
use zelll::cellgrid::{generate_points_random, CellGrid};

pub fn bench_cellgrid_random(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("CellGrid");
    group
        .sampling_mode(SamplingMode::Flat)
        .plot_config(plot_config.clone());

    for size in (0..=6).map(|exp| 10usize.pow(exp)) {
        let pointcloud = generate_points_random(size, 100.0, [0.0, 0.0, 0.0]);

        group.bench_with_input(
            BenchmarkId::new("::new()", size),
            &pointcloud,
            |b, pointcloud| {
                b.iter(|| {
                    CellGrid::new(pointcloud.iter(), 1.0);
                })
            },
        );
    }

    //TODO: Would expect ~350ms for 10e6 points but I'm getting 1s consistently...
    //TODO: compare iteration against Itertools::tuple_combinations() (as naive N^2 version)
    //TODO: Also: count interactions, benchmark other iteration variants
    //TODO: also: benchmark with constant density instead of constant volume
    for size in (0..=6).map(|exp| 10usize.pow(exp)) {
        let pointcloud = generate_points_random(size, 100.0, [0.0, 0.0, 0.0]);
        let cg = CellGrid::new(pointcloud.iter(), 1.0);

        group.bench_with_input(
            BenchmarkId::new("::for_each_point_pair()", size),
            &cg,
            |b, cg| {
                b.iter(|| {
                    cg.for_each_point_pair(|_, _| {
                        //std::hint::black_box(distance_squared(&pointcloud[i], &pointcloud[j]));
                    });
                })
            },
        );

        #[cfg(feature = "rayon")]
        group.bench_with_input(
            BenchmarkId::new("::par_for_each_point_pair()", size),
            &cg,
            |b, cg| {
                b.iter(|| {
                    cg.par_for_each_point_pair(|_, _| {});
                })
            },
        );

        let pointcloud = generate_points_random(size, 100.0, [0.0, 0.0, 0.0]);

        group.bench_with_input(BenchmarkId::new("::rebuild_mut()", size), &cg, |b, cg| {
            let mut cg = cg.clone();
            b.iter(|| {
                cg.rebuild_mut(pointcloud.iter(), None);
            })
        });
    }
    group.finish();
}

criterion_group!(benches, bench_cellgrid_random);
criterion_main!(benches);

