use criterion::{
    black_box, criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion,
    PlotConfiguration, SamplingMode,
};
use nalgebra::{distance_squared, Point3, Vector3};
use zelll::cellgrid::{generate_points_random, CellGrid};

pub fn bench_cellgrid_random(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut g_new = c.benchmark_group("CellGrid::new()");
    g_new
        .sampling_mode(SamplingMode::Flat)
        .plot_config(plot_config.clone());

    for size in (0..=6).map(|exp| 10usize.pow(exp)) {
        let pointcloud = generate_points_random(size, 100.0, [0.0, 0.0, 0.0]);

        g_new.bench_with_input(
            BenchmarkId::from_parameter(size),
            &pointcloud,
            |b, pointcloud| {
                b.iter(|| {
                    CellGrid::new(pointcloud.iter(), 1.0);
                })
            },
        );
    }
    g_new.finish();

    let mut g_iter = c.benchmark_group("CellGrid iteration");
    g_iter
        .sampling_mode(SamplingMode::Flat)
        .plot_config(plot_config);

    //TODO: Would expect ~350ms for 10e6 points but I'm getting 1s consistently...
    //TODO: compare iteration against Itertools::tuple_combinations() (as naive N^2 version)
    //TODO: Also: count interactions, benchmark other iteration variants
    //TODO: also: benchmark with constant density instead of constant volume
    for size in (0..=6).map(|exp| 10usize.pow(exp)) {
        let pointcloud = generate_points_random(size, 100.0, [0.0, 0.0, 0.0]);
        let cg = CellGrid::new(pointcloud.iter(), 1.0);

        g_iter.bench_with_input(BenchmarkId::from_parameter(size), &cg, |b, cg| {
            b.iter(|| {
                #[cfg(not(feature = "rayon"))]
                cg.for_each_point_pair(|_, _| {
                    //std::hint::black_box(distance_squared(&pointcloud[i], &pointcloud[j]));
                });
                #[cfg(feature = "rayon")]
                cg.par_for_each_point_pair(|_, _| {});
            })
        });
    }
}

criterion_group!(benches, bench_cellgrid_random);
criterion_main!(benches);

