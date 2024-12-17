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

pub fn bench_cellgrid_random(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("CellGrid");
    group
        .sampling_mode(SamplingMode::Flat)
        .plot_config(plot_config.clone());

    group.bench_with_input(BenchmarkId::new("Gonnet2007", 1000), &1000, |b, size| {
        let pointcloud = generate_points_random(*size, [3.166; 3], [0.0, 0.0, 0.0], None);
        let cg = CellGrid::new(pointcloud.iter().map(|p| p.coords.as_ref()), 1.0);
        b.iter(|| {
            cg.for_each_point_pair(|_, _| {});
        })
    });

    #[cfg(feature = "rayon")]
    group.bench_with_input(
        BenchmarkId::new("Gonnet2007_par", 1000),
        &1000,
        |b, size| {
            let pointcloud = generate_points_random(*size, [3.166; 3], [0.0, 0.0, 0.0], None);
            let cg = CellGrid::new(pointcloud.iter().map(|p| p.coords.as_ref()), 1.0);
            b.iter(|| {
                cg.par_for_each_point_pair(|_, _| {});
            })
        },
    );

    for size in (0..=5).map(|exp| 10usize.pow(exp)) {
        let pointcloud = generate_points_random(size, [100.0; 3], [0.0, 0.0, 0.0], None);

        group.bench_with_input(
            BenchmarkId::new("::new()", size),
            &pointcloud,
            |b, pointcloud| {
                b.iter(|| {
                    CellGrid::new(pointcloud.iter().map(|p| p.coords.as_ref()), 10.0);
                })
            },
        );
    }

    //TODO: Would expect ~350ms for 10e6 points but I'm getting 1s consistently...
    //TODO: compare iteration against Itertools::tuple_combinations() (as naive N^2 version)
    //TODO: Also: count interactions, benchmark other iteration variants
    //TODO: also: benchmark with constant density instead of constant volume
    for size in (0..=5).map(|exp| 10usize.pow(exp)) {
        let pointcloud = generate_points_random(size, [100.0; 3], [0.0, 0.0, 0.0], None);
        let cg = CellGrid::new(pointcloud.iter().map(|p| p.coords.as_ref()), 10.0);

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

        let pointcloud = generate_points_random(size, [100.0; 3], [0.0, 0.0, 0.0], None);

        group.bench_with_input(BenchmarkId::new("::rebuild_mut()", size), &cg, |b, cg| {
            let mut cg = cg.clone();
            b.iter(|| {
                cg.rebuild_mut(pointcloud.iter().map(|p| p.coords.as_ref()), None);
            })
        });
    }

    group.finish();
}

pub fn bench_cellgrid_concentration(c: &mut Criterion) {
    let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
    let mut group = c.benchmark_group("CellGrid");
    group
        .sampling_mode(SamplingMode::Flat)
        .plot_config(plot_config.clone());

    group.bench_with_input(BenchmarkId::new("Gonnet2007", 1000), &1000, |b, size| {
        let pointcloud = generate_points_random(*size, [3.166; 3], [0.0, 0.0, 0.0], None);
        let cg = CellGrid::new(pointcloud.iter().map(|p| p.coords.as_ref()), 1.0);
        b.iter(|| {
            //cg.for_each_point_pair(|_, _| {});
            cg.filter_point_pairs(
                |_, _| {},
                |i, j| {
                    /*UniformNorm.metric_distance(&pointcloud[i].coords, &pointcloud[j].coords) <= 1.0
                    && */
                    distance_squared(&pointcloud[i], &pointcloud[j]) <= 1.0
                },
            );
        })
    });

    #[cfg(feature = "rayon")]
    group.bench_with_input(
        BenchmarkId::new("Gonnet2007_par", 1000),
        &1000,
        |b, size| {
            let pointcloud = generate_points_random(*size, [3.166; 3], [0.0, 0.0, 0.0], None);
            let cg = CellGrid::new(pointcloud.iter().map(|p| p.coords.as_ref()), 1.0);
            b.iter(|| {
                //cg.par_for_each_point_pair(|_, _| {});
                cg.par_filter_point_pairs(
                    |_, _| {},
                    |i, j| {
                        /*UniformNorm.metric_distance(&pointcloud[i].coords, &pointcloud[j].coords)
                        <= 1.0 &&*/
                        distance_squared(&pointcloud[i], &pointcloud[j]) <= 1.0
                    },
                );
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
                    CellGrid::new(pointcloud.iter().map(|p| p.coords.as_ref()), cutoff);
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
        let cg = CellGrid::new(pointcloud.iter().map(|p| p.coords.as_ref()), cutoff);

        group.bench_with_input(BenchmarkId::new("::point_pairs()", size), &cg, |b, cg| {
            let cutoff_squared = cutoff.powi(2);
            b.iter(|| {
                //cg.for_each_point_pair(|_, _| {});
                /*cg.filter_point_pairs(
                    |_, _| {},
                    |i, j| {
                        /*UniformNorm
                        .metric_distance(&pointcloud[i].coords, &pointcloud[j].coords)
                        <= cutoff && */
                        distance_squared(&pointcloud[i], &pointcloud[j]) <= cutoff_squared
                    },
                );*/
                // cg.point_pairs()
                //     .filter(|&(i, j)| {
                //         distance_squared(&pointcloud[i], &pointcloud[j]) <= cutoff_squared
                //     })
                cg.point_pairs2()
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
                    //cg.par_for_each_point_pair(|_, _| {});
                    /*
                    cg.par_filter_point_pairs(
                        |_, _| {},
                        |i, j| {
                            /*UniformNorm
                            .metric_distance(&pointcloud[i].coords, &pointcloud[j].coords)
                            <= cutoff && */
                            distance_squared(&pointcloud[i], &pointcloud[j]) <= cutoff_squared
                        },
                    );
                    */
                    /*
                    cg.par_point_pairs().filter(|&(i, j)| {
                        distance_squared(&pointcloud[i], &pointcloud[j]) <= cutoff_squared
                    })
                    .for_each(|_| {});
                    */
                    // cg.par_point_pairs().for_each(|(i, j)| {
                    //     if distance_squared(&pointcloud[i], &pointcloud[j]) <= cutoff_squared {
                    //     } else {
                    //     }
                    // });
                    cg.par_point_pairs2().for_each(|((_i, p), (_j, q))| {
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

        group.bench_with_input(BenchmarkId::new("::rebuild_mut()", size), &cg, |b, cg| {
            let mut cg = cg.clone();
            b.iter(|| {
                cg.rebuild_mut(pointcloud.iter().map(|p| p.coords.as_ref()), None);
            })
        });
    }

    group.finish();
}

criterion_group!(benches, bench_cellgrid_concentration);
criterion_main!(benches);
