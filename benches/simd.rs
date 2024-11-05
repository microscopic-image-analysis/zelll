// #![feature(iter_array_chunks)]

// use criterion::{
//     criterion_group, criterion_main, AxisScale, BenchmarkId, Criterion,
//     PlotConfiguration, SamplingMode,
// };
// use itertools13::Itertools;
// use nalgebra::{distance_squared, Point, Point3, Vector3};
// use rand::distributions::Standard;
// use rand::prelude::*;
// use zelll::cellgrid::CellGrid;

// type PointCloud<const N: usize> = Vec<Point<f64, N>>;
// /// Generate a uniformly random 3D point cloud of size `n` in a cuboid of edge lengths `vol` centered around `origin`.
// fn generate_points_random(
//     n: usize,
//     vol: [f64; 3],
//     origin: [f64; 3],
//     seed: Option<u64>,
// ) -> PointCloud<3> {
    // with fixed seed for reproducability
//     let mut rng = StdRng::seed_from_u64(seed.unwrap_or(3079380797442975911));

//     std::iter::repeat_with(|| {
//         Point3::<f64>::from(
//             (Vector3::from_iterator((&mut rng).sample_iter(Standard))
//                 - Vector3::new(0.5, 0.5, 0.5)
//                 + Vector3::from(origin))
//             .component_mul(&Vector3::from(vol)),
//         )
//     })
//     .take(n)
//     .collect()
// }

// pub fn bench_simd(c: &mut Criterion) {
//     let plot_config = PlotConfiguration::default().summary_scale(AxisScale::Logarithmic);
//     let mut group = c.benchmark_group("Vectorization");
//     group
//         .sampling_mode(SamplingMode::Flat)
//         .plot_config(plot_config.clone());


//     for size in (2..=5).map(|exp| 10usize.pow(exp)) {
//         let cutoff: f64 = 10.0;
//         let conc = 10.0 / cutoff.powi(3); //i.e. 10mol per cutoff^3 volume units
//         let a = 3.0 * cutoff;
//         let b = 3.0 * cutoff;
//         let c = (size as f64 / conc) / a / b;

//         let pointcloud = generate_points_random(size, [a, b, c], [0.0, 0.0, 0.0], None);
//         let cg = CellGrid::new(pointcloud.iter().map(|p| p.coords.as_ref()), cutoff);

//         group.bench_with_input(BenchmarkId::new("No vectorization", size), &cg, |b, cg| {
//             let cutoff_squared = cutoff.powi(2);
//             b.iter(|| {
                // cg.point_pairs()
                //     .filter(|&(i, j)| {
                //         distance_squared(&pointcloud[i], &pointcloud[j]) <= cutoff_squared
                //     })
                //     .for_each(|_| {});
//                 cg.point_pairs()
//                     .for_each(|(i, j)| {
//                         let _distancecheck = distance_squared(&pointcloud[i], &pointcloud[j]) <= cutoff_squared;
//                     });
//             })
//         });

//         group.bench_with_input(BenchmarkId::new("Auto-vectorization", size), &cg, |b, cg| {
//             let cutoff_squared = cutoff.powi(2);
//             b.iter(|| {
//                 let mut chunked = cg.point_pairs().array_chunks::<4>();
//                 chunked.by_ref().for_each(|chunk|{
//                     chunk.iter().for_each(|&(i, j)| {
//                         let _distancecheck = distance_squared(&pointcloud[i], &pointcloud[j]) <= cutoff_squared;
//                     });
//                 });

//                 if let Some(rem) = chunked.into_remainder() {
//                     rem.for_each(|(i, j)| {
//                         let _distancecheck = distance_squared(&pointcloud[i], &pointcloud[j]) <= cutoff_squared;
//                     });
//                 }
                // TODO: with chunks I have to compute distances first, and map &(i,j) to [(&(i,j),d); N]
                // which I then have to flatten() + filter() or put filtering step into potential computation
//             })
//         });

//         group.bench_with_input(BenchmarkId::new("Itertools::chunk_by()", size), &cg, |b, cg| {
//             let cutoff_squared = cutoff.powi(2);
//             b.iter(|| {
//                 cg.point_pairs()
//                     .chunk_by(|&(i, _)| i)
//                     .into_iter()
//                     .for_each(|(i, chunk)| {
//                         let pi = &pointcloud[i];
//                         chunk.for_each(|(_, j)| {
//                             let _distancecheck = distance_squared(pi, &pointcloud[j]) <= cutoff_squared;
//                         })
                        //let _distancecheck = distance_squared(&pointcloud[i], &pointcloud[j]) <= cutoff_squared;
//                     });
//             })
//         });

        // group.bench_with_input(BenchmarkId::new("Explicit vectorization", size), &cg, |b, cg| {
        //     let cutoff_squared = cutoff.powi(2);
        //     b.iter(|| {
        //         cg.point_pairs()
        //             .filter(|&(i, j)| {
        //                 distance_squared(&pointcloud[i], &pointcloud[j]) <= cutoff_squared
        //             })
        //             .for_each(|_| {});
        //     })
        // });
//     }

//     group.finish();
// }


// criterion_group!(benches, bench_simd);
// criterion_main!(benches);