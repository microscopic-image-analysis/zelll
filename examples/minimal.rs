use hashbrown::{HashMap, HashSet};
//use nohash_hasher::BuildNoHashHasher;
use nalgebra::{Point, Point3, Vector3};
use rand::distributions::Standard;
use rand::prelude::*;
#[cfg(feature = "rayon")]
use rayon::prelude::ParallelIterator;
use std::hint::black_box;
use std::iter::FromIterator;
use zelll::cellgrid::CellGrid;

type PointCloud<const N: usize> = Vec<Point<f64, N>>;
/// Generate a uniformly random 3D point cloud of size `n` in a cuboid of edge lengths `vol` centered around `origin`.
fn generate_points_random(n: usize, vol: [f64; 3], origin: [f64; 3]) -> PointCloud<3> {
    std::iter::repeat_with(|| {
        Point3::<f64>::from(
            (Vector3::from_iterator(thread_rng().sample_iter(Standard))
                - Vector3::new(0.5, 0.5, 0.5)
                + Vector3::from(origin))
            .component_mul(&Vector3::from(vol)),
        )
    })
    .take(n)
    .collect()
}

fn main() {
    for size in (2..=6).map(|exp| 10usize.pow(exp)) {
        let cutoff: f64 = 10.0;
        let conc = 10.0 / cutoff.powi(3); //i.e. 100mol per 10^3 volume units
        let a = 3.0 * cutoff;
        let b = 3.0 * cutoff;
        let c = (size as f64 / conc) / a / b;
        let _vol_edges = (size as f64 / conc).cbrt();
        let pointcloud = generate_points_random(size, [a, b, c], [0.0, 0.0, 0.0]);

        let cg = CellGrid::new(pointcloud.iter().map(|p| p.coords.as_ref()), cutoff);
        println!("{:?}", cg.shape());
        let _cutoff_squared = cutoff.powi(2);

        //cg.for_each_point_pair(|_, _| black_box(()));
        let mut count: usize = 0;
        #[cfg(not(feature = "rayon"))]
        /*cg.filter_point_pairs(
            |_, _| {
                count += 1;
            },
            |_i, _j| true, //distance_squared(&pointcloud[i], &pointcloud[j]) <= cutoff_squared,
        );*/
        let count = cg.point_pairs().count();
        #[cfg(feature = "rayon")]
        cg.par_filter_point_pairs(
            |_, _| {
                //count += 1;
                black_box(());
            },
            |_i, _j| true, //distance_squared(&pointcloud[i], &pointcloud[j]) <= cutoff_squared,
        );
        //let count = cg.par_point_pairs().count();
        println!("{}", count);
    }

    /*let vals = 0..100usize;
    let vals = vec![
        443, 262, 52, 843, 917, 777, 323, 801, 73, 961, 780, 212, 220, 852, 619, 649, 297, 689,
        588, 617, 146, 105, 560, 481, 108, 442, 796, 299, 261, 257, 740, 541, 337, 435, 47, 197,
        545, 770, 564, 502, 753, 203, 524, 982, 29, 946, 600, 253, 558, 907, 886, 394, 491, 951,
        272, 151, 31, 904, 206, 863, 464, 772, 799, 610, 125, 501, 938, 399, 701, 44, 335, 543,
        631, 940, 759, 252, 351, 413, 978, 602, 972, 254, 193, 949, 553, 80, 324, 711, 851, 969,
        721, 308, 540, 576, 923, 875, 751, 14, 412, 27,
    ];
    let hs: HashSet<usize, BuildNoHashHasher<usize>> = HashSet::from_iter(vals);
    let collected: Vec<_> = hs.iter().collect();
    println!("{collected:?}");
    let collected: Vec<_> = hs
        .iter()
        .map(|i| (i, i % hs.raw_table().buckets()))
        .collect();
    println!("{collected:?}");
    println!(
        "{} {} {}",
        hs.len(),
        hs.capacity(),
        hs.raw_table().buckets()
    );*/
}
