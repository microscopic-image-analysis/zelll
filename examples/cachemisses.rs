//! # Usage
//!
//! ```sh
//! valgrind --tool=callgrind --cache-sim=yes --instr-atstart=no target/release/examples/cachemisses 100000 1
//! ```

use crabgrind::callgrind as valgrind;
use nalgebra::{Point, Point3, Vector3};
use rand::distributions::Standard;
use rand::prelude::*;
use zelll::CellGrid;

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
    let mut args = std::env::args();
    args.next();

    let size = args
        .next()
        .and_then(|arg| arg.parse::<usize>().ok())
        .unwrap_or(100);
    let repeat = args
        .next()
        .and_then(|arg| arg.parse::<usize>().ok())
        .unwrap_or(1);
    let sort = args
        .next()
        .and_then(|arg| arg.parse::<bool>().ok())
        .unwrap_or(false);
    let single_precision = args
        .next()
        .and_then(|arg| arg.parse::<bool>().ok())
        .unwrap_or(false);

    let cutoff: f64 = 10.0;
    let conc = 10.0 / cutoff.powi(3); //i.e. 100mol per 10^3 volume units
    let a = 3.0 * cutoff;
    let b = 3.0 * cutoff;
    let c = (size as f64 / conc) / a / b;
    let _vol_edges = (size as f64 / conc).cbrt();
    let mut pointcloud = generate_points_random(size, [a, b, c], [0.0, 0.0, 0.0]);

    if sort {
        pointcloud.sort_unstable_by(|p, q| p.z.partial_cmp(&q.z).unwrap());
    }

    if single_precision {
        let cutoff = cutoff as f32;
        let pointcloud: Vec<_> = pointcloud.iter().map(|p| p.cast::<f32>()).collect();

        valgrind::start_instrumentation();
        for _ in 0..repeat {
            let _cg = CellGrid::new(pointcloud.iter().map(|p| p.coords), cutoff);
        }
        valgrind::stop_instrumentation();
    } else {
        valgrind::start_instrumentation();
        for _ in 0..repeat {
            let _cg = CellGrid::new(pointcloud.iter().map(|p| p.coords), cutoff);
        }
        valgrind::stop_instrumentation();
    }
}
