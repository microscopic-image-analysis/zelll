//! # Usage
//!
//! ```sh
//! cargo run --release --example lmp-data -- <n> <seed> > atomsinabox.txt
//! ```
use nalgebra::{Point, Point3, Vector3};
use rand::distributions::Standard;
use rand::prelude::*;

type PointCloud<const N: usize> = Vec<Point<f64, N>>;
/// Generate a uniformly random 3D point cloud of size `n` in a cuboid of edge lengths `vol` centered around `origin`.
fn generate_points_random(
    n: usize,
    vol: [f64; 3],
    origin: [f64; 3],
    seed: Option<u64>,
) -> PointCloud<3> {
    // with fixed seed for reproducability
    let mut rng = StdRng::seed_from_u64(seed.unwrap_or(3079380797442975911));

    std::iter::repeat_with(|| {
        Point3::<f64>::from(
            (Vector3::from_iterator((&mut rng).sample_iter(Standard))
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
        .unwrap_or(1000);
    let seed = args
        .next()
        .and_then(|arg| arg.parse::<u64>().ok())
        .unwrap_or(3079380797442975911);

    let origin = [0.0; 3];

    let cutoff: f64 = 10.0;
    let conc = 10.0 / cutoff.powi(3); //i.e. 10mol per cutoff^3 volume units
    let a = 3.0 * cutoff;
    let b = 3.0 * cutoff;
    let c = (size as f64 / conc) / a / b;

    let pointcloud = generate_points_random(size, [a, b, c], origin, Some(seed));

    println!("# {size} random atom positions taken from zelll benchmarks:");
    println!(
        "# generate_points_random({}, {:?}, {:?}, Some({}));",
        size,
        [a, b, c],
        origin,
        seed,
    );

    println!("{size} atoms");
    println!("1 atom types");
    println!("-{} {} xlo xhi", 0.5 * a, 0.5 * a);
    println!("-{} {} ylo yhi", 0.5 * b, 0.5 * b);
    println!("-{} {} zlo zhi", 0.5 * c, 0.5 * c);
    println!("");
    println!("Atoms # atomic");

    for (i, atom) in pointcloud.iter().enumerate() {
        println!("{} 1 {} {} {}", i + 1, atom.x, atom.y, atom.z);
    }

    println!(
        "{} 1 {} {} {} # additional atom to make lammps read_data happy: https://docs.lammps.org/Errors_details.html#err0016",
        size + 1,
        origin[0],
        origin[1],
        origin[2],
    );
    println!("");
}
