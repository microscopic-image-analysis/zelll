#[allow(dead_code)]
pub mod atom;
#[allow(dead_code)]
pub mod io;
#[allow(dead_code)]
pub mod sdf;
#[allow(dead_code)]
pub mod surface;

// `nuts-rs` wants f64 anyway
pub type Angstrom = f64;

pub mod utils {
    use super::Angstrom;
    use crate::atom::Atom;
    use nalgebra::{Point, SVector, distance};

    pub fn approx_geodesic_dist(
        x: Atom,
        y: Atom,
        x_normal: [Angstrom; 3],
        y_normal: [Angstrom; 3],
    ) -> Angstrom {
        let x: Point<Angstrom, 3> = x.coords.into();
        let y: Point<Angstrom, 3> = y.coords.into();

        let mut x_normal: SVector<Angstrom, 3> = x_normal.into();
        let mut y_normal: SVector<Angstrom, 3> = y_normal.into();

        // just to be safe, normalize again, this code is not performance-critical.
        x_normal.normalize_mut();
        y_normal.normalize_mut();

        // distance(&x, &y)
        distance(&x, &y) * (2.0 - x_normal.dot(&y_normal))
        // distance(&x, &y) * (2.0 - x_normal.dot(&y_normal)).sqrt()
        // distance(&x, &y) * (2.0 - distance(&x_normal.into(), &y_normal.into()))
    }
}
