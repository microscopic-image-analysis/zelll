#[allow(dead_code)]
#[cfg(feature = "enzyme")]
mod enzyme;
#[cfg(not(feature = "enzyme"))]
mod numdual;

#[cfg(feature = "enzyme")]
#[allow(unused_imports)]
pub use enzyme::*;
#[cfg(not(feature = "enzyme"))]
#[allow(unused_imports)]
pub use numdual::*;

use crate::Angstrom;
use crate::atom::Atom;
use crate::io::PointCloud;
use zelll::CellGrid;

#[derive(Clone, Debug)]
pub struct SmoothDistanceField {
    inner: CellGrid<Atom, 3, Angstrom>,
    pub(crate) surface_radius: Angstrom,
    pub(crate) k_force: Angstrom, // this wouldn't survive dimensional analysis though
}

impl SmoothDistanceField {
    pub fn new(protein: &PointCloud, cutoff: Angstrom) -> Self {
        Self {
            inner: CellGrid::new(protein.points.iter().copied(), cutoff),
            surface_radius: 1.05,
            k_force: 10.0,
        }
    }

    /// Sets a custom radius defining the targeted iso-surface for sampling.
    pub fn with_surface_radius(self, surface_radius: Angstrom) -> Self {
        Self {
            surface_radius,
            ..self
        }
    }

    /// Sets a custom force constant for the iso-surface potential.
    /// Note that the actual physical dimension is not length.
    pub fn with_k_force(self, k_force: Angstrom) -> Self {
        Self { k_force, ..self }
    }
}
