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
}

impl SmoothDistanceField {
    pub fn new(protein: &PointCloud, cutoff: Angstrom) -> Self {
        Self {
            inner: CellGrid::new(protein.points.iter().copied(), cutoff),
            surface_radius: 1.05,
        }
    }

    /// Sets a custom radius defining the targeted iso-surface for sampling.
    pub fn with_surface_radius(self, surface_radius: Angstrom) -> Self {
        Self {
            surface_radius,
            ..self
        }
    }
}
