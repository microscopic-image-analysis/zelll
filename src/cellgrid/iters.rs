//TODO iterate over all neighboured cells (full/half space), pairs of particles
//TODO: perhaps move parallel iteration into separate submodule
use super::{CellGrid, CellNeighbors};
//TODO: hide this behind a featureflag?
use itertools::Itertools;
#[cfg(feature = "rayon")]
use ndarray::parallel::prelude::*;

#[derive(Debug, Clone, Copy)]
pub struct GridCell<'g, const N: usize> {
    //TODO: maybe provide proper accessors to these fields for neighbors.rs to use?
    pub(crate) grid: &'g CellGrid<N>,
    //TODO: I guess I could just store Option<usize> directly as well
    pub(crate) head: &'g Option<usize>,
}

impl<'g, const N: usize> GridCell<'g, N> {
    pub fn is_empty(&self) -> bool {
        self.head.is_none()
    }

    /// Return the (multi-)index of this (non-empty) `GridCell`.
    /// Returns `None` if the cell is empty.
    //TODO: However, the public API should not provide a way to address empty cells
    pub(crate) fn index(&self) -> Option<[usize; N]> {
        let idx = (*self.head)?;
        Some(self.grid.index.index[idx])
    }

    pub fn iter(&self) -> GridCellIterator<'g, N> {
        GridCellIterator {
            grid: self.grid,
            state: self.head,
        }
    }

    /// Check whether this `GridCell` is on the boundary of the [`CellGrid`].
    /// Returns None if the cell is empty.
    //TODO: Again, empty cells shouldn't be accessible via the public API
    //TODO: I don't think I need it but let's keep it anyway
    //TODO: Semantically, Result<> would make more sense
    pub fn on_boundary(&self) -> Option<bool> {
        let idx = self.index()?;
        let shape = self.grid.index.grid_info.shape;

        for (i, dim) in shape.iter().enumerate() {
            if idx[i] == 0 || idx[i] + 1 == *dim {
                return Some(true);
            }
        }
        Some(false)
    }

    /// Return [`CellNeighbors`], an iterator over all (currently half-space) non-empty neighboring cells.
    //TODO: currently only half-space and aperiodic boundaries
    //TODO: handle half-/full-space  and (a-)periodic boundary conditions
    pub fn neighbors(&self) -> CellNeighbors<N> {
        CellNeighbors::half_space(self)
    }

    /// Iterate over all unique pairs of points in this `GridCell`.
    fn intra_cell_pairs(&self) -> impl Iterator<Item = (usize, usize)> + 'g {
        self.iter().tuple_combinations::<(usize, usize)>()
    }

    /// Iterate over all unique pairs of points in this `GridCell` with points of the neighboring cells.
    fn inter_cell_pairs(&'g self) -> impl Iterator<Item = (usize, usize)> + 'g {
        self.iter()
            .cartesian_product(self.neighbors().flat_map(|neighbor| neighbor.iter()))
    }

    /// Iterate over all "relevant" pairs of points within in the neighborhood of this `GridCell`.
    //TODO: explain what "relevant" means here.
    //TODO: handle full-space as well
    pub fn point_pairs(&'g self) -> impl Iterator<Item = (usize, usize)> + 'g {
        self.intra_cell_pairs().chain(self.inter_cell_pairs())
    }
}
/// Iterates over all points (or rather their indices) in the [`GridCell`] this `GridCellIterator` was created from.
//TODO: impl FusedIterator
#[derive(Debug, Clone, Copy)]
#[must_use = "iterators are lazy and do nothing unless consumed"]
pub struct GridCellIterator<'g, const N: usize> {
    grid: &'g CellGrid<N>,
    //TODO: I guess I could just store Option<usize> directly as well
    state: &'g Option<usize>,
}

impl<const N: usize> Iterator for GridCellIterator<'_, N> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        //TODO: use Option.map()?
        if let Some(index) = self.state {
            self.state = &self.grid.cell_lists[*index];
            Some(*index)
        } else {
            None
        }
    }
}

impl<const N: usize> CellGrid<N> {
    /// Returns an iterator over all [`GridCell`]s in this `CellGrid`, excluding empty cells.
    ///
    /// # Examples
    ///
    /// ```
    /// # use zelll::cellgrid::CellGrid;
    /// # let points = [[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
    /// # let cell_grid = CellGrid::new(&points, 1.0);
    /// cell_grid.iter().filter(|cell| !cell.is_empty());
    /// ```
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn iter(&self) -> impl Iterator<Item = GridCell<N>> {
        self.cells
            .iter()
            // It seems a bit weird but I'm just moving a reference to self (if I'm not mistaken).
            .filter_map(move |head| {
                if head.is_some() {
                    Some(GridCell { grid: self, head })
                } else {
                    None
                }
            })
    }

    #[cfg(feature = "rayon")]
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn par_iter(&self) -> impl ParallelIterator<Item = GridCell<N>> {
        self.cells
            .par_iter()
            // It seems a bit weird but I'm just moving a reference to self (if I'm not mistaken).
            .filter_map(move |head| {
                if head.is_some() {
                    Some(GridCell { grid: self, head })
                } else {
                    None
                }
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cellgrid::util::generate_points;

    #[test]
    fn test_cellgrid_iter() {
        // Using 0-origin to avoid floating point errors
        let points = generate_points([3, 3, 3], 1.0, [0.0, 0.0, 0.0]);
        let cell_grid: CellGrid<3> = CellGrid::new(&points, 1.0);

        assert_eq!(cell_grid.iter().count(), 14, "testing iter()");

        #[cfg(feature = "rayon")]
        assert_eq!(cell_grid.par_iter().count(), 14, "testing par_iter()");
    }

    #[test]
    fn test_gridcell_iter() {
        // Using 0-origin to avoid floating point errors
        let points = generate_points([3, 3, 3], 1.0, [0.0, 0.0, 0.0]);
        let cell_grid: CellGrid<3> = CellGrid::new(&points, 1.0);

        assert_eq!(
            cell_grid.iter().flat_map(|cell| cell.iter()).count(),
            points.len(),
            "testing iter()"
        );

        //TODO: test GridCell.index() and GridCell.on_boundary()

        #[cfg(feature = "rayon")]
        assert_eq!(
            cell_grid
                .par_iter()
                .flat_map_iter(|cell| cell.iter())
                .count(),
            points.len(),
            "testing par_iter()"
        );
    }

    #[test]
    fn test_neighborcell_pointpairs() {
        // Using 0-origin to avoid floating point errors
        let points = generate_points([3, 3, 3], 1.0, [0.0, 0.0, 0.0]);
        let cell_grid: CellGrid<3> = CellGrid::new(&points, 1.0);

        //TODO: test intra and inter cell pairs

        assert_eq!(
            cell_grid.iter().flat_map(|cell| cell.iter()).count(),
            points.len(),
            "testing iter()"
        );
    }
}

