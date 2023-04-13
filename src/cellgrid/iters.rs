//TODO iterate over all neighboured cells (full/half space), pairs of particles
//TODO: perhaps move parallel iteration into separate submodule
use super::{CellGrid, CellNeighbors};
use core::iter::FusedIterator;
//TODO: hide this behind a featureflag?
use itertools::Itertools;
#[cfg(feature = "rayon")]
use rayon::prelude::ParallelIterator;

#[derive(Debug, Clone, Copy)]
pub struct GridCell<'g, const N: usize> {
    //TODO: maybe provide proper accessors to these fields for neighbors.rs to use?
    pub(crate) grid: &'g CellGrid<N>,
    //TODO: don't really need Option<> here but it makes it easier to distinguish non-empty from empty cells
    //TODO: although I could do this at another place using HashMaps Entry API
    pub(crate) head: Option<usize>,
}

impl<'g, const N: usize> GridCell<'g, N> {
    //TODO: remove or emulate by storing Option<usize> in GridCell
    pub fn is_empty(&self) -> bool {
        self.head.is_none()
    }

    /// Return the (multi-)index of this (non-empty) `GridCell`.
    /// Returns `None` if the cell is empty.
    //TODO: However, the public API should not provide a way to address empty cells
    pub(crate) fn index(&self) -> Option<[usize; N]> {
        let idx = self.head?;
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
    fn intra_cell_pairs(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.iter().tuple_combinations::<(usize, usize)>()
    }

    /// Iterate over all unique pairs of points in this `GridCell` with points of the neighboring cells.
    fn inter_cell_pairs(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.iter()
            .cartesian_product(self.neighbors().flat_map(|neighbor| neighbor.iter()))
    }

    /// Iterate over all "relevant" pairs of points within in the neighborhood of this `GridCell`.
    //TODO: explain what "relevant" means here.
    //TODO: handle full-space as well
    pub fn point_pairs(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.intra_cell_pairs().chain(self.inter_cell_pairs())
    }
}
/// Iterates over all points (or rather their indices) in the [`GridCell`] this `GridCellIterator` was created from.
#[derive(Debug, Clone, Copy)]
#[must_use = "iterators are lazy and do nothing unless consumed"]
pub struct GridCellIterator<'g, const N: usize> {
    grid: &'g CellGrid<N>,
    //TODO: I guess I could just store Option<usize> directly as well
    state: Option<usize>,
}

impl<const N: usize> Iterator for GridCellIterator<'_, N> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        self.state.map(|index| {
            self.state = self.grid.cell_lists[index];
            index
        })
    }
}

impl<const N: usize> FusedIterator for GridCellIterator<'_, N> {}

impl<const N: usize> CellGrid<N> {
    /// Returns an iterator over all [`GridCell`]s in this `CellGrid`, excluding empty cells.
    /// A particular iteration order is not guaranteed.
    ///
    /// # Examples
    ///
    //TODO: this example should still work but it's nonsensical
    /// ```
    /// # use zelll::cellgrid::CellGrid;
    /// # let points = [[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
    /// # let cell_grid = CellGrid::new(&points, 1.0);
    /// cell_grid.iter().filter(|cell| !cell.is_empty());
    /// ```
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn iter(&self) -> impl Iterator<Item = GridCell<N>> {
        self.cells
            .values()
            // It seems a bit weird but I'm just moving a reference to self (if I'm not mistaken).
            .map(move |&head| GridCell {
                grid: self,
                head: Some(head),
            })
    }
    //TODO: parallel iteration is broken, now that we use std::collections::HashMap instead of ndarray::ArrayD...
    //TODO: However, if we instead switch to the crate hashbrown, we can still have that
    //TODO: (also gives us AHash with potentially better performance than SipHash?)
    //TODO: hashbrown does seem to spend less time hashing (using AHash)
    //TODO: but slightly increases cache misses...
    #[cfg(feature = "rayon")]
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn par_iter(&self) -> impl ParallelIterator<Item = GridCell<N>> {
        self.cells
            .par_values()
            // It seems a bit weird but I'm just moving a reference to self (if I'm not mistaken).
            .map(move |&head| GridCell {
                grid: self,
                head: Some(head),
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
        let points = generate_points([2, 2, 2], 1.0, [0.0, 0.0, 0.0]);
        let cell_grid: CellGrid<3> = CellGrid::new(&points, 1.0);

        assert_eq!(
            cell_grid
                .iter()
                .map(|cell| cell.intra_cell_pairs().count())
                .sum::<usize>(),
            4,
            "testing intra_cell_pairs()"
        );

        assert_eq!(
            cell_grid
                .iter()
                .map(|cell| cell.inter_cell_pairs().count())
                .sum::<usize>(),
            24,
            "testing inter_cell_pairs()"
        );
    }
}

