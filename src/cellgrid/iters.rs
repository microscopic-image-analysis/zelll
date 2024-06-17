//TODO iterate over all neighboured cells (full/half space), pairs of particles
//TODO: perhaps move parallel iteration into separate submodule
use super::{CellGrid, CellNeighbors};
use core::iter::FusedIterator;
use itertools::Itertools;
#[cfg(feature = "rayon")]
use rayon::prelude::ParallelIterator;

/// `GridCell` represents a non-empty (by construction) cell of a `CellGrid`.
#[derive(Debug, Clone, Copy)]
pub struct GridCell<'g, const N: usize> {
    //TODO: maybe provide proper accessors to these fields for neighbors.rs to use?
    pub(crate) grid: &'g CellGrid<N>,
    pub(crate) index: i32,
}

impl<'g, const N: usize> GridCell<'g, N> {
    /// Return the (flat) index of this (non-empty) `GridCell`.
    pub(crate) fn index(&self) -> i32 {
        self.index
    }

    pub fn iter(&self) -> std::slice::Iter<'g, usize> {
        self.grid
            .cell_lists
            .cell_slice(&self.grid.cells[&self.index])
            .iter()
    }
    /*
    /// Check whether this `GridCell` is on the boundary of the [`CellGrid`].
    //TODO: I don't think I need it but let's keep it anyway
    pub fn on_boundary(&self) -> bool {
        let idx = self.index();
        let shape = self.grid.index.grid_info.shape;

        for (i, dim) in shape.iter().enumerate() {
            if idx[i] == 0 || idx[i] + 1 == *dim {
                return true;
            }
        }
        false
    }*/

    /// Return [`CellNeighbors`], an iterator over all (currently half-space) non-empty neighboring cells.
    #[deprecated(
        note = "cell_neighbors() returns CellNeighbors<N> which will be phased out soon. Use neighbors() instead."
    )]
    pub fn cell_neighbors(&self) -> CellNeighbors<N> {
        CellNeighbors::half_space(self)
    }

    /// Return an iterator over all (currently half-space) non-empty neighboring cells.
    //TODO: currently only half-space and aperiodic boundaries
    //TODO: handle half-/full-space  and (a-)periodic boundary conditions
    pub fn neighbors(&self) -> impl FusedIterator<Item = GridCell<'g, N>> + Clone + '_ {
        //TODO: could pre-compute the neighbor indices and store in self.grid.cells?
        //todo: OR: this should be able to use SIMD efficiently, right?
        //TODO: just chunk neighbor_indices and add onto it (however, it requires another allocation then?)
        //TODO: or do the iteration in chunks
        //TODO: OR: in principle could just increment neighbor_indices +1 to get abs. neighbours for the next cell
        //TODO: However, makes parallel iteration ugly and I'd have to make sure to reset relative_indices afterwards
        //TODO: also, doesn't really matter if I do +=1 or +self.index()
        self.grid
            .index
            .neighbor_indices
            .iter()
            .filter_map(|rel| {
                let index = rel + self.index();
                //TODO: I mean I could also store the slice since I'm already accessing its metadata?
                //TODO: would save me one lookup into the hashmap self.grid.cells in ::iter()
                //TODO: tested it, I gain ~2pairs/kcycle from it (32->34) might think about it again
                //TODO: ::iter() is nicer this way, other parts less so
                self.grid.cells.get(&index).map(|_| GridCell {
                    grid: self.grid,
                    index,
                })
            })
    }

    /// Iterate over all unique pairs of points in this `GridCell`.
    fn intra_cell_pairs(
        &self,
    ) -> impl FusedIterator<Item = (usize, usize)> + Clone + '_ {
        self.iter().copied().tuple_combinations::<(usize, usize)>()
    }

    /// Iterate over all unique pairs of points in this `GridCell` with points of the neighboring cells.
    fn inter_cell_pairs(
        &self,
    ) -> impl FusedIterator<Item = (usize, usize)> + Clone + '_ {
        //TODO: storing neighboring slices in a temporary allocation improves sequential iteration
        //TODO: but negatively impacts parallel iteration (perhaps due to multiple threads wanting to allocate concurrently?)
        /*let others: Vec<_> = self
            .neighbors()
            .flat_map(|cell| cell.iter().copied())
            .collect();
        self.iter().copied().cartesian_product(others)*/

        //TODO: might help instead of having a closure in the iterator
        //TODO: if I wanted introduce a type alias?
        /*
        fn copied_cell_iter<const N: usize>(cell: GridCell<N>) -> std::iter::Copied<std::slice::Iter<'_, usize>> {
            cell.iter().copied()
        }
        */
        self.iter()
            .copied()
            .cartesian_product(self.neighbors().flat_map(|cell| cell.iter().copied()))
            //.cartesian_product(self.neighbors().flat_map(copied_cell_iter))

        //TODO: itertools >= "0.12" cartesian_product() doesn't optimize well?
        //TODO: possibly because of nested Option in https://github.com/rust-itertools/itertools/pull/800
        //TODO: which is likely the correct thing to do but it's unfortunate to lose 10-15% of performance
        //TODO: For itertools >= "0.12", two nested flat_maps perform just as well as cartesian_product()
        /*self.iter().copied().flat_map(move |i| {
            self.neighbors()
                .flat_map(|cell| cell.iter().copied())
                .map(move |j| (i, j))
        })*/
    }

    /// Iterate over all "relevant" pairs of points within in the neighborhood of this `GridCell`.
    //TODO: explain what "relevant" means here.
    //TODO: handle full-space as well
    pub fn point_pairs(&self) -> impl FusedIterator<Item = (usize, usize)> + Clone + '_ {
        self.intra_cell_pairs().chain(self.inter_cell_pairs())
    }
}

impl<const N: usize> CellGrid<N> {
    /// Returns an iterator over all [`GridCell`]s in this `CellGrid`, excluding empty cells.
    /// A particular iteration order is not guaranteed.
    ///
    /// # Examples
    ///
    //TODO: this example should still work but it's nonsensical
    /// ```
    /// # use zelll::cellgrid::CellGrid;
    /// # use nalgebra::Point;
    /// # let points = [Point::from([0.0, 0.0, 0.0]), Point::from([1.0,2.0,0.0]), Point::from([0.0, 0.1, 0.2])];
    /// # let points = [[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
    /// # let cell_grid = CellGrid::new(points.iter(), 1.0);
    /// cell_grid.iter().flat_map(|cell| cell.iter()).count();
    /// ```
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn iter(&self) -> impl FusedIterator<Item = GridCell<N>> + Clone {
        self.cells
            .keys()
            .map(|&index| GridCell { grid: self, index })
    }

    #[cfg(feature = "rayon")]
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn par_iter(&self) -> impl ParallelIterator<Item = GridCell<N>> {
        self.cells
            .par_keys()
            .map(|&index| GridCell { grid: self, index })
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
        let cell_grid: CellGrid<3> = CellGrid::new(points.iter(), 1.0);

        assert_eq!(cell_grid.iter().count(), 14, "testing iter()");

        #[cfg(feature = "rayon")]
        assert_eq!(cell_grid.par_iter().count(), 14, "testing par_iter()");
    }

    #[test]
    fn test_gridcell_iter() {
        // Using 0-origin to avoid floating point errors
        let points = generate_points([3, 3, 3], 1.0, [0.0, 0.0, 0.0]);
        let cell_grid: CellGrid<3> = CellGrid::new(points.iter(), 1.0);

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
        let cell_grid: CellGrid<3> = CellGrid::new(points.iter(), 1.0);

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

