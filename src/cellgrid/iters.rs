//TODO iterate over all neighboured cells (full/half space), pairs of particles
//TODO: perhaps move parallel iteration into separate submodule
use super::CellGrid;
use core::iter::FusedIterator;
use core::slice::Iter;
use itertools::Itertools;
use nalgebra::Point;
use num_traits::{AsPrimitive, ConstOne, Float, NumAssignOps};
#[cfg(feature = "rayon")]
use rayon::prelude::ParallelIterator;

/// `GridCell` represents a non-empty (by construction) cell of a `CellGrid`.
#[derive(Debug, Clone, Copy)]
pub struct GridCell<'g, const N: usize = 3, F: Float = f64>
where
    F: NumAssignOps + ConstOne + AsPrimitive<i32> + std::fmt::Debug,
{
    //TODO: maybe provide proper accessors to these fields for neighbors.rs to use?
    //TODO: is there a better way than having a reference to the containing CellGrid?
    pub(crate) grid: &'g CellGrid<N, F>,
    pub(crate) index: i32,
}

impl<'g, const N: usize, F> GridCell<'g, N, F>
where
    F: Float + NumAssignOps + ConstOne + AsPrimitive<i32> + Send + Sync + std::fmt::Debug,
{
    /// Return the (flat) index of this (non-empty) `GridCell`.
    pub(crate) fn index(&self) -> i32 {
        self.index
    }

    // TODO: should probably rather impl IntoIterator to match consuming/copy behaviour of neighbors()/point_pairs()?
    pub fn iter(&self) -> Iter<'g, (usize, Point<F, N>)> {
        self.grid
            .cell_lists
            .cell_slice(
                &self
                    .grid
                    .cells
                    .get(&self.index)
                    .expect("This GridCell should be contained in the CellGrid but it is not."),
            )
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

    /// Return an iterator over all (currently half-space) non-empty neighboring cells.
    //TODO: currently only half-space and aperiodic boundaries
    //TODO: handle half-/full-space  and (a-)periodic boundary conditions
    //TODO: document that we're relying on GridCell impl'ing Copy here (so we can safely consume `self`)
    pub fn neighbors(self) -> impl FusedIterator<Item = GridCell<'g, N, F>> + Clone {
        self.grid
            .index
            .neighbor_indices
            .iter()
            .filter_map(move |rel| {
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
    #[inline]
    fn intra_cell_pairs(self) -> impl FusedIterator<Item = (usize, usize)> + Clone + 'g {
        // this is equivalent to
        // self.iter().copied().tuple_combinations::<(usize, usize)>()
        // but faster for our specific case (pairs from slice of `Copy` values)
        self.iter()
            .copied()
            .enumerate()
            .flat_map(move |(n, i)| self.iter().copied().skip(n + 1).map(move |j| (i, j)))
            .map(|((i, _p), (j, _q))| (i, j))
    }

    /// Iterate over all unique pairs of points in this `GridCell` with points of the neighboring cells.
    #[inline]
    fn inter_cell_pairs(self) -> impl FusedIterator<Item = (usize, usize)> + Clone + 'g {
        self.iter()
            .copied()
            .cartesian_product(self.neighbors().flat_map(|cell| cell.iter().copied()))
            .map(|((i, _p), (j, _q))| (i, j))
    }

    #[inline]
    fn intra_cell_pairs2(
        self,
    ) -> impl FusedIterator<Item = ((usize, Point<F, N>), (usize, Point<F, N>))> + Clone + 'g {
        self.iter()
            .copied()
            .enumerate()
            .flat_map(move |(n, i)| self.iter().copied().skip(n + 1).map(move |j| (i, j)))
    }

    /// Iterate over all unique pairs of points in this `GridCell` with points of the neighboring cells.
    #[inline]
    fn inter_cell_pairs2(
        self,
    ) -> impl FusedIterator<Item = ((usize, Point<F, N>), (usize, Point<F, N>))> + Clone + 'g {
        self.iter()
            .copied()
            .cartesian_product(self.neighbors().flat_map(|cell| cell.iter().copied()))
    }

    pub fn point_pairs2(
        self,
    ) -> impl FusedIterator<Item = ((usize, [F; N]), (usize, [F; N]))> + Clone + Send + Sync + 'g
    {
        self.intra_cell_pairs2()
            .chain(self.inter_cell_pairs2())
            .map(|((i, p), (j, q))| ((i, p.into()), (j, q.into())))
    }

    /// Iterate over all "relevant" pairs of points within in the neighborhood of this `GridCell`.
    //TODO: explain what "relevant" means here.
    //TODO: handle full-space as well
    //TODO: document that we're relying on GridCell impl'ing Copy here (so we can safely consume `self`)
    pub fn point_pairs(
        self,
    ) -> impl FusedIterator<Item = (usize, usize)> + Clone + Send + Sync + 'g {
        self.intra_cell_pairs().chain(self.inter_cell_pairs())
    }
}

impl<const N: usize, F> CellGrid<N, F>
where
    F: Float + NumAssignOps + ConstOne + AsPrimitive<i32> + Send + Sync + std::fmt::Debug,
{
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
    pub fn iter(&self) -> impl FusedIterator<Item = GridCell<N, F>> + Clone {
        // note that ::keys() does not keep a stable iteration order!
        self.cells.keys().map(|&index| GridCell {
            grid: self,
            index: index,
        })
    }

    #[cfg(feature = "rayon")]
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn par_iter(&self) -> impl ParallelIterator<Item = GridCell<N, F>> {
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
