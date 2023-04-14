//TODO: - allow both &[Point<f64, N>] and impl Iterator<Item = Point<f64, N>> (or IntoIterator)?
//TODO: - possible approach: require impl Iterator and then .take(usize::MAX) to ensure we have only finite point clouds
#[allow(dead_code)]
pub mod iters;
#[allow(dead_code)]
pub mod multiindex;
#[allow(dead_code)]
pub mod neighbors;
#[allow(dead_code)]
pub mod util;

use hashbrown::HashMap;
pub use iters::*;
pub use multiindex::*;
use nalgebra::Point;
pub use neighbors::*;
#[cfg(feature = "rayon")]
use rayon::prelude::ParallelIterator;
pub use util::*;

#[derive(Debug, Default)]
pub struct CellGrid<const N: usize> {
    cells: HashMap<[usize; N], usize>,
    cell_lists: Vec<Option<usize>>,
    index: MultiIndex<N>,
}

impl<const N: usize> CellGrid<N> {
    pub fn new(points: &[Point<f64, N>], cutoff: f64) -> Self {
        CellGrid::default().rebuild(points, Some(cutoff))
    }

    //TODO: Documentation: rebuild does not really update because it does not make a lot of sense:
    //TODO: If MultiIndex did change, we have to re-allocate (or re-initialize) almost everything anyway;
    //TODO: If MultiIndex did not change, we don't need to update.
    //TODO: Therefore we check for that and make CellGrid::new() just a wrapper around CellGrid::rebuild (with an initially empty MultiIndex)
    pub fn rebuild(self, points: &[Point<f64, N>], cutoff: Option<f64>) -> Self {
        let cutoff = cutoff.unwrap_or(self.index.grid_info.cutoff);
        let index = MultiIndex::from_points(points, cutoff);

        if index == self.index {
            self
        } else {
            let mut cells = HashMap::with_capacity(points.len());

            let cell_lists = index
                .index
                .iter()
                .enumerate()
                .map(|(i, cell)| cells.insert(*cell, i))
                .collect();

            Self {
                cells,
                cell_lists,
                index,
            }
        }
    }

    #[must_use = "rebuild_mut() consumes `self` and returns the mutated `CellGrid`"]
    pub fn rebuild_mut(mut self, points: &[Point<f64, N>], cutoff: Option<f64>) -> Self {
        if self.index.rebuild_mut(points, cutoff) {
            self.cell_lists.resize(points.len(), None);

            self.cells.clear();
            // We're not `reserve()`ing or `shrink_to_fit()`ing here.
            // HashMap should be sufficiently smart about reallocating in chunks.
            // Also this does not happen that often.

            // move out of `self`
            let index = self.index;
            let mut cell_lists = self.cell_lists;
            let mut cells = self.cells;

            index.index.iter().enumerate().for_each(|(i, cell)| {
                // see https://docs.rs/hashbrown/latest/hashbrown/struct.HashMap.html#method.insert
                // HashMap::insert() returns the old value or None, therefore we don't have to clear
                // cell_lists before refilling as long as cells is cleared.
                cell_lists[i] = cells.insert(*cell, i);
            });

            Self {
                cells,
                cell_lists,
                index,
            }
        } else {
            self
        }
    }

    pub fn shape(&self) -> &[usize; N] {
        &self.index.grid_info.shape
    }

    pub fn bounding_box(&self) -> &Aabb<N> {
        &self.index.grid_info.aabb
    }

    /// Iterate over all relevant (i.e. within cutoff threshold + some extra) unique pairs of points in this `CellGrid`.
    /// This has some overhead due to internal use of [`flat_map`].
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn point_pairs(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.iter()
            .flat_map(|cell| cell.point_pairs().collect::<Vec<_>>())
    }

    /// Call a closure on each relevant (i.e. within cutoff threshold + some extra) unique pair of points in this `CellGrid`.
    /// This should be preferred over [`CellGrid::point_pairs()`].
    /// `for_each_point_pair()` is equivalent to a nested loop:
    /// ```ignore
    /// for cell in cell_grid.iter() {
    ///     for (borb, other) in cell.point_pairs() {
    ///         ...
    ///     }
    /// }
    /// ```
    pub fn for_each_point_pair<F>(&self, mut f: F)
    where
        F: FnMut(usize, usize),
    {
        for cell in self.iter() {
            for (i, j) in cell.point_pairs() {
                f(i, j);
            }
        }
    }

    #[cfg(feature = "rayon")]
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn par_point_pairs(&self) -> impl ParallelIterator<Item = (usize, usize)> + '_ {
        //TODO: Find a way to handle cell lifetimes instead of collecting into a Vec?
        //TODO: seems to be related to flat_map()
        self.par_iter()
            .flat_map(|cell| cell.point_pairs().collect::<Vec<_>>())
    }
}

#[cfg(test)]
mod tests {
    //use super::*;

    #[test]
    fn test_cellgrid() {
        todo!()
    }
}

