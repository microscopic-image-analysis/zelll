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

pub use iters::*;
pub use multiindex::*;
use nalgebra::Point;
use ndarray::ArrayD;
pub use neighbors::*;
pub use util::*;

//TODO: I don't like this so far but a builder pattern is a bit overkill right now
#[derive(Debug, Default)]
pub struct CellGrid<const N: usize> {
    cells: ArrayD<Option<usize>>,
    //TODO: see https://crates.io/crates/stable-vec and https://crates.io/crates/slab
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
            let mut cell_lists: Vec<Option<usize>> = vec![None; points.len()];
            let mut cells: ArrayD<Option<usize>> =
                ArrayD::default(index.grid_info.shape.as_slice());

            index.index.iter().enumerate().for_each(|(i, cell)| {
                if let Some(head) = cells[cell.as_slice()] {
                    cell_lists[i] = Some(head);
                }
                cells[cell.as_slice()] = Some(i);
            });

            Self {
                cells,
                cell_lists,
                index,
            }
        }
    }

    pub fn rebuild_mut(mut self, points: &[Point<f64, N>], cutoff: Option<f64>) -> Self {
        if self.index.rebuild_mut(points, cutoff) {
            self.cell_lists.clear();
            self.cell_lists.resize(points.len(), None);
            self.cell_lists.shrink_to_fit();

            //TODO: see https://docs.rs/ndarray/latest/ndarray/type.ArrayView.html#method.from_shape
            //TODO: and https://docs.rs/ndarray/latest/ndarray/struct.ArrayBase.html#method.from_shape_vec
            //TODO: and https://docs.rs/ndarray/latest/ndarray/struct.ArrayBase.html#method.into_raw_vec
            //TODO: and https://github.com/rust-ndarray/ndarray/issues/433#issuecomment-377288065
            //TODO: I need to do rebuild_mut(mut self, ...) (vs. &mut self) in this case and return the consumed struct
            //TODO: this should be documented.
            let mut cells = self.cells.into_raw_vec();

            cells.clear();
            cells.resize(self.index.grid_info.shape.iter().product(), None);
            cells.shrink_to_fit();

            // Panicking is okay here. But this *should* not happen 😬
            let mut cells = ArrayD::from_shape_vec(self.index.grid_info.shape.as_slice(), cells)
                .expect("ArrayD::from_shape_vec() failed");

            // complete the partial move out of `self`
            let index = self.index;
            let mut cell_lists = self.cell_lists;

            index.index.iter().enumerate().for_each(|(i, cell)| {
                if let Some(head) = cells[cell.as_slice()] {
                    cell_lists[i] = Some(head);
                }
                cells[cell.as_slice()] = Some(i);
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

    /// Iterate over all relevant (i.e. within cutoff threshold + some extra) unique pairs of points in this `CellGrid`
    pub fn point_pairs(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        //TODO: Find a way to handle cell lifetimes instead of collecting into a Vec?
        //TODO: seems to be related to flat_map()
        self.iter()
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

