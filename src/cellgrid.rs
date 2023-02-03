//TODO: handle updating of CellGrid and inner fields with point cloud data of different sizes
//TODO: - currently assuming particle indices are stable
//TODO: - handle removal/addition of particles -> handle update or enforce new construction of CellGrid
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
#[derive(Debug)]
pub struct CellGrid<const N: usize> {
    cells: ArrayD<Option<usize>>,
    //TODO: see https://crates.io/crates/stable-vec and https://crates.io/crates/slab
    cell_lists: Vec<Option<usize>>,
    index: MultiIndex<N>,
}

impl<const N: usize> CellGrid<N> {
    pub fn new(points: &[Point<f64, N>], cutoff: f64) -> Self {
        let index = MultiIndex::from_points(points, cutoff);

        let mut cell_lists: Vec<Option<usize>> = [None].repeat(points.len());
        let mut cells: ArrayD<Option<usize>> = ArrayD::default(index.grid_info.shape.as_slice());

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

    //TODO: either &mut self or move unchanged fields to returned new struct
    //TODO: see https://doc.rust-lang.org/book/ch05-01-defining-structs.html#creating-instances-from-other-instances-with-struct-update-syntax
    fn update(self) -> Self {
        todo!()
    }

    /// Iterate over all relevant (i.e. within cutoff threshold + some extra) unique pairs of points in this `CellGrid`
    pub fn point_pairs(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        //TODO: Find a way to handle cell lifetimes instead of collecting into a Vec?
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

