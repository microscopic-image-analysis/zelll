#[allow(dead_code)]
pub mod iters;
#[allow(dead_code)]
pub mod multiindex;
#[allow(dead_code)]
pub mod util;

pub use iters::*;
pub use multiindex::*;
use ndarray::Array3;
pub use util::*;

//TODO: I don't like this so far but a builder pattern is a bit overkill right now
pub struct CellGrid {
    points: PointCloud,
    cells: Array3<Option<usize>>,
    //TODO: see https://crates.io/crates/stable-vec and https://crates.io/crates/slab
    cell_lists: Vec<Option<usize>>,
    index: MultiIndex,
}

impl CellGrid {
    pub fn new(points: &[[f64; 3]], cutoff: f64) -> Self {
        let points = PointCloud::from_points(points);
        let index = MultiIndex::from_pointcloud(&points, cutoff);

        let mut cell_lists: Vec<Option<usize>> = Vec::default();
        let mut cells: Array3<Option<usize>> = Array3::default(index.grid_info.shape);

        index.index.iter().enumerate().for_each(|(i, cell)| {
            if let Some(head) = cells[*cell] {
                cell_lists[i] = Some(head);
            }
            cells[*cell] = Some(i);
        });

        Self {
            points,
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
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cellgrid() {
        todo!()
    }
}

