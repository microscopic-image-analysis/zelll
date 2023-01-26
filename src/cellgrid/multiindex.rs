//TODO: impl Deref/AsRef? to index
//TODO: maybe MultiIndex should know about a point cloud (i.e. store a reference?)
//TODO: Also: currently assuming that the order of points in point cloud does not change
//TODO: i.e. index in multiindex corresponds to index in point cloud
//TODO: maybe CellGrid should have this (immutable reference), therefore enforcing that
//TODO: changes to a pointcloud can only be done if CellGrid is out of scope (all references are dropped)
use crate::cellgrid::util::*;
use nalgebra::Point;

#[derive(Debug)]
pub struct MultiIndex<const N: usize> {
    pub(crate) grid_info: GridInfo<N>,
    pub(crate) index: Vec<[usize; N]>,
}

impl<const N: usize> MultiIndex<N> {
    pub fn with_capacity(info: GridInfo<N>, capacity: usize) -> Self {
        Self {
            grid_info: info,
            index: Vec::with_capacity(capacity),
        }
    }

    pub fn from_points(points: &[Point<f64, N>], cutoff: f64) -> Self {
        let aabb = Aabb::from_points(points);
        let grid_info = GridInfo::new(aabb, cutoff);

        let index = points
            .iter()
            .map(|point| grid_info.cell_index(point))
            .collect();

        Self { grid_info, index }
    }

    //TODO: should I allow to change the cutoff in this or a similar method?
    //TODO: might be nice and save some allocations
    //TODO: but one could argue changing the cutoff would justify constructing a new multi index
    //TODO: but I copy grid info any way, so might as well change it
    pub fn update_indices(&mut self, points: &[Point<f64, N>]) {
        //`GridInfo` is `Copy`
        // we need to copy here since partial borrowing is not possible
        //TODO: this really depends on which struct should handle cell_index() computation
        //TODO: but I probably shouldn't think about it too much, rustc will probably optimize it away
        //TODO: but we might have to benchmark it at some point
        //TODO: also see: https://www.forrestthewoods.com/blog/should-small-rust-structs-be-passed-by-copy-or-by-borrow/
        let info = self.grid_info;
        self.index
            .iter_mut()
            .zip(points.iter())
            .for_each(|(idx, point)| *idx = info.cell_index(point));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_multiindex() {
        // using 0-origin for simplicity and to avoid floating point errors
        let points = generate_points([3, 3, 3], 1.0, [0.0, 0.0, 0.0]);

        let mut idx = Vec::with_capacity(points.len());

        for x in 0..3 {
            for y in 0..3 {
                for z in 0..3 {
                    if (x + y + z) % 2 == 0 {
                        idx.push([x, y, z]);
                        idx.push([x, y, z]);
                    }
                }
            }
        }

        let index = MultiIndex::from_points(&points, 1.0);

        assert_eq!(index.index, idx, "testing MultiIndex::from_points()")
    }
}

