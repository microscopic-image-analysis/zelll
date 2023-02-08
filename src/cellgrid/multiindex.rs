//TODO: impl Deref/AsRef? to index
//TODO: maybe MultiIndex should know about a point cloud (i.e. store a reference?)
//TODO: Also: currently assuming that the order of points in point cloud does not change
//TODO: i.e. index in multiindex corresponds to index in point cloud
use crate::cellgrid::util::*;
use nalgebra::Point;

#[derive(Debug, PartialEq, Default)]
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
    // there is no rebuild(), named it rebuild_mut() to match CellGrid::rebuild_mut()
    //TODO: Documentation: return bool indicating whether the index changed at all (in length or any individual entry)
    pub fn rebuild_mut(&mut self, points: &[Point<f64, N>], cutoff: Option<f64>) -> bool {
        let cutoff = cutoff.unwrap_or(self.grid_info.cutoff);
        let aabb = Aabb::from_points(points);
        let grid_info = GridInfo::new(aabb, cutoff);

        self.index.resize(points.len(), [0; N]);

        let new_index = points.iter().map(|point| grid_info.cell_index(point));

        // Using `|` here because `||` is lazy and we always need to run the iterator
        let index_changed = (points.len() != self.index.len())
            | self
                .index
                .iter_mut()
                .zip(new_index)
                .fold(false, |has_changed, (old, new)| {
                    if *old != new {
                        *old = new;
                        true
                    } else {
                        has_changed
                    }
                });
        index_changed
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

