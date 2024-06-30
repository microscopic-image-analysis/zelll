//TODO: impl GridIndex trait
//TODO: impl Deref/AsRef? to index
//TODO: maybe FlatIndex should know about a point cloud (i.e. store a reference?)
//TODO: Also: currently assuming that the order of points in point cloud does not change
//TODO: i.e. index in flatindex corresponds to index in point cloud
use crate::cellgrid::util::*;
use itertools::Itertools;
use std::borrow::Borrow;

#[derive(Debug, PartialEq, Default, Clone)]
pub struct FlatIndex<const N: usize> {
    pub(crate) grid_info: GridInfo<N>,
    pub(crate) index: Vec<i32>,
    pub(crate) neighbor_indices: Vec<i32>,
}

impl<const N: usize> FlatIndex<N> {
    pub fn with_capacity(info: GridInfo<N>, capacity: usize) -> Self {
        Self {
            grid_info: info,
            index: Vec::with_capacity(capacity),
            neighbor_indices: FlatIndex::neighbor_indices(info, 1),
        }
    }
    //TODO: this is a candidate for SIMD AoSoA
    //TODO: see https://www.rustsim.org/blog/2020/03/23/simd-aosoa-in-nalgebra/#using-simd-aosoa-for-linear-algebra-in-rust-ultraviolet-and-nalgebra
    //TODO: or can I chunk iterators such that rustc auto-vectorizes?
    //TODO: see https://www.nickwilcox.com/blog/autovec/
    pub fn from_points(
        points: impl IntoIterator<Item = impl Borrow<[f64; N]>> + Clone,
        cutoff: f64,
    ) -> Self {
        let aabb = Aabb::from_points(points.clone().into_iter());
        let grid_info = GridInfo::new(aabb, cutoff);
        let index = points
            .into_iter()
            .take(i32::MAX as usize)
            .map(|point| grid_info.flat_cell_index(point.borrow()))
            .collect();

        Self {
            grid_info,
            index,
            neighbor_indices: FlatIndex::neighbor_indices(grid_info, 1),
        }
    }
    // there is no rebuild(), named it rebuild_mut() to match CellGrid::rebuild_mut()
    //TODO: Documentation: return bool indicating whether the index changed at all (in length or any individual entry)
    pub fn rebuild_mut(
        &mut self,
        points: impl IntoIterator<Item = impl Borrow<[f64; N]>> + Clone,
        cutoff: Option<f64>,
    ) -> bool {
        let cutoff = cutoff.unwrap_or(self.grid_info.cutoff);
        let aabb = Aabb::from_points(points.clone().into_iter());
        let grid_info = GridInfo::new(aabb, cutoff);

        let size = points.clone().into_iter().take(i32::MAX as usize).count();
        self.index.resize(size, 0);

        let new_index = points
            .into_iter()
            .take(i32::MAX as usize)
            .map(|point| grid_info.flat_cell_index(point.borrow()));
        self.grid_info = grid_info;

        // FIXME: there seems to be a Bug in shape/stride computation, causing redundant indices for small shapes e.g. (2,2,2)?
        // FIXME: at least it used to, with RelativeNeighborIndices
        self.neighbor_indices = FlatIndex::neighbor_indices(grid_info, 1);

        let index_changed =
            self.index
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

// TODO: maybe make it part of the public API to allow changing the rank of the neighborhood
// TODO: but then we should make sure to re-scale the cell edge lengths, i.e. GridInfo needs to know about the neighborhood
// TODO: or store not the actual cutoff but rescaled. The public API should at most expose the option to set neighborhood rank
// TODO: this could easily handle full space as well (just have to ignore center of the neighborhood)
impl<const N: usize> FlatIndex<N> {
    fn neighbor_indices(grid_info: GridInfo<N>, rank: i32) -> Vec<i32> {
        // TODO: not sure if it makes sense to use u32 or NonZero<u32> here?
        assert!(rank > 0, "rank should be positive");

        (0..N)
            .map(|_| (-rank..rank + 1))
            .multi_cartesian_product()
            .map(|idx| grid_info.flatten_index(TryInto::<[i32; N]>::try_into(idx).unwrap()))
            //.take_while(|idx| *idx != 0) // not sure which one I like better
            .take((2 * rank + 1).pow(N as u32) as usize / 2) // equivalent to .div_euclid()
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_flatindex() {
        // using 0-origin for simplicity and to avoid floating point errors
        let points = generate_points([3, 3, 3], 1.0, [0.0, 0.0, 0.0]);
        let index = FlatIndex::from_points(points.iter(), 1.0);
        let mut idx = Vec::with_capacity(points.len());

        for x in 0..3 {
            for y in 0..3 {
                for z in 0..3 {
                    if (x + y + z) % 2 == 0 {
                        idx.push(index.grid_info.flatten_index([x, y, z]));
                        idx.push(index.grid_info.flatten_index([x, y, z]));
                    }
                }
            }
        }

        assert_eq!(index.index, idx, "testing FlatIndex::from_points()")
    }
}
