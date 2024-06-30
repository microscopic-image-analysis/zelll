//TODO: impl GridIndex trait
//TODO: impl Deref/AsRef? to index
//TODO: maybe FlatIndex should know about a point cloud (i.e. store a reference?)
//TODO: Also: currently assuming that the order of points in point cloud does not change
//TODO: i.e. index in flatindex corresponds to index in point cloud
use crate::cellgrid::neighbors::RelativeNeighborIndices;
use crate::cellgrid::util::*;
use std::borrow::Borrow;
//use itertools::Itertools;

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
            neighbor_indices: RelativeNeighborIndices::half_space()
                .map(|idx| info.flatten_index(idx))
                .collect(),
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
            //TODO: it's a bit annoying that we're handling half-/full-space here, i.e. we currently would need to re-compute the FlatIndex to switch. This is not very nice behaviour.
            neighbor_indices: RelativeNeighborIndices::half_space()
                .map(|idx| grid_info.flatten_index(idx))
                .collect(),
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

        // TODO: this actually boils down to just enumerating from 1 to 13
        // TODO: and should also work for larger neighbor "kernels" (e.g. order 5x5 instead of 3x3)
        // TODO: just have to compute neighbors = floor(order^N/2)
        // TODO: (probably don't need floor(), could just do +1 somewhere since order is odd)
        // TODO: and could maybe just store the range 1..=neighbors (although I'd prefer a exclusive range)
        // FIXME: there seems to be a Bug in this map()/flatten_index() or RelativeNeighborIndices?
        self.neighbor_indices = RelativeNeighborIndices::half_space()
            .map(|idx| grid_info.flatten_index(idx))
            .collect();
        // FIXME: it doesn't produce the same as this and sometimes even redundant values
        // eprintln!("{:?}", self.neighbor_indices);
        // self.neighbor_indices = [-2, -1, 0, 1, 2].into_iter()
        //     .map(|_| [-2, -1, 0, 1, 2].into_iter() )
        //     .multi_cartesian_product()
        //     .map(|idx| grid_info.flatten_index(TryInto::<[i32; N]>::try_into(idx).unwrap()))
        //     .filter(|idx| idx.is_positive())
        //     .collect();
        // eprintln!("{:?}", self.neighbor_indices);
        // [1, 1, 2, 3, 3, 4, 5, 5, 6, 7, 7, 8]
        // [3, 5, 1, 7, 4, 6, 2, 8, 5, 1, 7, 3, 9]
        //
        // [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
        // [5, 8, 2, 11, 6, 9, 3, 12, 7, 1, 10, 4, 13]
        // FIXME: first, RelativeNeighbors omits the last relative neighbor? check for 2D
        // FIXME: second, redundant flat indices indicate a bug in computing shape/strides
        // FIXME: probably my +1 hack/trick. so it should be revisited
        // FIXME: especially in light of allowing higher order neighborhoods (with cells of edge length cutoff/(order-1/2))

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
