//TODO: impl GridIndex trait
//TODO: impl Deref/AsRef? to index
//TODO: maybe FlatIndex should know about a point cloud (i.e. store a reference?)
//TODO: Also: currently assuming that the order of points in point cloud does not change
//TODO: i.e. index in flatindex corresponds to index in point cloud
use crate::cellgrid::neighbors::RelativeNeighborIndices;
use crate::cellgrid::util::*;

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
    pub fn from_points<'p>(
        points: impl Iterator<Item = &'p [f64; N]> + Clone,
        cutoff: f64,
    ) -> Self {
        let aabb = Aabb::from_points(points.clone());
        let grid_info = GridInfo::new(aabb, cutoff);
        let index = points
            .take(i32::MAX as usize)
            .map(|point| grid_info.flat_cell_index(point))
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
    pub fn rebuild_mut<'p>(
        &mut self,
        points: impl Iterator<Item = &'p [f64; N]> + Clone,
        cutoff: Option<f64>,
    ) -> bool {
        let cutoff = cutoff.unwrap_or(self.grid_info.cutoff);
        let aabb = Aabb::from_points(points.clone());
        let grid_info = GridInfo::new(aabb, cutoff);

        //TODO:
        self.index
            .resize(points.clone().take(i32::MAX as usize).count(), 0);

        let new_index = points
            .take(i32::MAX as usize)
            .map(|point| grid_info.flat_cell_index(point));
        self.grid_info = grid_info;

        self.neighbor_indices = RelativeNeighborIndices::half_space()
            .map(|idx| grid_info.flatten_index(idx))
            .collect();

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
        let index = FlatIndex::from_points(points.iter().map(|p| p.coords.as_ref()), 1.0);
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


