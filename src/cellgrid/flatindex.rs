//TODO: currently assuming that the order of points in point cloud does not change
//TODO: i.e. index in flatindex corresponds to index in point cloud, this should be documented
use super::util::{Aabb, GridInfo};
use crate::Particle;
use itertools::Itertools;
use nalgebra::SimdPartialOrd;
use num_traits::{AsPrimitive, ConstOne, ConstZero, Float, NumAssignOps};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::borrow::Borrow;

#[derive(Debug, PartialEq, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub(crate) struct FlatIndex<const N: usize = 3, F: Float = f64>
where
    F: std::fmt::Debug + 'static,
{
    pub(crate) grid_info: GridInfo<N, F>,
    pub(crate) index: Vec<i32>,
    pub(crate) neighbor_indices: Vec<i32>,
}

impl<const N: usize, F> Default for FlatIndex<N, F>
where
    F: Float + std::fmt::Debug + Default + NumAssignOps + AsPrimitive<i32> + ConstOne,
{
    fn default() -> Self {
        FlatIndex::with_capacity(GridInfo::default(), 0)
    }
}

impl<const N: usize, F> FlatIndex<N, F>
where
    F: Float + std::fmt::Debug,
{
    pub fn with_capacity(info: GridInfo<N, F>, capacity: usize) -> Self {
        Self {
            grid_info: info,
            index: Vec::with_capacity(capacity),
            neighbor_indices: FlatIndex::neighbor_indices(info),
        }
    }

    // TODO: maybe make it part of the public API to allow changing the rank of the neighborhood
    // TODO: but then we should make sure to re-scale the cell edge lengths, i.e. GridInfo needs to know about the neighborhood
    // TODO: GridInfo should store rank and cutoff, so neighbor_indices can access it.
    // TODO: However, for rank > 1, HashMap is not the best choice anymore. for N=3: 62 vs 13 neighboring cells (half-space)
    // TODO: that many random lookups add up. Also, more non-empty cells makes HashMap construction more expensive
    // TODO: for higher ranks, we'd sth. with more spatial locality
    // TODO:
    // FIXME: handling full space is a bit more expensive (due to filtering out center of the neighborhood)
    // FIXME: in this case (-rank..rank+1) is not quite ideal. sth. like (0..rank+1).chain(-rank..0)
    // FIXME: and then skip first element of cartesian product
    // FIXME: but this is not noticeable for non-trivially sized point clouds
    fn neighbor_indices(grid_info: GridInfo<N, F>) -> Vec<i32> {
        // this is the rank of the neighborhood, 1 -> 3^N, 2 -> 5^N
        const RANK: i32 = 1;

        (0..N)
            .map(|_| -RANK..RANK + 1)
            .multi_cartesian_product()
            .map(|idx| grid_info.flatten_index(TryInto::<[i32; N]>::try_into(idx).unwrap()))
            .filter(|idx| *idx != 0)
            .collect()
    }
}

impl<const N: usize, F> FlatIndex<N, F>
where
    F: Float + ConstZero + NumAssignOps + SimdPartialOrd + AsPrimitive<i32> + std::fmt::Debug,
{
    //TODO: this is a candidate for SIMD AoSoA
    //TODO: see https://www.rustsim.org/blog/2020/03/23/simd-aosoa-in-nalgebra/#using-simd-aosoa-for-linear-algebra-in-rust-ultraviolet-and-nalgebra
    //TODO: or can I chunk iterators such that rustc auto-vectorizes?
    //TODO: see https://www.nickwilcox.com/blog/autovec/
    pub fn from_particles<P: Particle<[F; N]>>(
        particles: impl IntoIterator<Item = impl Borrow<P>> + Clone,
        cutoff: F,
    ) -> Self {
        // TODO: We could actually use a fixed Aabb (cf. util::GridInfo::new_fixed())
        // TODO: However, computing the bounding box is cheap enough
        // TODO: We might reconsider this after experimenting with hashbrown::HashTable
        let aabb = Aabb::from_particles(particles.clone().into_iter());
        let grid_info = GridInfo::new(aabb, cutoff);
        let index = particles
            .into_iter()
            .take(i32::MAX as usize)
            .map(|p| grid_info.flat_cell_index(p.borrow().coords()))
            .collect();

        // TODO: this does seem to have a *small* effect due to autovectorization?
        // TODO: (although often buried by cache effects)
        // TODO: should examine more closely
        // let mut it = particles.into_iter()
        //     .map(|p| grid_info.flat_cell_index(p.borrow().coords()));
        // let mut index: Vec<i32> = Vec::with_capacity(it.size_hint().0);

        // while let Some(chunk) = it.next_array::<16>() {
        //     index.extend_from_slice(&chunk);
        // }

        // index.extend(it);

        Self {
            grid_info,
            index,
            neighbor_indices: FlatIndex::neighbor_indices(grid_info),
        }
    }
    // there is no rebuild(), named it rebuild_mut() to match CellGrid::rebuild_mut()
    // TODO: Documentation: return bool indicating whether the index changed at all (in length or any individual entry)
    // TODO: benchmark with changing point iterators
    pub fn rebuild_mut<P: Particle<[F; N]>>(
        &mut self,
        particles: impl IntoIterator<Item = impl Borrow<P>> + Clone,
        cutoff: Option<F>,
    ) -> bool {
        let cutoff = cutoff.unwrap_or(self.grid_info.cutoff);
        // TODO: see TODO notes for ::from_points()
        let aabb = Aabb::from_particles(particles.clone().into_iter());
        let grid_info = GridInfo::new(aabb, cutoff);

        let size = particles
            .clone()
            .into_iter()
            .take(i32::MAX as usize)
            .count();
        // TODO: should benchmark this
        // let size_changed = size != self.index.len();
        self.index.resize(size, 0);

        let new_index = particles
            .into_iter()
            .take(i32::MAX as usize)
            .map(|p| grid_info.flat_cell_index(p.borrow().coords()));

        self.grid_info = grid_info;
        self.neighbor_indices = FlatIndex::neighbor_indices(grid_info);

        self.index
            .iter_mut()
            .zip(new_index)
            .fold(false, |has_changed, (old, new)| {
                if
                /* size_changed || */
                *old != new {
                    *old = new;
                    true
                } else {
                    has_changed
                }
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::CellGrid;
    use crate::cellgrid::util::generate_pointcloud;

    #[test]
    fn test_neighbor_indices() {
        // cf. `util.rs::GridInfo::new()`
        // to see why this corresponds to a 8x8 chessboard
        let points = vec![[0.0, 0.0], [3.0, 3.0]];
        let cg = CellGrid::new(points.iter().copied(), 1.0);
        let indices = FlatIndex::neighbor_indices(*cg.info());

        assert_eq!(indices, vec![-9, -1, 7, -8, 8, -7, 1, 9]);
    }

    #[test]
    fn test_flatindex() {
        // using 0-origin for simplicity and to avoid floating point errors
        let points = generate_pointcloud([3, 3, 3], 1.0, [0.0, 0.0, 0.0]);
        let index = FlatIndex::from_particles::<[_; 3]>(points.iter(), 1.0);
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
