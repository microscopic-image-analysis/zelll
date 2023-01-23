//TODO: impl Deref/AsRef? to index
//TODO: maybe MultiIndex should own a point cloud and I should provide methods to deref to point cloud
//TODO: and/or deconstruct to underlying point cloud
//TODO: implementing From/Into is not trivial here because point cloud has no knowledge about the cutoff
//TODO: Also: this currently assumes that the order of points in point cloud does not change
//TODO: i.e. index in multiindex corresponds to index in point cloud
use crate::cellgrid::util::*;

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

    pub fn from_pointcloud(point_cloud: &PointCloud<N>, cutoff: f64) -> Self {
        let aabb = Aabb::from_pointcloud(point_cloud);
        let grid_info = GridInfo::new(aabb, cutoff);

        let index = point_cloud
            .0
            .iter()
            .map(|point| grid_info.cell_index(point))
            .collect();

        Self { grid_info, index }
    }

    //TODO: should I allow to change the cutoff in this or a similar method?
    //TODO: might be nice and save some allocations
    //TODO: but one could argue changing the cutoff would justify constructing a new multi index
    //TODO: but I copy grid info any way, so might as well change it
    pub fn update_indices(&mut self, point_cloud: &PointCloud<N>) {
        //`GridInfo` is `Copy`
        // we need to copy here since partial borrowing is not possible
        //TODO: this really depends on which struct should handle cell_index() computation
        //TODO: but I probably shouldn't think about it too much, rustc will probably optimize it away
        //TODO: but we might have to benchmark it at some point
        //TODO: also see: https://www.forrestthewoods.com/blog/should-small-rust-structs-be-passed-by-copy-or-by-borrow/
        let info = self.grid_info;
        self.index
            .iter_mut()
            .zip(point_cloud.0.iter())
            .for_each(|(idx, point)| *idx = info.cell_index(point));
    }
}

