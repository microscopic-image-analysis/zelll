//TODO: make everything a Point3 or [;3] for now?
//TODO: mixing it seems wrong but somehow it makes sense to use either of both types
//TODO: in different situations. [;3] -> Point3 is easy, reverse requires some boilerplate

#![allow(dead_code)]
use nalgebra::*;

//TODO: For now we're just copying into our own struct for simplicity
//TODO: There's not much use to this yet
//TODO: maybe I should just keep this a type alias
//TODO: which would allow me to just use &[Point3<f64>] in function signatures
//TODO: also keep in mind https://rust-unofficial.github.io/patterns/anti_patterns/deref.html
//TODO: rather impl AsRef https://doc.rust-lang.org/std/convert/trait.AsRef.html
//TODO: #[repr(transparent)]?
//TODO: see https://doc.rust-lang.org/nomicon/other-reprs.html#reprtransparent
#[derive(Debug)]
pub struct PointCloud<const N: usize>(pub(crate) Vec<Point<f64, N>>);

impl<const N: usize> PointCloud<N> {
    //TODO: we'll want a more generic way for this (maybe impl Deref/AsRef)
    pub fn from_points(points: &[[f64; N]]) -> Self {
        Self(points.iter().map(|p| Point::from(*p)).collect())
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

//TODO: at some point, if I mutate points in a point cloud and do not rebuild the cell grid entirely,
//TODO: I'd probably need to check against the concrete Aabb instance to make sure points stay in the grid
//TODO: or rather: if points stay in the grid, I won't have to rebuild the cell grid
#[derive(Clone, Copy, Debug)]
pub struct Aabb<const N: usize> {
    inf: Point<f64, N>,
    sup: Point<f64, N>,
}

impl<const N: usize> Aabb<N> {
    pub fn from_pointcloud(point_cloud: &PointCloud<N>) -> Self {
        let init = if point_cloud.is_empty() {
            Point::<f64, N>::default()
        } else {
            point_cloud.0[0]
        };

        let (inf, sup) = point_cloud
            .0
            .iter()
            .fold((init, init), |(i, s), point| (i.inf(point), s.sup(point)));

        Self { inf, sup }
    }
}

/// The grid described by `GridInfo` may be slightly larger than the underlying bounding box `aabb`.
#[derive(Clone, Copy, Debug)]
pub struct GridInfo<const N: usize> {
    aabb: Aabb<N>,
    cutoff: f64,
    //TODO: probably should implement a method instead of using pub/pub(crate)
    pub(crate) shape: [usize; N],
}

impl<const N: usize> GridInfo<N> {
    pub fn new(aabb: Aabb<N>, cutoff: f64) -> Self {
        // TODO: not sure yet if I want shape to be a Point3
        let mut shape = [0; N];
        // TODO: This is not very nice yet. We'll figure the precise types out later
        shape.copy_from_slice(
            ((aabb.sup - aabb.inf) / cutoff)
                .map(|coord| coord.floor() as usize + 1)
                .as_slice(),
        );

        Self {
            aabb,
            cutoff,
            shape,
        }
    }

    pub fn origin(&self) -> &Point<f64, N> {
        &self.aabb.inf
    }

    //TODO: not sure where it fits better
    //TODO: GridInfo knows enough to do compute the cell index for an arbitrary point
    //TODO: but MultiIndex seems more fitting semantically?
    pub fn cell_index(&self, point: &Point<f64, N>) -> [usize; N] {
        let mut idx = [0; N];

        idx.copy_from_slice(
            ((point - self.origin()) / self.cutoff)
                .map(|coord| coord.floor() as usize)
                .as_slice(),
        );

        idx
    }
}

#[cfg(test)]
mod tests {
    //use super::*;

    #[test]
    fn test_gridinfo() {
        /*let point_cloud = PointCloud::from_points(&POINTS);
        let aabb = Aabb::from_pointcloud(&point_cloud);
        let grid_info = GridInfo::new(aabb, 5.0);
        println!("{:?}", grid_info);*/
        todo!("test");
    }
}

