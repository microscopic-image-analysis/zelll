#![allow(dead_code)]
use nalgebra::*;

pub type PointCloud<const N: usize> = Vec<Point<f64, N>>;

#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct Aabb<const N: usize> {
    pub inf: Point<f64, N>,
    pub sup: Point<f64, N>,
}

impl<const N: usize> Aabb<N> {
    pub fn from_points<'p>(mut points: impl Iterator<Item = &'p Point<f64, N>>) -> Self {
        let init = points.next().copied().unwrap_or_default();

        let (inf, sup) = points
            .take(i32::MAX as usize) //TODO: this works but maybe explicit try_into() would be better?
            .fold((init, init), |(i, s), point| (i.inf(point), s.sup(point)));

        Self { inf, sup }
    }
}

/// The grid described by `GridInfo` may be slightly larger than the underlying bounding box `aabb`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GridInfo<const N: usize> {
    pub(crate) aabb: Aabb<N>,
    pub(crate) cutoff: f64,
    //TODO: probably should implement a method instead of using pub/pub(crate)
    pub(crate) shape: [i32; N],
    pub(crate) strides: [i32; N],
}

impl<const N: usize> GridInfo<N> {
    pub fn new(aabb: Aabb<N>, cutoff: f64) -> Self {
        // TODO: not sure yet if I want shape to be a Point<N>
        let mut shape = [0; N];
        // TODO: This is not very nice yet. We'll figure the precise types out later
        shape.copy_from_slice(
            ((aabb.sup - aabb.inf) / cutoff)
                .map(|coord| coord.floor() as i32 + 1)
                .as_slice(),
        );

        let mut strides = shape;
        strides.iter_mut().fold(1, |prev, curr| {
            //TODO: simulating larger shape to increase strides; this allows for relative negative indices,
            //TODO: at least those represented by BalancedTernary<N>
            //TODO: could also simply increase the shape, maybe that's less confusing
            //TODO: and it doesn't affect memory since we're using a hash map anyway
            //TODO: the better approach would be:
            //TODO: padded shape (i.e. (2,3,4) -> (4,5,6))
            //TODO: compute strides from padded shape (i.e. instead of (5,5)->(1,5) or (1,5+1), do (7,7)->(1,7))
            //TODO: compute cell_index() from 1-based multi-index instead of 0-based, i.e. lower left corner is [1; N] instead of [0; N]
            let next = prev * (*curr + 1);
            *curr = prev;
            next
        });

        Self {
            aabb,
            cutoff,
            shape,
            strides,
        }
    }

    pub fn origin(&self) -> &Point<f64, N> {
        &self.aabb.inf
    }

    //TODO: not sure where it fits better
    //TODO: GridInfo knows enough to compute the cell index for an arbitrary point
    //TODO: but might make more sense in FlatIndex?
    //TODO: sth. like Lattice trait maybe
    pub fn cell_index(&self, point: &Point<f64, N>) -> [i32; N] {
        let mut idx = [0; N];

        idx.copy_from_slice(
            ((point - self.origin()) / self.cutoff)
                .map(|coord| coord.floor() as i32)
                .as_slice(),
        );

        assert!(
            idx < self.shape,
            "cell index ({:?}) out of bounds ({:?})",
            idx,
            self.shape
        );

        idx
    }

    pub fn flat_cell_index(&self, point: &Point<f64, N>) -> i32 {
        self.flatten_index(self.cell_index(point))
    }

    pub fn flatten_index(&self, idx: [i32; N]) -> i32 {
        //TODO: benchmark if nalgebra dot product might be faster (due to SIMD)
        idx.iter()
            .zip(self.strides)
            .map(|(i, s)| *i * s)
            .sum::<i32>()
    }
}

impl<const N: usize> Default for GridInfo<N> {
    fn default() -> Self {
        GridInfo::new(Aabb::default(), 1.0)
    }
}

/// Generate generate 3-dimensional point arrays for testing purposes in the following fashion:
/// In a grid with cells of length `cutoff` only cells with even linear index contain points (chessboard pattern).
/// These non-empty cells contain two points each:
/// - the first at the origin of the cell (equivalent to the cell's multi-index + the origin of the grid)
/// - the second at the center of the cell
// We'll stay in 3D for simplicity here
pub fn generate_points(shape: [usize; 3], cutoff: f64, origin: [f64; 3]) -> PointCloud<3> {
    let mut points = Vec::with_capacity(((shape.iter().product::<usize>() + 1) / 2) * 2);

    for x in 0..shape[0] {
        for y in 0..shape[1] {
            for z in 0..shape[2] {
                if (x + y + z) % 2 == 0 {
                    points.push(Point::from([
                        cutoff.mul_add(x as f64, origin[0]),
                        cutoff.mul_add(y as f64, origin[1]),
                        cutoff.mul_add(z as f64, origin[2]),
                    ]));
                    points.push(Point::from([
                        cutoff.mul_add(x as f64, cutoff.mul_add(0.5, origin[0])),
                        cutoff.mul_add(y as f64, cutoff.mul_add(0.5, origin[1])),
                        cutoff.mul_add(z as f64, cutoff.mul_add(0.5, origin[2])),
                    ]));
                }
            }
        }
    }

    points
}

/// Generate a uniformly random 3D point cloud of size `n` in a cuboid of edge lengths `vol` centered around `origin`.
//TODO: Matrix::new_random() requires nalgebra with feature `rand`
//TODO: I'm not sure why cargo check and cargo clippy are not using dev dependencies?
//TODO: but I probably should set some feature in Cargo.toml if I want to keep this public
#[cfg(feature = "dep:nalgebra/rand")]
pub(crate) fn generate_points_random(n: usize, vol: [f64; 3], origin: [f64; 3]) -> PointCloud<3> {
    std::iter::repeat_with(|| {
        Point3::<f64>::from(
            (Vector3::new_random() - Vector3::new(0.5, 0.5, 0.5) + Vector3::from(origin))
                .component_mul(&Vector3::from(vol)),
        )
    })
    .take(n)
    .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_points() {
        let points = vec![
            Point::from([0.0, 0.0, 0.0]),
            Point::from([0.5, 0.5, 0.5]),
            Point::from([0.0, 0.0, 2.0]),
            Point::from([0.5, 0.5, 2.5]),
            Point::from([0.0, 1.0, 1.0]),
            Point::from([0.5, 1.5, 1.5]),
            Point::from([0.0, 2.0, 0.0]),
            Point::from([0.5, 2.5, 0.5]),
            Point::from([0.0, 2.0, 2.0]),
            Point::from([0.5, 2.5, 2.5]),
            Point::from([1.0, 0.0, 1.0]),
            Point::from([1.5, 0.5, 1.5]),
            Point::from([1.0, 1.0, 0.0]),
            Point::from([1.5, 1.5, 0.5]),
            Point::from([1.0, 1.0, 2.0]),
            Point::from([1.5, 1.5, 2.5]),
            Point::from([1.0, 2.0, 1.0]),
            Point::from([1.5, 2.5, 1.5]),
            Point::from([2.0, 0.0, 0.0]),
            Point::from([2.5, 0.5, 0.5]),
            Point::from([2.0, 0.0, 2.0]),
            Point::from([2.5, 0.5, 2.5]),
            Point::from([2.0, 1.0, 1.0]),
            Point::from([2.5, 1.5, 1.5]),
            Point::from([2.0, 2.0, 0.0]),
            Point::from([2.5, 2.5, 0.5]),
            Point::from([2.0, 2.0, 2.0]),
            Point::from([2.5, 2.5, 2.5]),
        ];
        assert_eq!(points, generate_points([3, 3, 3], 1.0, [0.0, 0.0, 0.0]));
    }

    #[test]
    fn test_utils() {
        let points = generate_points([3, 3, 3], 1.0, [0.2, 0.25, 0.3]);
        assert_eq!(points.len(), 28, "testing PointCloud.len()");

        let aabb = Aabb::from_points(points.iter());
        assert_eq!(
            aabb,
            Aabb {
                inf: [0.2, 0.25, 0.3].into(),
                sup: [2.7, 2.75, 2.8].into()
            },
            "testing Aabb::from_pointcloud()"
        );

        let grid_info = GridInfo::new(aabb, 1.0);
        assert_eq!(
            grid_info.origin(),
            &Point::from([0.2, 0.25, 0.3]),
            "testing GridInfo.origin()"
        );
        assert_eq!(grid_info.shape, [3, 3, 3], "testing GridInfo.shape");
        //TODO: note that these are the strides for grid_info.shape + [1, 1, 1]
        //TODO: this allows us to use negative relative indices representable by BalancedTernary<N>
        assert_eq!(grid_info.strides, [1, 4, 16], "testing GridInfo.strides");

        // Intuitively you'd expect [2, 2, 2] for this
        // but we're having floating point imprecision:
        // 2.3 - 0.3 = 1.9999999999999998
        // This is not an issue though because the index is still uniquely determined
        assert_eq!(
            grid_info.cell_index(&Point::from([2.7, 2.75, 2.3])),
            [2, 2, 1],
            "testing GridInfo.cell_index()"
        );
        assert_eq!(
            grid_info.flat_cell_index(&Point::from([2.7, 2.75, 2.3])),
            26,
            "testing GridInfo.flat_cell_index()"
        );
        assert_eq!(
            grid_info.cell_index(&Point::from([2.7, 2.75, 2.8])),
            [2, 2, 2],
            "testing GridInfo.cell_index()"
        );
        assert_eq!(
            grid_info.flat_cell_index(&Point::from([2.7, 2.75, 2.8])),
            42,
            "testing GridInfo.flat_cell_index()"
        );
    }
}

