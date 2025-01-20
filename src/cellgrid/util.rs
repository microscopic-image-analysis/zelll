#![allow(dead_code)]
use crate::Particle;
use nalgebra::*;
use num_traits::{AsPrimitive, ConstOne, ConstZero, Float, NumAssignOps};
use std::borrow::Borrow;

// TODO: remove this alias?
pub type PointCloud<const N: usize> = Vec<[f64; N]>;

//TODO: rename fields, infimum/supremum might be confusing outside of a lattice context
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct Aabb<const N: usize = 3, F: Float = f64>
where
    F: std::fmt::Debug + 'static,
{
    inf: Point<F, N>,
    sup: Point<F, N>,
}

pub type Aabb64<const N: usize> = Aabb<N, f64>;
pub type Aabb32<const N: usize> = Aabb<N, f32>;

impl<const N: usize, F> Aabb<N, F>
where
    F: Float + std::fmt::Debug + SimdPartialOrd + ConstZero,
{
    pub fn from_points<P: Particle<[F; N]>>(
        mut points: impl Iterator<Item = impl Borrow<P>>,
    ) -> Self {
        let init = points
            .next()
            .map(|p| p.borrow().coords())
            .unwrap_or([F::ZERO; N]);
        let init = Point::from(init);

        let (inf, sup) = points
            .take(i32::MAX as usize)
            .fold((init, init), |(i, s), point| {
                let point = Point::from(point.borrow().coords());
                (i.inf(&point), s.sup(&point))
            });

        Self { inf, sup }
    }

    //TODO: could also pass iterators here (single point could be wrapped by std::iter::once or Option::iter())
    fn update<P: Particle<[F; N]>>(&mut self, point: impl Borrow<P>) {
        let point = Point::from(point.borrow().coords());
        self.inf = point.inf(&self.inf);
        self.sup = point.sup(&self.sup);
    }

    pub fn inf(&self) -> [F; N] {
        self.inf.into()
    }

    pub fn sup(&self) -> [F; N] {
        self.sup.into()
    }
}

/// The grid described by `GridInfo` may be slightly larger than the underlying bounding box `aabb`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GridInfo<const N: usize = 3, F: Float = f64>
where
    F: std::fmt::Debug + 'static,
{
    pub(crate) aabb: Aabb<N, F>,
    pub(crate) cutoff: F,
    shape: SVector<i32, N>,
    strides: SVector<i32, N>,
}

pub type GridInfo64<const N: usize> = GridInfo<N, f64>; // = GridInfo<N> would suffice but let's be explicit here
pub type GridInfo32<const N: usize> = GridInfo<N, f32>;

impl<const N: usize, F> GridInfo<N, F>
where
    F: Float + std::fmt::Debug,
{
    pub fn origin(&self) -> [F; N] {
        self.aabb.inf.into()
    }

    pub fn shape(&self) -> [i32; N] {
        self.shape.into()
    }

    pub fn strides(&self) -> [i32; N] {
        self.strides.into()
    }

    pub fn flatten_index(&self, idx: impl Borrow<[i32; N]>) -> i32 {
        let i = Vector::from(*idx.borrow());
        i.dot(&self.strides)
    }
}

impl<const N: usize, F> GridInfo<N, F>
where
    F: Float + NumAssignOps + AsPrimitive<i32> + std::fmt::Debug,
{
    pub fn new(aabb: Aabb<N, F>, cutoff: F) -> Self {
        let shape = ((aabb.sup - aabb.inf) / cutoff).map(|coord| coord.floor().as_() + 1);

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
            //FIXME: "attempt to multiply with overflow"
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
    //TODO: not sure where it fits better
    //TODO: GridInfo knows enough to compute the cell index for an arbitrary point
    //TODO: but might make more sense in FlatIndex?
    //TODO: sth. like Lattice trait maybe
    pub fn cell_index(&self, point: impl Borrow<[F; N]>) -> [i32; N] {
        let point = Point::from(*point.borrow());

        let idx = ((point - self.aabb.inf) / self.cutoff).map(|coord| coord.floor().as_());

        // FIXME: cell index ([2, -1, 0]) out of bounds ([1, 1, 1])
        assert!(
            idx < self.shape,
            "cell index ({:?}) out of bounds ({:?})",
            idx,
            self.shape
        );

        idx.into()
    }

    pub fn flat_cell_index(&self, point: impl Borrow<[F; N]>) -> i32 {
        let point = Point::from(*point.borrow());

        ((point - self.aabb.inf) / self.cutoff)
            .map(|coord| coord.floor().as_())
            .dot(&self.strides)

        //TODO: bounds checks are a bit unintuitive now. might defer it to hashmap lookup?
        // note that the following line is not as efficient:
        // self.flatten_index(self.cell_index(point))
    }
}

impl<const N: usize, F> Default for GridInfo<N, F>
where
    F: Float + std::fmt::Debug + Default + NumAssignOps + AsPrimitive<i32> + ConstOne,
{
    fn default() -> Self {
        GridInfo::new(Aabb::default(), F::ONE)
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
                    points.push([
                        cutoff.mul_add(x as f64, origin[0]),
                        cutoff.mul_add(y as f64, origin[1]),
                        cutoff.mul_add(z as f64, origin[2]),
                    ]);
                    points.push([
                        cutoff.mul_add(x as f64, cutoff.mul_add(0.5, origin[0])),
                        cutoff.mul_add(y as f64, cutoff.mul_add(0.5, origin[1])),
                        cutoff.mul_add(z as f64, cutoff.mul_add(0.5, origin[2])),
                    ]);
                }
            }
        }
    }

    points
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_points() {
        let points = vec![
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5],
            [0.0, 0.0, 2.0],
            [0.5, 0.5, 2.5],
            [0.0, 1.0, 1.0],
            [0.5, 1.5, 1.5],
            [0.0, 2.0, 0.0],
            [0.5, 2.5, 0.5],
            [0.0, 2.0, 2.0],
            [0.5, 2.5, 2.5],
            [1.0, 0.0, 1.0],
            [1.5, 0.5, 1.5],
            [1.0, 1.0, 0.0],
            [1.5, 1.5, 0.5],
            [1.0, 1.0, 2.0],
            [1.5, 1.5, 2.5],
            [1.0, 2.0, 1.0],
            [1.5, 2.5, 1.5],
            [2.0, 0.0, 0.0],
            [2.5, 0.5, 0.5],
            [2.0, 0.0, 2.0],
            [2.5, 0.5, 2.5],
            [2.0, 1.0, 1.0],
            [2.5, 1.5, 1.5],
            [2.0, 2.0, 0.0],
            [2.5, 2.5, 0.5],
            [2.0, 2.0, 2.0],
            [2.5, 2.5, 2.5],
        ];
        assert_eq!(points, generate_points([3, 3, 3], 1.0, [0.0, 0.0, 0.0]));
    }

    #[test]
    fn test_utils() {
        let points = generate_points([3, 3, 3], 1.0, [0.2, 0.25, 0.3]);
        assert_eq!(points.len(), 28, "testing PointCloud.len()");

        let aabb = Aabb::from_points::<[_; 3]>(points.iter());
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
            [0.2, 0.25, 0.3],
            "testing GridInfo.origin()"
        );
        assert_eq!(grid_info.shape(), [3, 3, 3], "testing GridInfo.shape");
        //TODO: note that these are the strides for grid_info.shape + [1, 1, 1]
        //TODO: this allows us to use negative relative indices representable by BalancedTernary<N>
        assert_eq!(grid_info.strides(), [1, 4, 16], "testing GridInfo.strides");

        // Intuitively you'd expect [2, 2, 2] for this
        // but we're having floating point imprecision:
        // 2.3 - 0.3 = 1.9999999999999998
        // This is not an issue though because the index is still uniquely determined
        assert_eq!(
            grid_info.cell_index(&[2.7, 2.75, 2.3]),
            [2, 2, 1],
            "testing GridInfo.cell_index()"
        );
        assert_eq!(
            grid_info.flat_cell_index(&[2.7, 2.75, 2.3]),
            26,
            "testing GridInfo.flat_cell_index()"
        );
        assert_eq!(
            grid_info.cell_index(&[2.7, 2.75, 2.8]),
            [2, 2, 2],
            "testing GridInfo.cell_index()"
        );
        assert_eq!(
            grid_info.flat_cell_index(&[2.7, 2.75, 2.8]),
            42,
            "testing GridInfo.flat_cell_index()"
        );
    }
}
