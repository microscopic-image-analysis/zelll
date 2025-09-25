//! Several utility items that might be useful but usually do not need be interacted with.
use crate::Particle;
use nalgebra::{Point, SVector, SimdPartialOrd};
use num_traits::{AsPrimitive, ConstOne, ConstZero, Float, NumAssignOps};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::borrow::Borrow;

/// A `Vec<[f64; N]>`, aliased for convenience and testing purposes.
pub type PointCloud<const N: usize> = Vec<[f64; N]>;

/// An axis-aligned bounding box constructed from particle data and used internally.
///
/// An `Aabb` is conveniently expressed by two points, infimum and supremum,
/// where each particle datum is partially ordered between those points.
///
/// See also [`GridInfo::bounding_box()`].
//TODO: rename fields, infimum/supremum might be confusing outside of a lattice context
#[derive(Clone, Copy, Debug, PartialEq, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Aabb<const N: usize = 3, F: Float = f64>
where
    F: std::fmt::Debug + 'static,
{
    inf: Point<F, N>,
    sup: Point<F, N>,
}

impl<const N: usize, F> Aabb<N, F>
where
    F: Float + std::fmt::Debug + SimdPartialOrd + ConstZero,
{
    /// Computes the componentwise minimum and maximum from the coordinates
    /// of the supplied particle data.
    pub fn from_particles<P: Particle<[F; N]>>(
        mut particles: impl Iterator<Item = impl Borrow<P>>,
    ) -> Self {
        let init = particles
            .next()
            .map(|p| p.borrow().coords())
            .unwrap_or([F::ZERO; N]);
        let init = Point::from(init);

        let (inf, sup) = particles
            .take(i32::MAX as usize)
            .fold((init, init), |(i, s), p| {
                let p = Point::from(p.borrow().coords());
                (i.inf(&p), s.sup(&p))
            });

        Self { inf, sup }
    }

    //TODO: could also pass iterators here (single point could be wrapped by std::iter::once or Option::iter())
    fn update<P: Particle<[F; N]>>(&mut self, particle: impl Borrow<P>) {
        let p = Point::from(particle.borrow().coords());
        self.inf = p.inf(&self.inf);
        self.sup = p.sup(&self.sup);
    }

    /// Returns a copy of the infimum of this `Aabb`.
    pub fn inf(&self) -> [F; N] {
        self.inf.into()
    }

    /// Returns a copy of the supremum of this `Aabb`.
    pub fn sup(&self) -> [F; N] {
        self.sup.into()
    }
}

/// A type containing geometry information necessary to compute cell indices
/// of a [`CellGrid`](crate::CellGrid).
///
/// <div class="warning">
///
/// The grid described by `GridInfo` may be slightly larger than the underlying bounding box.
///
/// </div>
#[derive(Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct GridInfo<const N: usize = 3, F: Float = f64>
where
    F: std::fmt::Debug + 'static,
{
    pub(crate) aabb: Aabb<N, F>,
    pub(crate) cutoff: F,
    shape: SVector<i32, N>,
    strides: SVector<i32, N>,
}

// impl<const N: usize, F> Default for Aabb<N, F>
// where
//     F: Float + std::fmt::Debug,
// {
//     fn default() -> Self {
//         Self {
//             inf: Point::from([F::min_value(); N]),
//             sup: Point::from([F::max_value(); N]),
//         }
//     }
// }

// TODO: Experiment with a fixed shape/strides
// TODO: technically, we don't even need a bounding box, arbitrarily big strides would suffice
// TODO: not quite like this:
// TODO: https://matthias-research.github.io/
// TODO: which would require re-computing hashes to get neighbor cells, instead of just adding
// TODO: however we can just arbitrarily big prime strides (and use ADD instead of XOR)
// TODO: as long as point clouds fit into that shape
// FIXME: this does improve cell grid construction/rebuild but penalizes particle_pairs()
// impl<const N: usize, F> GridInfo<N, F>
// where
//     F: Float + NumAssignOps + AsPrimitive<i32> + std::fmt::Debug + Default,
// {
//     pub(crate) fn new_fixed(cutoff: F) -> Self {
//         let shape = SVector::from([65536; N]);
//         let mut strides = shape;
//         strides.iter_mut().fold(1, |prev, curr| {
//             let next = prev * (*curr + 1);
//             *curr = prev;
//             next
//         });
//
//         Self {
//             aabb: Aabb::default(),
//             cutoff,
//             shape,
//             strides,
//         }
//     }
// }

impl<const N: usize, F> GridInfo<N, F>
where
    F: Float + std::fmt::Debug,
{
    /// Returns the origin corner of this grid.
    pub fn origin(&self) -> [F; N] {
        self.aabb.inf.into()
    }

    /// Returns the shape this grid is partitioned into.
    pub fn shape(&self) -> [i32; N] {
        self.shape.into()
    }

    /// Returns the strides this grid is using to compute flat cell indices.
    pub fn strides(&self) -> [i32; N] {
        self.strides.into()
    }

    /// Returns a reference to the internal axis-aligned bounding box.\
    /// See [`Aabb`] for details.
    pub fn bounding_box(&self) -> &Aabb<N, F> {
        &self.aabb
    }

    /// Computes a flat index from cell "coordinates".
    ///
    /// <div class="warning">
    ///
    /// Although negative coordinate values are technically allowed, only values `>= -1i32`
    /// are actually valid. We need this to compute cell neighborhoods easily.
    ///
    /// </div>
    ///
    /// # Panics
    ///
    /// Panics on **debug** builds if at least one component of `idx` is strictly less than `-1`.
    pub fn flatten_index(&self, idx: impl Borrow<[i32; N]>) -> i32 {
        debug_assert!(*idx.borrow() >= [-1i32; N]);

        let i = SVector::from(*idx.borrow());
        i.dot(&self.strides)
    }

    /// Returns the cutoff radius (ie. the cell edge length) of this grid.
    #[inline]
    pub fn cutoff(&self) -> F {
        self.cutoff
    }
}

impl<const N: usize, F> GridInfo<N, F>
where
    F: Float + NumAssignOps + AsPrimitive<i32> + std::fmt::Debug,
{
    /// Constructs a new instance of `GridInfo`.
    /// Usually, this does not need to be used manually.
    pub fn new(aabb: Aabb<N, F>, cutoff: F) -> Self {
        // +1 to cover all cells fitting inside the bounding box `aabb`
        // FIXME: do I really need this? I'll add two-layer padding to the strides and
        // FIXME: in try_cell_index() make sure to return None if the first layer is exceeded
        // // +2 to explicitly allow one layer of empty cells around the bounding box
        // // (so they can be queried but we won't store those in our hash map)
        // // so in total +3
        let shape = ((aabb.sup - aabb.inf) / cutoff).map(|coord| coord.floor().as_() + 1);

        let mut strides = shape;
        strides.iter_mut().fold(1, |prev, curr| {
            // let next = prev * (*curr + 1);
            // +1 would suffice to allow one layer of negative relative flat indices
            // +2 allows to query one layer outside all around `shape`
            // +4 adds two layers of implicit padding around `shape`
            // such that querying cells in the first padding layer
            // produces _unique_ neighbor flat indices.
            // This is important because we would get helical boundaries otherwise
            let next = prev * (*curr + 4);
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

    /// Computes integer "coordinates" of the cell the given coordinates belong to,
    /// ie. a _cell index_.
    ///
    /// # Panics
    ///
    /// Panics if the computed cell index does not fit the shape of this grid.
    /// See [`try_cell_index()`](GridInfo::try_cell_index()) for details.
    pub fn cell_index(&self, coordinates: impl Borrow<[F; N]>) -> [i32; N] {
        self.try_cell_index(coordinates)
            .expect("cell index is out of bounds")
    }

    /// Computes integer "coordinates" of the cell the given coordinates belong to,
    /// ie. a _cell index_.
    /// Returns `None` if the computed index is _too far_ outside of this `CellGrid`.
    ///
    /// <div class="warning">
    ///
    /// Querying locations within a boundary of thickness `cutoff` around the shape of
    /// the cell grid is allowed because their neighborhoods might reach into the
    /// grid's bounding box.
    ///
    /// </div>
    pub fn try_cell_index(&self, coordinates: impl Borrow<[F; N]>) -> Option<[i32; N]> {
        let p = Point::from(*coordinates.borrow());
        let idx = ((p - self.aabb.inf) / self.cutoff).map(|coord| coord.floor().as_());

        // cf. GridInfo::new()
        // we allow querying indices in the first implicit padding layer
        if SVector::from([-1; N]) <= idx && idx <= self.shape {
            Some(idx.into())
        } else {
            None
        }
    }

    /// Computes a flat cell index directly.
    ///
    /// <div class="warning">
    ///
    /// In contrast to [`flatten_index()`](GridInfo::flatten_index()) and
    /// [`cell_index()`](GridInfo::cell_index()), this does **not** panic although
    /// the same invariants are expected to hold:
    ///
    /// </div>
    ///
    /// ```
    /// # use zelll::CellGrid;
    /// # let data = vec![[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
    /// let cell_grid = CellGrid::new(data.iter().copied(), 1.0);
    /// # let info = cell_grid.info();
    /// let p = [-1.0; 3];
    /// assert_eq!(info.flat_cell_index(p), info.flatten_index(info.cell_index(p)));
    /// ```
    /// ```should_panic
    /// # use zelll::CellGrid;
    /// # let data = vec![[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
    /// # let cell_grid = CellGrid::new(data.iter().copied(), 1.0);
    /// # let info = cell_grid.info();
    /// let p = [-2.0; 3];
    /// // this is fine
    /// info.flat_cell_index(p);
    /// // this will panic
    /// info.cell_index(p);
    /// ```
    // the reason this does not do any bounds checks is twofold:
    // 1. particles used for constructing a specific cell grid always compute to a valid index
    // 2. if the resulting flat cell index is invalid (ie. empty cell), its lookup in our internal hash map
    //    will fail, which is fine
    pub fn flat_cell_index(&self, coordinates: impl Borrow<[F; N]>) -> i32 {
        let p = Point::from(*coordinates.borrow());

        ((p - self.aabb.inf) / self.cutoff)
            .map(|coord| coord.floor().as_())
            .dot(&self.strides)
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

/// Generates a [`Vec`] of 3-dimensional points for testing purposes.
///
/// In a grid of `shape` with cells of length `cutoff`, only cells with even linear index contain points (chessboard pattern).
/// These non-empty cells contain two points each:
///
/// - the first at the origin of the cell (equivalent to the cell's multi-index + the `origin` of the grid)
/// - the second at the center of the cell
// We'll stay in 3D for simplicity here
pub fn generate_pointcloud(shape: [usize; 3], cutoff: f64, origin: [f64; 3]) -> PointCloud<3> {
    let mut pointcloud = Vec::with_capacity(shape.iter().product::<usize>().div_ceil(2));

    for x in 0..shape[0] {
        for y in 0..shape[1] {
            for z in 0..shape[2] {
                if (x + y + z) % 2 == 0 {
                    pointcloud.push([
                        cutoff.mul_add(x as f64, origin[0]),
                        cutoff.mul_add(y as f64, origin[1]),
                        cutoff.mul_add(z as f64, origin[2]),
                    ]);
                    pointcloud.push([
                        cutoff.mul_add(x as f64, cutoff.mul_add(0.5, origin[0])),
                        cutoff.mul_add(y as f64, cutoff.mul_add(0.5, origin[1])),
                        cutoff.mul_add(z as f64, cutoff.mul_add(0.5, origin[2])),
                    ]);
                }
            }
        }
    }

    pointcloud
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_generate_pointcloud() {
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
        assert_eq!(points, generate_pointcloud([3, 3, 3], 1.0, [0.0, 0.0, 0.0]));
    }

    #[test]
    fn test_utils() {
        let points = generate_pointcloud([3, 3, 3], 1.0, [0.2, 0.25, 0.3]);
        assert_eq!(points.len(), 28, "testing PointCloud.len()");

        let aabb = Aabb::from_particles::<[_; 3]>(points.iter());
        assert_eq!(
            aabb,
            Aabb {
                inf: [0.2, 0.25, 0.3].into(),
                sup: [2.7, 2.75, 2.8].into()
            },
            "testing Aabb::from_particles()"
        );

        let grid_info = GridInfo::new(aabb, 1.0);
        assert_eq!(
            grid_info.origin(),
            [0.2, 0.25, 0.3],
            "testing GridInfo.origin()"
        );
        assert_eq!(grid_info.shape(), [3, 3, 3], "testing GridInfo.shape");
        //TODO: note that these are the strides for grid_info.shape + [4, 4, 4]
        assert_eq!(grid_info.strides(), [1, 7, 49], "testing GridInfo.strides");

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
            65,
            "testing GridInfo.flat_cell_index()"
        );
        assert_eq!(
            grid_info.cell_index(&[2.7, 2.75, 2.8]),
            [2, 2, 2],
            "testing GridInfo.cell_index()"
        );
        assert_eq!(
            grid_info.flat_cell_index(&[2.7, 2.75, 2.8]),
            114,
            "testing GridInfo.flat_cell_index()"
        );
    }
}
