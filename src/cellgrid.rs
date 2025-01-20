//TODO: - allow both &[Point<f64, N>] and impl Iterator<Item = Point<f64, N>> (or IntoIterator)?
//TODO: - possible approach: require impl Iterator and then .take(usize::MAX) to ensure we have only finite point clouds
#[allow(dead_code)]
pub mod flatindex;
#[allow(dead_code)]
pub mod iters;
#[allow(dead_code)]
pub mod storage;
#[allow(dead_code)]
pub mod util;

use crate::Particle;
pub use flatindex::*;
use hashbrown::HashMap;
pub use iters::*;
use nalgebra::SimdPartialOrd;
use num_traits::{AsPrimitive, ConstOne, ConstZero, Float, NumAssignOps};
#[cfg(feature = "rayon")]
//TODO: should do a re-export of rayon?
pub use rayon::prelude::ParallelIterator;
pub use storage::*;
pub use util::*;
//TODO: make CellGrid and related stuff generic over (internally) used numeric types: https://docs.rs/num-traits/latest/num_traits/
//TODO: expose type aliases for practical combinations of numeric primitive types
//TODO: decide what's getting re-exported at the crate root
//TODO: I wonder if I could benefit from https://docs.rs/hashbrown/latest/hashbrown/struct.HashTable.html
//TODO: https://docs.rs/indexmap/latest/indexmap/ is essentially using HashTable but it has too much overhead for us
#[derive(Debug, Default, Clone)]
pub struct CellGrid<P, const N: usize = 3, T: Float = f64>
where
    T: NumAssignOps + ConstOne + AsPrimitive<i32> + std::fmt::Debug,
{
    cells: HashMap<i32, CellSliceMeta>,
    //cells: DenseMap,
    // TODO: could make this CellStorage<(usize, [f64; N])>
    // cell_lists: CellStorage<(usize, Point<T, N>)>,
    cell_lists: CellStorage<(usize, P)>,
    index: FlatIndex<N, T>,
}

// pub type ArrayParticleGrid<T, const N: usize> = CellGrid<[T; N], N, T>;

impl<P: Particle<[T; N]>, const N: usize, T> CellGrid<P, N, T>
where
    T: Float
        + NumAssignOps
        + ConstOne
        + ConstZero
        + AsPrimitive<i32>
        + SimdPartialOrd
        + Send
        + Sync
        + std::fmt::Debug
        + Default,
{
    // TODO: document that P: Copy is important here
    // TODO: intended usage: CellGrid::new(points.iter().copied())
    // TODO: previously we required <I as IntoIterator>::Item: Borrow<P>
    // TODO: which provided more flexibility (accepting P and &P)
    // TODO: but this required lots of type annotations for the user
    // TODO: since Particle<T>: Copy anyway, we'll sacrifices this flexibility
    // TODO: in favor of not cluttering user code with type annotations
    pub fn new<I>(points: I, cutoff: T) -> Self
    where
        I: IntoIterator<Item = P> + Clone,
        P: Default,
    {
        CellGrid::default().rebuild(points, Some(cutoff))
    }

    //TODO: Documentation: rebuild does not really update because it does not make a lot of sense:
    //TODO: If FlatIndex did change, we have to re-allocate (or re-initialize) almost everything anyway;
    //TODO: If FlatIndex did not change, we don't need to update.
    //TODO: Therefore we check for that and make CellGrid::new() just a wrapper around CellGrid::rebuild (with an initially empty FlatIndex)
    #[must_use = "rebuild() consumes `self` and returns the rebuilt `CellGrid`"]
    pub fn rebuild<I>(self, points: I, cutoff: Option<T>) -> Self
    where
        I: IntoIterator<Item = P> + Clone,
        P: Default,
    {
        let cutoff = cutoff.unwrap_or(self.index.grid_info.cutoff);
        let index = FlatIndex::from_points(points.clone(), cutoff);

        if index == self.index {
            self
        } else {
            let mut cell_lists = CellStorage::with_capacity(index.index.len());

            // FIXME: This should be HashMap<i32, Either<usize, CellSliceMeta>> or CellSliceMeta an enum
            let mut cells: HashMap<i32, CellSliceMeta> =
                index.index.iter().fold(HashMap::new(), |mut map, idx| {
                    map.entry(*idx).or_default().move_cursor(1);
                    map
                });

            cells.iter_mut().for_each(|(_, slice)| {
                *slice = cell_lists.reserve_cell(slice.cursor());
            });

            index
                .index
                .iter()
                .zip(points)
                .enumerate()
                //TODO: clean this up, this could be nicer since we know cells.get_mut() won't fail?
                .for_each(|(i, (cell, point))| {
                    cell_lists.push(
                        (i, point),
                        cells
                            .get_mut(cell)
                            .expect("cell grid should contain every cell in the grid index"),
                    )
                });

            Self {
                cells,
                cell_lists,
                index,
            }
        }
    }

    pub fn rebuild_mut<I>(&mut self, points: I, cutoff: Option<T>)
    where
        I: IntoIterator<Item = P> + Clone,
        P: Default,
    {
        if self.index.rebuild_mut(points.clone(), cutoff) {
            self.cells.clear();
            self.cell_lists.clear();

            // FIXME: This should be HashMap<i32, Either<usize, CellSliceMeta>> or CellSliceMeta an enum
            self.index.index.iter().for_each(|idx| {
                self.cells.entry(*idx).or_default().move_cursor(1);
            });

            self.cells.iter_mut().for_each(|(_, slice)| {
                *slice = self.cell_lists.reserve_cell(slice.cursor());
            });

            //TODO: we'll re-evaluate this once benchmarks with sparse data have been added
            //TODO: However, even without shrinking, self.cells.len() has an upper bound of  self.index.index.len()
            //TODO: (if updated point cloud did not shrink in length)
            self.cells.shrink_to_fit();

            // TODO: https://docs.rs/hashbrown/latest/hashbrown/struct.HashMap.html#method.get_many_mut
            // TODO: maybe could use get_many_mut here, but unfortunately we'd have to handle
            // TODO: duplicate keys (i.e. cells)
            // TODO: for that, we could chunk the index iter(), sort & count the chunks and then
            // TODO: get the cell once and push each particle index with a single lookup into the hashmap
            // TODO: perhaps this would be autovectorized?
            self.index
                .index
                .iter()
                .zip(points)
                .enumerate()
                //TODO: clean this up, this could be nicer since we know cells.get_mut() won't fail?
                //TODO: could use unwrap_unchecked() since this can't/shouldn't fail
                .for_each(|(i, (cell, point))| {
                    self.cell_lists.push(
                        (i, point),
                        self.cells
                            .get_mut(cell)
                            .expect("cell grid should contain every cell in the grid index"),
                    )
                });
        }
    }

    pub fn shape(&self) -> [i32; N] {
        self.index.grid_info.shape()
    }

    pub fn bounding_box(&self) -> &Aabb<N, T> {
        &self.index.grid_info.aabb
    }
}

impl<P: Particle<[T; N]>, const N: usize, T> CellGrid<P, N, T>
where
    T: Float + ConstOne + AsPrimitive<i32> + std::fmt::Debug + NumAssignOps + Send + Sync,
    P: Send + Sync,
{
    /// Iterate over all relevant (i.e. within cutoff threshold + some extra) unique pairs of points in this `CellGrid`.
    /// ```ignore
    /// cell_grid.point_pairs()
    ///     .filter(|&(i, j)| {
    ///         distance_squared(&points[i], &points[j]) <= cutoff_squared
    ///     })
    ///     .for_each(|&(i, j)| {
    ///         ...
    ///     });
    /// ```
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn point_pairs<'p>(
        &'p self,
    ) -> impl Iterator<Item = ((usize, P), (usize, P))> + Clone + 'p {
        self.iter().flat_map(|cell| cell.point_pairs())
    }

    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn pair_indices<'p>(&'p self) -> impl Iterator<Item = (usize, usize)> + Clone + 'p {
        self.iter()
            .flat_map(|cell| cell.point_pairs())
            .map(|((i, _p), (j, _q))| (i, j))
    }
}

#[cfg(feature = "rayon")]
impl<P, const N: usize, T> CellGrid<P, N, T>
where
    T: Float + NumAssignOps + ConstOne + AsPrimitive<i32> + Send + Sync + std::fmt::Debug,
    P: Particle<[T; N]> + Send + Sync,
{
    /// Iterate in parallel over all relevant (i.e. within cutoff threshold + some extra) unique pairs of points in this `CellGrid`.
    /// Try to avoid filtering this [`ParallelIterator`] to avoid significant overhead:
    /// ```ignore
    /// cell_grid.par_point_pairs()
    ///     .for_each(|(i, j)| {
    ///         if distance_squared(&pointcloud[i], &pointcloud[j]) <= cutoff_squared {
    ///             ...
    ///         } else {
    ///             ...
    ///         }
    ///     });
    /// ```
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn par_point_pairs(&self) -> impl ParallelIterator<Item = ((usize, P), (usize, P))> + '_ {
        self.par_iter().flat_map_iter(|cell| cell.point_pairs())
    }

    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn par_pair_indices(&self) -> impl ParallelIterator<Item = (usize, usize)> + '_ {
        self.par_iter()
            //.flat_map(|cell| cell.point_pairs().collect::<Vec<_>>())
            //.flat_map_iter(|cell| cell.point_pairs().collect::<Vec<_>>())
            .flat_map_iter(|cell| cell.point_pairs())
            .map(|((i, _p), (j, _q))| (i, j))
    }
}
