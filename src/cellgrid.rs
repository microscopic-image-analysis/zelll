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
//TODO: decide what's getting re-exported at the crate root
#[derive(Debug, Default, Clone)]
pub struct CellGrid<P, const N: usize = 3, T: Float = f64>
where
    T: NumAssignOps + ConstOne + AsPrimitive<i32> + std::fmt::Debug,
{
    cells: HashMap<i32, CellSliceMeta>,
    cell_lists: CellStorage<(usize, P)>,
    index: FlatIndex<N, T>,
}

// TODO: maybe should provide ::rebuild_from_internal()? that feeds cellgrid with iterator over internal CellStorage
// TODO: then we should also provide point_pairs_mut() that directly allows to overwrite points in cell storage?
// TODO: although that complicates leapfrog integration
// TODO: instead provide convenience method to chain GridCell::iter()'s or directly iterate over CellStorage
// TODO: however for sequential data (biomolecules) this may destroy implicit order (so that would have to be modelled implicitley)
// TODO: also, with reordered data, leapfrog integration buffers have to be shuffled accordingly
// TODO: which is not that nice
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
            let mut cells: HashMap<i32, CellSliceMeta> = HashMap::new();
            index.index.iter().for_each(|idx| {
                // FIXME: this will trigger debug_assert!() in CellSliceMeta::move_cursor()
                cells.entry(*idx).or_default().move_cursor(1);
            });

            // TODO: could also use Itertools::next_array() here
            // TODO: could count in chunk using Itertools::counts()
            // TODO: but then I'd have to put those into an array of chunk size (with purposefully invalid/non-existing keys)
            // let mut sliced = &index.index[..];
            // while let Some((&chunk, tail)) = sliced.split_first_chunk::<8>() {
            // for idx in chunk[..].iter() {
            //     cells.entry(*idx).or_default().move_cursor(1);
            // }

            //     use itertools::Itertools;
            //     let occurrences = chunk.iter().counts();
            //     occurrences.into_iter().for_each(|(idx, count)| {
            //         cells.entry(*idx).or_default().move_cursor(count);
            //     });

            //     sliced = tail;
            // }

            // sliced.iter().for_each(|idx| {
            //     cells.entry(*idx).or_default().move_cursor(1);
            // });

            // use itertools::Itertools;
            // index.index.iter().chunks(8).into_iter().for_each(|chunk| {
            //     let occurrences = chunk.counts();

            //     for (idx, count) in occurrences {
            //         cells.entry(*idx).or_default().move_cursor(count);
            //     }
            // });

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
                    // FIXME: in principle could have multiple &mut slices into CellStorage (for parallel pushing)
                    // FIXME: would just have to make sure that cell is always unique when operating on chunks
                    // FIXME: (pretty much the same issue as with counting cell sizes concurrently)
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
    // TODO: rebuild() could simply do this but rebuild_mut() on
    // TODO: an empty CellGrid does have some overhead due to FlatIndex::rebuild_mut()
    // {
    //     let mut cellgrid = self;
    //     cellgrid.rebuild_mut(points, cutoff);
    //     cellgrid
    // }

    // TODO: documentation
    // TODO: currently, rebuild_mut() will usually be more expensive than rebuild() even though it does not
    // TODO: allocate.
    // TODO: it makes mostly sense to use this when it's unlikely that the index changed
    // TODO: when a simulation reaches a stable state or cutoff is larger than the max. interaction distance
    // TODO: (cutoff-interaction distance >= max. distance particles might travel in one simulation step (similar to verlet lists))
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
                // FIXME: this will trigger debug_assert!() in CellSliceMeta::move_cursor()
                self.cells.entry(*idx).or_default().move_cursor(1);
            });

            // let mut sliced = &self.index.index[..];
            // while let Some((&chunk, tail)) = sliced.split_first_chunk::<4>() {
            //     self.cells.entry(chunk[0]).or_default().move_cursor(1);
            //     self.cells.entry(chunk[1]).or_default().move_cursor(1);
            //     self.cells.entry(chunk[2]).or_default().move_cursor(1);
            //     self.cells.entry(chunk[3]).or_default().move_cursor(1);
            //     sliced = tail;
            // }

            // sliced.iter().for_each(|idx| {
            //     self.cells.entry(*idx).or_default().move_cursor(1);
            // });

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
    /// ```
    /// # use zelll::cellgrid::CellGrid;
    /// use nalgebra::distance_squared;
    /// # let points = [[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
    /// # let cell_grid = CellGrid::new(points.iter().copied(), 1.0);
    /// cell_grid.point_pairs()
    ///     // usually, .filter_map() is preferable (so distance computations can be re-used)
    ///     .filter(|&((_i, p), (_j, q))| {
    ///         distance_squared(&p.into(), &q.into()) <= 1.0
    ///     })
    ///     .for_each(|((_i, p), (_j, q))| {
    ///         /* do some work */
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
    /// TODO: fact-check the statement below:
    /// Try to avoid filtering this [`ParallelIterator`] to avoid significant overhead:
    /// ```
    /// # // TODO: re-export ParallelIterator
    /// # use crate::zelll::cellgrid::ParallelIterator;
    /// # use zelll::cellgrid::CellGrid;
    /// use nalgebra::distance_squared;
    /// # let points = [[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
    /// # let cell_grid = CellGrid::new(points.iter().copied(), 1.0);
    /// cell_grid.par_point_pairs()
    ///     .for_each(|((_i, p), (_j, q))| {
    ///         if distance_squared(&p.into(), &q.into()) <= 1.0 {
    ///             /* do some work */
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
            .flat_map_iter(|cell| cell.point_pairs())
            .map(|((i, _p), (j, _q))| (i, j))
    }
}
