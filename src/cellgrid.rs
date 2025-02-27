//! Primary module of this crate.
//!
//! The most important items are available from the crate root.
//! Refer to items and submodules for more information.
#[allow(dead_code)]
mod flatindex;
#[allow(dead_code)]
mod iters;
#[allow(dead_code)]
mod storage;
#[allow(dead_code)]
pub mod util;

#[cfg(feature = "rayon")]
use crate::rayon::ParallelIterator;
use crate::Particle;
use flatindex::FlatIndex;
use hashbrown::HashMap;
#[doc(inline)]
pub use iters::GridCell;
use nalgebra::SimdPartialOrd;
use num_traits::{AsPrimitive, ConstOne, ConstZero, Float, NumAssignOps};
use storage::{CellSliceMeta, CellStorage};
#[doc(inline)]
pub use util::{Aabb, GridInfo};

/// The central type representing a grid of cells that provides an implementation of the _cell lists_ algorithm.
///
/// # TODO: Algorithm sketch
/// 1. compute (axis-aligned) bounding box
/// 2. compute cell index _`i`_ for each particle
/// 3. pre-allocate storage buffer for _`n`_ particles
/// 4. count number of particles _`nᵢ`_ in cell _`i`_
/// 5. slice storage buffer according to cell sizes
/// 6. copy particles into cell slices and store cell slice information in hash map
///
/// # Examples
/// TODO: examples with default/explicitly set type parameters (2D, 3D, f32/f64)\
/// TODO: examples with different types impl'ing `Particle`
#[derive(Debug, Default, Clone)]
pub struct CellGrid<P, const N: usize = 3, T: Float = f64>
where
    T: NumAssignOps + ConstOne + AsPrimitive<i32> + std::fmt::Debug,
{
    // TODO: experiment with hashbrown::HashTable
    // FIXME should expose type parameter S: BuildHasher publically
    cells: HashMap<i32, CellSliceMeta>,
    // TODO: rebuild from CellStorage iterator (boils down to counting/bucket sorting)
    // TODO: iterate (mutably) over cell storage, iterate mutably over particle pairs
    // TODO: make it responsibility of user to associate some index/label with P: Particle?
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

            // let estimated_cap = (index.grid_info.shape().iter().product::<i32>() as usize).min(index.index.len());
            // FIXME: This should be HashMap<i32, Either<usize, CellSliceMeta>> or CellSliceMeta an enum
            let mut cells: HashMap<i32, CellSliceMeta> = HashMap::new(); // HashMap::with_capacity(estimated_cap);
            index.index.iter().for_each(|idx| {
                // FIXME: this will trigger debug_assert!() in CellSliceMeta::move_cursor()
                cells.entry(*idx).or_default().move_cursor(1);
            });

            // cells.shrink_to_fit();

            // FIXME: is this why observed runtime does not match O(n) for large hash maps? (doesn't look like it is though)
            // FIXME: https://github.com/rust-lang/rust/pull/97215
            // FIXME: sorting benchmark data does not yield linear runtime (even though it massively improves cache miss rate)
            // FIXME: For rebuild() it's not that important but for rebuild_mut() it is!
            cells.iter_mut().for_each(|(_, slice)| {
                *slice = cell_lists.reserve_cell(slice.cursor());
            });
            // FIXME: what happens (below) if I reserve cells sorted by their size (above)?

            // FIXME: this seems to be more likely to be the cache miss culprit
            // FIXME: can we do something clever here? use an LRU cache?
            // FIXME: use sth. like itertools::tree_reduce() to somehow deal with
            // FIXME: the random access pattern of cell_lists' slices?
            // FIXME: maybe should not store cell indices in Vec but compute them *again* from points
            // FIXME: computation should be cheap enough and this way we might save some precious cache lines
            // FIXME: turns out it is in fact cheap enough (but does not improve cache behavior)
            // FIXME: so we should in fact remove FlatIndex (or make it a BTreeMap?)
            index
                .index
                .iter()
                .zip(points)
                .enumerate()
                // TODO: clean this up, this could be nicer since we know cells.get_mut() won't fail?
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
    // pub fn rebuild<I>(self, points: I, cutoff: Option<T>) -> Self
    // where
    //     I: IntoIterator<Item = P> + Clone,
    //     P: Default,
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

            // let estimated_cap = self.index.grid_info.shape().iter().product::<i32>().min(self.index.index.len() as i32);
            // self.cells.reserve(estimated_cap as usize);
            // self.cells.shrink_to(estimated_cap as usize);

            // FIXME: This should be HashMap<i32, Either<usize, CellSliceMeta>> or CellSliceMeta an enum
            self.index.index.iter().for_each(|idx| {
                // FIXME: this will trigger debug_assert!() in CellSliceMeta::move_cursor()
                self.cells.entry(*idx).or_default().move_cursor(1);
            });

            // TODO: Since hashmap iteration is `O(capacity)` not `O(len)` we want to make sure
            // TODO: that the load factor does not degenerate (resize policy says ~ 0.5-0.85)
            // TODO: however this means potential re-allocation
            self.cells.shrink_to_fit();

            self.cells.iter_mut().for_each(|(_, slice)| {
                *slice = self.cell_lists.reserve_cell(slice.cursor());
            });

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

    pub fn info(&self) -> &GridInfo<N, T> {
        &self.index.grid_info
    }
}

impl<P: Particle<[T; N]>, const N: usize, T> CellGrid<P, N, T>
where
    T: Float + ConstOne + AsPrimitive<i32> + std::fmt::Debug + NumAssignOps + Send + Sync,
    P: Send + Sync,
{
    /// Iterate over all relevant (i.e. within cutoff threshold + some extra) unique pairs of points in this `CellGrid`.
    ///
    /// # Examples
    /// ```
    /// # use zelll::CellGrid;
    /// use nalgebra::distance_squared;
    /// # let points = [[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
    /// # let cell_grid = CellGrid::new(points.iter().copied(), 1.0);
    /// cell_grid.particle_pairs()
    ///     // usually, .filter_map() is preferable (so distance computations can be re-used)
    ///     .filter(|&((_i, p), (_j, q))| {
    ///         distance_squared(&p.into(), &q.into()) <= 1.0
    ///     })
    ///     .for_each(|((_i, p), (_j, q))| {
    ///         /* do some work */
    ///     });
    /// ```
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn particle_pairs(&self) -> impl Iterator<Item = ((usize, P), (usize, P))> + Clone + '_ {
        self.iter().flat_map(|cell| cell.particle_pairs())
    }

    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn pair_indices(&self) -> impl Iterator<Item = (usize, usize)> + Clone + '_ {
        self.iter()
            .flat_map(|cell| cell.particle_pairs())
            .map(|((i, _p), (j, _q))| (i, j))
    }
}

#[cfg(feature = "rayon")]
impl<P, const N: usize, T> CellGrid<P, N, T>
where
    T: Float + NumAssignOps + ConstOne + AsPrimitive<i32> + Send + Sync + std::fmt::Debug,
    P: Particle<[T; N]> + Send + Sync,
{
    /// Iterate in parallel over all relevant (i.e. within cutoff threshold + some extra)
    /// unique pairs of points in this `CellGrid`.
    ///
    /// # Examples
    /// ```
    /// # // TODO: re-export ParallelIterator
    /// # use zelll::{CellGrid, rayon::ParallelIterator};
    /// use nalgebra::distance_squared;
    /// # let points = [[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
    /// # let cell_grid = CellGrid::new(points.iter().copied(), 1.0);
    /// cell_grid.par_particle_pairs()
    //      TODO: fact-check the statement below:
    ///     // Try to avoid filtering this ParallelIterator to avoid significant overhead:
    ///     .for_each(|((_i, p), (_j, q))| {
    ///         if distance_squared(&p.into(), &q.into()) <= 1.0 {
    ///             /* do some work */
    ///         }
    ///     });
    /// ```
    pub fn par_particle_pairs(
        &self,
    ) -> impl ParallelIterator<Item = ((usize, P), (usize, P))> + '_ {
        self.par_iter().flat_map_iter(|cell| cell.particle_pairs())
    }

    pub fn par_pair_indices(&self) -> impl ParallelIterator<Item = (usize, usize)> + '_ {
        self.par_iter()
            .flat_map_iter(|cell| cell.particle_pairs())
            .map(|((i, _p), (j, _q))| (i, j))
    }
}
