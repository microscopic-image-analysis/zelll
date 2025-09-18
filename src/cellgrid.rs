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

use crate::Particle;
#[cfg(feature = "rayon")]
use crate::rayon::ParallelIterator;
use flatindex::FlatIndex;
use hashbrown::HashMap;
#[doc(inline)]
pub use iters::GridCell;
pub use iters::neighborhood;
use nalgebra::SimdPartialOrd;
use num_traits::{AsPrimitive, ConstOne, ConstZero, Float, NumAssignOps};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use storage::{CellSliceMeta, CellStorage};
#[doc(inline)]
pub use util::{Aabb, GridInfo};

/// The central type representing a grid of cells that provides an implementation of the _cell lists_ algorithm.
///
/// While iterating over cells, particles and particle pairs is the most important functionality of `CellGrid`,
/// there is not much that could go wrong using the appropriate methods of this type.
///
/// In contrast, there are some intricacies with setting up the necessary data structures internal to `CellGrid,
/// despite it being mostly straight-forward.
/// There are many variations of this algorithm, ours is a relatively simple one:
///
/// # Algorithm sketch
/// 1. compute (axis-aligned) bounding box
/// 2. compute cell index _`i`_ for each particle
/// 3. pre-allocate storage buffer for _`n`_ particles
/// 4. count number of particles _`nᵢ`_ in cell _`i`_
/// 5. slice storage buffer according to cell sizes
/// 6. copy particles into cell slices and store cell slice information in hash map
///
/// <div class="warning">
///
/// Note that step 1 is not strictly necessary, any bounding box will do,
/// so we could just use the largest possible one and save some time.
/// However, we can afford this.\
/// Also, there are many ways to do step 2.
/// </div>
///
/// Essentially, our `CellGrid` construction is an instance of
/// [_counting sort_](https://en.wikipedia.org/wiki/Counting_sort)
/// with the number of buckets _`k`_ (here: _non-empty_ cells) being bounded by _`n`_
/// due to using sparse storage.
///
/// Because we are sorting data, construction of a cell grid is _cache aware_.
/// However, iterating over particle pairs after the fact benefits from this.
/// Unfortunately, there is not much we can do about this.
///
/// In some settings, input data already has some structure reducing this effect.
/// Pre-sorting data can be helpful but is not always an option, as is the case for some
/// types of simulation.
///
/// # Examples
///
/// The [prototypical example](crate#examples) does not require explicit type annotations:
///
/// ```
/// # use zelll::CellGrid;
/// let data = vec![[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
/// let mut cg = CellGrid::new(data.iter().copied(), 1.0);
/// ```
/// Equivalently:
/// ```
/// # use zelll::CellGrid;
/// # let data = vec![[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
/// let mut cg: CellGrid<[f64; 3], 3, f64> = CellGrid::new(data.iter().copied(), 1.0);
/// // usually, type parameters can be elided:
/// let mut cg: CellGrid<_, 3, f64> = CellGrid::new(data.iter().copied(), 1.0);
/// ```
/// Often, `f32` coordinates are sufficient for simulations:
/// ```
/// # use zelll::CellGrid;
/// let data = vec![[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
/// let mut cg: CellGrid<_, 3, f32> = CellGrid::new(data.iter().copied(), 1.0);
/// ```
/// If you do not need to spell out the specific type of the cell grid, the type checker is happy
/// with you specifying the input data explicitly:
/// ```
/// # use zelll::CellGrid;
/// let data: Vec<[f32; 3]> = vec![[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
/// let mut cg = CellGrid::new(data.iter().copied(), 1.0);
/// // similarly:
/// let data: Vec<[f32; 2]> = vec![[0.0, 0.0], [1.0,2.0], [0.0, 0.1]];
/// let mut cg = CellGrid::new(data.iter().copied(), 1.0);
/// ```
/// Any type implementing [`Particle`] can be used:
/// ```
/// # use zelll::CellGrid;
/// use nalgebra::SVector;
///
/// let data: Vec<SVector<f32, 2>> = vec![[0.0, 0.0].into(), [1.0,2.0].into(), [0.0, 0.1].into()];
/// let mut cg = CellGrid::new(data.iter().copied(), 1.0);
/// ```
#[derive(Debug, Default, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
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
    /// Constructs a new `CellGrid` from particle data and a cutoff distance.
    ///
    /// # Examples
    /// ```
    /// # use zelll::CellGrid;
    /// let data = vec![[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
    /// let cell_grid = CellGrid::new(data.iter().copied(), 1.0);
    /// ```
    ///
    /// <div class="warning">
    ///
    /// Note that the intended usage includes `.iter().copied()`
    /// because we require the items of `particles` to implement `Particle` and we do not
    /// provide a blanket implementation for references to types implementing `Particle`.
    ///
    /// Previously we required `<I as IntoIterator>::Item: Borrow<P>`
    /// which provided more flexibility by accepting both `P` and `&P`
    /// but this forces lots of type annotations on the user.\
    /// Since `Particle<T>: Copy` anyway, we sacrifice this flexibility
    /// in favor of not cluttering user code with type annotations.
    ///
    /// </div>
    pub fn new<I>(particles: I, cutoff: T) -> Self
    where
        I: IntoIterator<Item = P> + Clone,
        P: Default,
    {
        CellGrid::default().rebuild(particles, Some(cutoff))
    }

    /// Consumes `self` and rebuilds the cell grid from the supplied iterator.
    /// This method can be used to add or remove particles but a full rebuild will happen anyway,
    /// with the exception of every particle still belonging to the same cell as before.
    ///
    /// # Examples
    /// ```
    /// # use zelll::CellGrid;
    /// # let data = vec![[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
    /// # let mut cell_grid = CellGrid::new(data.iter().copied(), 1.0);
    /// // rebuild `cell_grid` with reversed input data and unchanged cutoff distance
    /// cell_grid = cell_grid.rebuild(data.iter().rev().copied(), None);
    /// ```
    #[must_use = "rebuild() consumes `self` and returns the rebuilt `CellGrid`"]
    pub fn rebuild<I>(self, particles: I, cutoff: Option<T>) -> Self
    where
        I: IntoIterator<Item = P> + Clone,
        P: Default,
    {
        let cutoff = cutoff.unwrap_or(self.index.grid_info.cutoff);
        let index = FlatIndex::from_particles(particles.clone(), cutoff);

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

            // technically, this is O(capacity) not O(n)
            // cf. https://github.com/rust-lang/rust/pull/97215
            cells.iter_mut().for_each(|(_, slice)| {
                *slice = cell_lists.reserve_cell(slice.cursor());
            });
            // FIXME: what happens (below) if I reserve cells sorted by their size (above)?
            // FIXME: this seems to be more likely to be the cache miss culprit
            // FIXME: can we do something clever here? use an LRU cache?
            // FIXME: use sth. like itertools::tree_reduce() to somehow deal with
            // FIXME: the random access pattern of cell_lists' slices?
            index
                .index
                .iter()
                .zip(particles)
                .enumerate()
                // TODO: we know cells.get_mut() won't fail by construction
                // TODO: but maybe use try_for_each() instead
                .for_each(|(i, (cell, particle))| {
                    // FIXME: in principle could have multiple &mut slices into CellStorage (for parallel pushing)
                    // FIXME: would just have to make sure that cell is always unique when operating on chunks
                    // FIXME: (pretty much the same issue as with counting cell sizes concurrently)
                    cell_lists.push(
                        (i, particle),
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
    // pub fn rebuild<I>(self, particles: I, cutoff: Option<T>) -> Self
    // where
    //     I: IntoIterator<Item = P> + Clone,
    //     P: Default,
    // {
    //     let mut cellgrid = self;
    //     cellgrid.rebuild_mut(particles, cutoff);
    //     cellgrid
    // }

    /// Rebuilds `self` but does so mutably. Internally allocated memory will be re-used but
    /// re-allocations may happen.
    /// In some settings, this method is preferred over [`rebuild()`](CellGrid::rebuild()),
    /// e.g. in simulations that are expected to converge towards a stable configuration.
    ///
    /// # Examples
    /// ```
    /// # use zelll::CellGrid;
    /// # let data = vec![[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
    /// # let mut cell_grid = CellGrid::new(data.iter().copied(), 1.0);
    /// // rebuild `cell_grid` with reversed input data and unchanged cutoff distance
    /// cell_grid.rebuild_mut(data.iter().rev().copied(), None);
    /// ```
    pub fn rebuild_mut<I>(&mut self, particles: I, cutoff: Option<T>)
    where
        I: IntoIterator<Item = P> + Clone,
        P: Default,
    {
        if self.index.rebuild_mut(particles.clone(), cutoff) {
            self.cells.clear();
            self.cell_lists.clear();

            // let estimated_cap = self.index.grid_info.shape().iter().product::<i32>().min(self.index.index.len() as i32);
            // self.cells.reserve(estimated_cap as usize);
            // self.cells.shrink_to(estimated_cap as usize);

            // FIXME: This should be HashMap<i32, Either<usize, CellSliceMeta>> or CellSliceMeta an enum
            self.index.index.iter().for_each(|idx| {
                // FIXME: see `storage.rs::CellStorage::move_cursor()` for details
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
                .zip(particles)
                .enumerate()
                //TODO: see `::rebuild()`
                .for_each(|(i, (cell, particle))| {
                    self.cell_lists.push(
                        (i, particle),
                        self.cells
                            .get_mut(cell)
                            .expect("cell grid should contain every cell in the grid index"),
                    )
                });
        }
    }
}

impl<P: Particle<[T; N]>, const N: usize, T> CellGrid<P, N, T>
where
    T: Float + ConstOne + AsPrimitive<i32> + std::fmt::Debug + NumAssignOps + Send + Sync,
    P: Send + Sync,
{
    /// Returns an iterator over all relevant (i.e. within cutoff threshold + some extra) unique
    /// pairs of particles in this `CellGrid`.
    ///
    /// <div class="warning">
    ///
    /// The `Item` type is currently `((usize, P), (usize, P))`,
    /// where each particle is labelled with its position it was inserted into this `CellGrid`.
    /// A future breaking change might remove this label (think _entity ID_), putting the responsibility to associate `P`
    /// with additional data (eg. velocities, momenta) onto implementors of [`Particle`].
    /// This will likely also deprecate [`pair_indices()`](CellGrid::pair_indices()).
    ///
    /// </div>
    ///
    /// # Examples
    /// ```
    /// # use zelll::CellGrid;
    /// use nalgebra::distance_squared;
    /// # let data = [[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
    /// # let cell_grid = CellGrid::new(data.iter().copied(), 1.0);
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
    pub fn particle_pairs(&self) -> impl Iterator<Item = ((usize, P), (usize, P))> + Clone {
        self.iter().flat_map(|cell| cell.particle_pairs())
    }

    /// Returns an iterator over all relevant (i.e. within cutoff threshold + some extra) unique
    /// pairs of particle indices in this `CellGrid`.
    ///
    /// A particle index is its position in the iterator that was used for constructing or rebuilding a `CellGrid`.
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn pair_indices(&self) -> impl Iterator<Item = (usize, usize)> + Clone {
        self.iter()
            .flat_map(|cell| cell.particle_pairs())
            .map(|((i, _p), (j, _q))| (i, j))
    }

    /// Returns spatial information about this cell grid, as well as auxiliary functionality
    /// facilitated by this information.
    ///
    /// See [`GridInfo`] for details.
    pub fn info(&self) -> &GridInfo<N, T> {
        &self.index.grid_info
    }

    /// Returns the [`GridCell`] the queried particle belongs into.
    /// The particle does not have to be present in the cell grid.
    ///
    /// <div class="warning">
    ///
    /// Note that the queried particle's cell must be within `cutoff` of this `CellGrid`'s bounding box.
    /// If that is not the case, `None` is returned.
    /// This restriction might be lifted in the future.
    ///
    /// </div>
    pub fn query<Q: Particle<[T; N]>>(&self, particle: Q) -> Option<GridCell<'_, P, N, T>> {
        self.info()
            .try_cell_index(particle.coords())
            .map(|index| self.info().flatten_index(index))
            .map(|index| GridCell { grid: self, index })
    }

    /// Returns an iterator over all relevant (i.e. within cutoff threshold + some extra)
    /// neighbor particles of the queried particle.
    /// This may include `particle` itself if it is part of this `CellGrid`.
    ///
    /// This is a convenience wrapper around [`query()`](CellGrid::query()).
    /// See also [`neighborhood`].
    ///
    /// ```
    /// # use zelll::CellGrid;
    /// use nalgebra::distance_squared;
    /// # let data = [[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
    /// # let cell_grid = CellGrid::new(data.iter().copied(), 1.0);
    /// let p = [0.5, 1.0, 0.1];
    /// cell_grid.query_neighbors(p)
    ///     .expect("the queried particle should be within `cutoff` of this grid's shape")
    ///     // usually, .filter_map() is preferable (so distance computations can be re-used)
    ///     .filter(|&(_j, q)| {
    ///         distance_squared(&p.into(), &q.into()) <= 1.0
    ///     })
    ///     .for_each(|(_j, q)| {
    ///         /* do some work */
    ///     });
    /// ```
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn query_neighbors<Q: Particle<[T; N]>>(
        &self,
        particle: Q,
    ) -> Option<impl Iterator<Item = (usize, P)> + Clone> {
        self.query(particle).map(|this| {
            this.iter().copied().chain(
                this.neighbors::<neighborhood::Full>()
                    .flat_map(|cell| cell.iter().copied()),
            )
        })
    }
}

#[cfg(feature = "rayon")]
impl<P, const N: usize, T> CellGrid<P, N, T>
where
    T: Float + NumAssignOps + ConstOne + AsPrimitive<i32> + Send + Sync + std::fmt::Debug,
    P: Particle<[T; N]> + Send + Sync,
{
    /// Returns a parallel iterator over all relevant (i.e. within cutoff threshold + some extra)
    /// unique pairs of particles in this `CellGrid`.
    ///
    /// <div class="warning">
    ///
    /// See also [`particle_pairs()`](CellGrid::particle_pairs()).
    ///
    /// </div>
    ///
    /// # Examples
    /// ```
    /// # use zelll::{CellGrid, rayon::ParallelIterator};
    /// use nalgebra::distance_squared;
    /// # let data = [[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
    /// # let cell_grid = CellGrid::new(data.iter().copied(), 1.0);
    /// cell_grid.par_particle_pairs()
    //      TODO: fact-check the statement below:
    ///     // Try to avoid filtering this ParallelIterator to avoid significant overhead:
    ///     .for_each(|((_i, p), (_j, q))| {
    ///         if distance_squared(&p.into(), &q.into()) <= 1.0 {
    ///             /* do some work */
    ///         }
    ///     });
    /// ```
    pub fn par_particle_pairs(&self) -> impl ParallelIterator<Item = ((usize, P), (usize, P))> {
        // TODO: ideally, we would schedule 2 threads for cell.particle_pairs() with the same CPU affinity
        // TODO: so they can share their resources
        self.par_iter().flat_map_iter(|cell| cell.particle_pairs())
    }

    /// Returns a parallel iterator over all relevant (i.e. within cutoff threshold + some extra) unique
    /// pairs of particle indices in this `CellGrid`.
    ///
    /// A particle index is its position in the iterator that was used for constructing or rebuilding a `CellGrid`.
    pub fn par_pair_indices(&self) -> impl ParallelIterator<Item = (usize, usize)> {
        self.par_iter()
            .flat_map_iter(|cell| cell.particle_pairs())
            .map(|((i, _p), (j, _q))| (i, j))
    }
}
