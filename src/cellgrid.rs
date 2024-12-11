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

pub use flatindex::*;
use hashbrown::HashMap;
pub use iters::*;
use nalgebra::SimdPartialOrd;
use num_traits::{AsPrimitive, ConstOne, ConstZero, Float, NumAssignOps};
#[cfg(feature = "rayon")]
//TODO: should do a re-export of rayon?
pub use rayon::prelude::ParallelIterator;
use std::borrow::Borrow;
pub use storage::*;
pub use util::*;
//TODO: make CellGrid and related stuff generic over (internally) used numeric types: https://docs.rs/num-traits/latest/num_traits/
//TODO: expose type aliases for practical combinations of numeric primitive types
//TODO: decide what's getting re-exported at the crate root
//TODO: I wonder if I could benefit from https://docs.rs/hashbrown/latest/hashbrown/struct.HashTable.html
//TODO: https://docs.rs/indexmap/latest/indexmap/ is essentially using HashTable but it has too much overhead for us
#[derive(Debug, Default, Clone)]
pub struct CellGrid<const N: usize = 3, T: Float = f64>
where
    T: NumAssignOps + ConstOne + AsPrimitive<i32> + std::fmt::Debug,
{
    cells: HashMap<i32, CellSliceMeta>,
    //cells: DenseMap,
    // TODO: could make this CellStorage<(usize, [f64; N])>
    cell_lists: CellStorage<usize>,
    index: FlatIndex<N, T>,
}

impl<const N: usize, T: Float + Send + Sync> CellGrid<N, T>
where
    T: NumAssignOps
        + ConstOne
        + ConstZero
        + AsPrimitive<i32>
        + SimdPartialOrd
        + std::fmt::Debug
        + Default,
{
    // TODO: whereever I need impl Iterator<>
    // TODO: I should probably use impl IntoIterator<Item = &'p [f64; N]> (+ Clone or + Borrow?)
    // TODO: usually this means, that the type impl'ing IntoIterator is a reference so it can be safely copied and iterated over
    // TODO: also maybe Item = impl Deref<Target = [f64; N]>? e.g. DashMap's iterators iterate over custom reference types
    // TODO: cf. https://users.rust-lang.org/t/declaring-associated-item-whose-references-implement-intoiterator/17103/6
    // TODO: cf. https://doc.rust-lang.org/stable/std/primitive.reference.html
    // TODO: cf. https://docs.rs/pyo3/latest/pyo3/prelude/trait.PyAnyMethods.html#implementors
    // TODO: impl IntoIterator for T: PyAnyMethods (kann ich das? brauch ich wrapper/super trait?)
    // TODO: Bound<'py, PyAny> impl's PyAnyMethods + Clone
    pub fn new(points: impl IntoIterator<Item = impl Borrow<[T; N]>> + Clone, cutoff: T) -> Self {
        CellGrid::default().rebuild(points, Some(cutoff))
    }

    //TODO: Documentation: rebuild does not really update because it does not make a lot of sense:
    //TODO: If FlatIndex did change, we have to re-allocate (or re-initialize) almost everything anyway;
    //TODO: If FlatIndex did not change, we don't need to update.
    //TODO: Therefore we check for that and make CellGrid::new() just a wrapper around CellGrid::rebuild (with an initially empty FlatIndex)
    #[must_use = "rebuild() consumes `self` and returns the rebuilt `CellGrid`"]
    pub fn rebuild(
        self,
        points: impl IntoIterator<Item = impl Borrow<[T; N]>> + Clone,
        cutoff: Option<T>,
    ) -> Self {
        let cutoff = cutoff.unwrap_or(self.index.grid_info.cutoff);
        let index = FlatIndex::from_points(points, cutoff);

        if index == self.index {
            self
        } else {
            let mut cell_lists = CellStorage::with_capacity(index.index.len());

            //TODO: there might be a better way for this temporary multiset/hashbag/histogram
            //TODO: e.g. have sth. like cells: HashMap<i32, Either<CellSliceMeta, usize>>
            //TODO: or divert CellSliceMeta from it's intended use and use CellSliceMeta::cursor as a counter
            //TODO: and initialize CellSliceMeta::range later once all counts are known
            let cell_sizes: HashMap<i32, usize> =
                index.index.iter().fold(HashMap::new(), |mut map, idx| {
                    *map.entry(*idx).or_default() += 1;
                    map
                });

            let mut cells: HashMap<i32, CellSliceMeta> = cell_sizes
                //let mut cells: DenseMap = cell_sizes
                .iter()
                .map(|(idx, count)| (*idx, cell_lists.reserve_cell(*count)))
                .collect();

            index
                .index
                .iter()
                .enumerate()
                //TODO: clean this up, this could be nicer since we know cells.get_mut() won't fail?
                .for_each(|(i, cell)| {
                    cell_lists.push(
                        i,
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

    pub fn rebuild_mut(
        &mut self,
        points: impl IntoIterator<Item = impl Borrow<[T; N]>> + Clone,
        cutoff: Option<T>,
    ) {
        if self.index.rebuild_mut(points, cutoff) {
            self.cells.clear();
            self.cell_lists.clear();

            let cell_sizes: HashMap<i32, usize> =
                self.index
                    .index
                    .iter()
                    .fold(HashMap::new(), |mut map, idx| {
                        *map.entry(*idx).or_default() += 1;
                        map
                    });

            cell_sizes.iter().for_each(|(idx, count)| {
                // SAFETY:
                // `insert_unique_unchecked()` is safe because `idx` is unique in `cell_sizes`
                // and `self.cells` has been cleared before.
                unsafe {
                    self.cells
                        .insert_unique_unchecked(*idx, self.cell_lists.reserve_cell(*count));
                }
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
                .enumerate()
                //TODO: clean this up, this could be nicer since we know cells.get_mut() won't fail?
                //TODO: could use unwrap_unchecked() since this can't/shouldn't fail
                .for_each(|(i, cell)| {
                    self.cell_lists.push(
                        i,
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
    pub fn point_pairs<'p>(&'p self) -> impl Iterator<Item = (usize, usize)> + Clone + 'p {
        self.iter().flat_map(|cell| cell.point_pairs())
    }

    /// Call a closure on each relevant (i.e. within cutoff threshold + some extra) unique pair of points in this `CellGrid`.
    /// `for_each_point_pair()` is equivalent to a nested loop:
    /// ```ignore
    /// for cell in cell_grid.iter() {
    ///     for (borb, other) in cell.point_pairs() {
    ///         ...
    ///     }
    /// }
    /// ```
    #[deprecated(note = "Please use [`point_pairs()`](#method.point_pairs) instead.")]
    pub fn for_each_point_pair<F>(&self, mut f: F)
    where
        F: FnMut(usize, usize),
    {
        self.iter().for_each(|cell| {
            cell.point_pairs().for_each(|(i, j)| {
                f(i, j);
            })
        });
    }

    #[deprecated(note = "Please use [`point_pairs()`](#method.point_pairs) instead.")]
    pub fn filter_point_pairs<F, P>(&self, mut f: F, p: P)
    where
        F: FnMut(usize, usize),
        P: Fn(usize, usize) -> bool,
    {
        //TODO: array_chunks() could be nice for autovectorization?
        //TODO: see https://doc.rust-lang.org/std/iter/trait.Iterator.html#method.array_chunks
        self.iter().for_each(|cell| {
            cell.point_pairs()
                .filter(|(i, j)| p(*i, *j))
                .for_each(|(i, j)| {
                    f(i, j);
                });
        });
    }
}

#[cfg(feature = "rayon")]
impl<const N: usize, T: Float + Send + Sync> CellGrid<N, T>
where
    T: NumAssignOps + ConstOne + AsPrimitive<i32> + std::fmt::Debug,
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
    pub fn par_point_pairs(&self) -> impl ParallelIterator<Item = (usize, usize)> + '_ {
        self.par_iter()
            //.flat_map(|cell| cell.point_pairs().collect::<Vec<_>>())
            //.flat_map_iter(|cell| cell.point_pairs().collect::<Vec<_>>())
            .flat_map_iter(|cell| cell.point_pairs())
    }

    #[deprecated(note = "Please use [`par_point_pairs()`](#method.par_point_pairs) instead.")]
    pub fn par_for_each_point_pair<F>(&self, f: F)
    where
        F: Fn(usize, usize) + Send + Sync,
    {
        self.par_iter().for_each(|cell| {
            cell.point_pairs().for_each(|(i, j)| {
                f(i, j);
            })
        });
    }

    #[deprecated(note = "Please use [`par_point_pairs()`](#method.par_point_pairs) instead.")]
    pub fn par_filter_point_pairs<F, P>(&self, f: F, p: P)
    where
        F: Fn(usize, usize) + Send + Sync,
        P: Fn(usize, usize) -> bool + Send + Sync,
    {
        self.par_iter().for_each(|cell| {
            cell.point_pairs()
                .filter(|(i, j)| p(*i, *j))
                .for_each(|(i, j)| {
                    f(i, j);
                })
        });
    }
}
