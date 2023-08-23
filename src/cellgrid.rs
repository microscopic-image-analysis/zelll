//TODO: - allow both &[Point<f64, N>] and impl Iterator<Item = Point<f64, N>> (or IntoIterator)?
//TODO: - possible approach: require impl Iterator and then .take(usize::MAX) to ensure we have only finite point clouds
#[allow(dead_code)]
pub mod flatindex;
#[allow(dead_code)]
pub mod iters;
#[allow(dead_code)]
pub mod neighbors;
#[allow(dead_code)]
pub mod storage;
#[allow(dead_code)]
pub mod util;

pub use flatindex::*;
use hashbrown::HashMap;
pub use iters::*;
use nalgebra::Point;
pub use neighbors::*;
#[cfg(feature = "rayon")]
//TODO: should do a re-export of rayon?
use rayon::prelude::ParallelIterator;
pub use storage::*;
pub use util::*;
//TODO: crate-global type alias for [i32/isize; N] (or [usize; N] if I revert back)
#[derive(Debug, Default, Clone)]
pub struct CellGrid<const N: usize> {
    cells: HashMap<i32, CellSliceMeta>,
    cell_lists: CellStorage<usize>,
    index: FlatIndex<N>,
}

impl<const N: usize> CellGrid<N> {
    pub fn new<'p>(points: impl Iterator<Item = &'p Point<f64, N>> + Clone, cutoff: f64) -> Self {
        CellGrid::default().rebuild(points, Some(cutoff))
    }

    //TODO: Documentation: rebuild does not really update because it does not make a lot of sense:
    //TODO: If FlatIndex did change, we have to re-allocate (or re-initialize) almost everything anyway;
    //TODO: If FlatIndex did not change, we don't need to update.
    //TODO: Therefore we check for that and make CellGrid::new() just a wrapper around CellGrid::rebuild (with an initially empty FlatIndex)
    #[must_use = "rebuild() consumes `self` and returns the rebuilt `CellGrid`"]
    pub fn rebuild<'p>(
        self,
        points: impl Iterator<Item = &'p Point<f64, N>> + Clone,
        cutoff: Option<f64>,
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
                .iter()
                .map(|(idx, count)| (*idx, cell_lists.reserve_cell(*count)))
                .collect();

            index
                .index
                .iter()
                .enumerate()
                //TODO: clean this up, this could be nicer since we know cells.get_mut() won't fail?
                .for_each(|(i, cell)| cell_lists.push(i, cells.get_mut(cell).unwrap()));

            Self {
                cells,
                cell_lists,
                index,
            }
        }
    }

    pub fn rebuild_mut<'p>(
        &mut self,
        points: impl Iterator<Item = &'p Point<f64, N>> + Clone,
        cutoff: Option<f64>,
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
                self.cells
                    .insert_unique_unchecked(*idx, self.cell_lists.reserve_cell(*count));
            });

            //TODO: we'll re-evaluate this once benchmarks with sparse data have been added
            //TODO: However, even without shrinking, self.cells.len() has an upper bound of  self.index.index.len()
            //TODO: (if updated point cloud did not shrink in length)
            self.cells.shrink_to_fit();

            for (i, cell) in self.index.index.iter().enumerate() {
                //TODO: clean this up, this could be nicer since we know cells.get_mut() won't fail?
                //TODO: could use unwrap_unchecked() since this can't/shouldn't fail
                self.cell_lists.push(
                    i,
                    self.cells
                        .get_mut(cell)
                        .expect("cell grid should contain every cell in the grid index"),
                );
            }
        }
    }

    pub fn shape(&self) -> &[i32; N] {
        &self.index.grid_info.shape
    }

    pub fn bounding_box(&self) -> &Aabb<N> {
        &self.index.grid_info.aabb
    }

    /// Iterate over all relevant (i.e. within cutoff threshold + some extra) unique pairs of points in this `CellGrid`.
    /// This has some overhead due to internal use of [`flat_map`].
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn point_pairs(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.iter()
            .flat_map(|cell| cell.point_pairs().collect::<Vec<_>>())
    }

    /// Call a closure on each relevant (i.e. within cutoff threshold + some extra) unique pair of points in this `CellGrid`.
    /// This should be preferred over [`CellGrid::point_pairs()`].
    /// `for_each_point_pair()` is equivalent to a nested loop:
    /// ```ignore
    /// for cell in cell_grid.iter() {
    ///     for (borb, other) in cell.point_pairs() {
    ///         ...
    ///     }
    /// }
    /// ```
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

    pub fn filter_point_pairs<F, P>(&self, mut f: F, p: P)
    where
        F: FnMut(usize, usize),
        P: Fn(usize, usize) -> bool,
    {
        self.iter().for_each(|cell| {
            cell.point_pairs()
                .filter(|(i, j)| p(*i, *j))
                .for_each(|(i, j)| {
                    f(i, j);
                })
        });
    }
}

#[cfg(feature = "rayon")]
impl<const N: usize> CellGrid<N> {
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn par_point_pairs(&self) -> impl ParallelIterator<Item = (usize, usize)> + '_ {
        self.par_iter()
            .flat_map(|cell| cell.point_pairs().collect::<Vec<_>>())
    }

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

