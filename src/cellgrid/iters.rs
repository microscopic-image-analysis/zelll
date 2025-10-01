//TODO iterate over all neighbored cells (full/half space), pairs of particles
//TODO: perhaps move parallel iteration into separate submodule
use crate::cellgrid::storage::CellSliceMeta;
#[cfg(feature = "rayon")]
use crate::rayon::ParallelIterator;
use crate::{CellGrid, Particle};
use core::iter::FusedIterator;
use core::slice::Iter;
use itertools::Itertools;
use num_traits::{AsPrimitive, ConstOne, Float, NumAssignOps};

// cf. https://predr.ag/blog/definitive-guide-to-sealed-traits-in-rust/#sealing-traits-via-method-signatures
mod private {
    pub struct Token;
}

/// "Marker" trait indicating a type being a valid neighborhood configuration.
pub trait SpaceConfig {
    #[doc(hidden)]
    #[inline]
    fn neighbors_as_slice(neighbors: &[i32], _: private::Token) -> &[i32] {
        neighbors
    }
}

/// _Full-space_ neighborhood.
pub struct Full;
/// _Half-space_ neighborhood.
pub struct Half;

impl SpaceConfig for Full {}
impl SpaceConfig for Half {
    #[inline]
    fn neighbors_as_slice(neighbors: &[i32], _: private::Token) -> &[i32] {
        &neighbors[..neighbors.len() / 2]
    }
}

pub mod neighborhood {
    //! This module enables advanced querying of [`GridCell`](super::GridCell) neighborhoods.
    //!
    //! In `zelll`, cells are represented by _flat_ indices, ie. scalar values.
    //!
    //! Similar to the directions in which the **♚** chess piece can move[^kingdirections],
    //! relative neighbor cell indices are represented by signed integers:
    //!
    // Unfortunately github-flavored markdown tables require a head row, so we include some HTML below
    // |    |        |        |    |    |
    // | -: | -----: | -----: | -: | -: |
    // |    |  **7** |    8   |  9 |    |
    // |    | **-1** | **♚** |  1 |    |
    // |    | **-9** | **-8** | -7 |    |
    // |    |        |        |    |    |
    //!
    //! <style type="text/css">
    //!     tr {text-align: center;}
    //      // hacky way to approximate squares
    //!     td {width: 2.5em; height: 2.5em;}
    //! </style>
    //! <table style="margin-left: auto; margin-right: auto;">
    //!     <tbody>
    //!         <tr>
    //!             <td style="background-color: var(--main-background-color)"><strong> +7 </strong></td>
    //!             <td style="background-color: var(--table-alt-row-background-color)"> +8 </td>
    //!             <td style="background-color: var(--main-background-color)"> +9 </td>
    //!         </tr>
    //!         <tr>
    //!             <td style="background-color: var(--table-alt-row-background-color)"><strong> -1 </strong></td>
    //!             <td style="background-color: var(--main-background-color)"><strong> ♚ </strong></td>
    //!             <td style="background-color: var(--table-alt-row-background-color)"> +1 </td>
    //!         </tr>
    //!         <tr>
    //!             <td style="background-color: var(--main-background-color)"><strong> -9 </strong></td>
    //!             <td style="background-color: var(--table-alt-row-background-color)"><strong> -8 </strong></td>
    //!             <td style="background-color: var(--main-background-color)"> -7 </td>
    //!         </tr>
    //!     </tbody>
    //! </table>
    //!
    //! The _[`Full`]-space_ neighborhood for **♚** simply consists of cells _`{-9, -1, +7, -8, +8, -7, +1, +9}`_
    //! (in this order).\
    //! A _[`Half`]-space_ neighborhood for **♚** consists of cells _`{-9, -1, +7, -8}`_.
    //!
    //! Note that there are multiple valid _half-space_ neighborhoods, this specific sequence
    //! is merely an implementation artifact.
    //! Also, the exact index values depend on the cell grid's size and shape.
    //!
    //! [^kingdirections]: [https://www.chessprogramming.org/Direction#Ray_Directions](https://www.chessprogramming.org/Direction#Ray_Directions)
    #[allow(unused_imports)]
    pub use super::{Full, Half, SpaceConfig};
}

/// `GridCell` represents a possibly empty (by construction) cell of a [`CellGrid`].
#[derive(Debug, Clone, Copy)]
pub struct GridCell<'g, P, const N: usize = 3, F: Float = f64>
where
    F: NumAssignOps + ConstOne + AsPrimitive<i32> + std::fmt::Debug,
    P: Particle<[F; N]>,
{
    //TODO: maybe provide proper accessors to these fields for neighbors.rs to use?
    //TODO: is there a better way than having a reference to the containing CellGrid?
    pub(crate) grid: &'g CellGrid<P, N, F>,
    pub(crate) index: i32,
}

impl<'g, P, const N: usize, F> GridCell<'g, P, N, F>
where
    F: Float + NumAssignOps + ConstOne + AsPrimitive<i32> + Send + Sync + std::fmt::Debug,
    P: Particle<[F; N]> + Send + Sync,
{
    /// Returns the (flat) cell index of this (possibly empty) `GridCell`.
    pub(crate) fn index(&self) -> i32 {
        self.index
    }

    /// Returns an iterator over all particles in this `GridCell`.
    ///
    /// The item type is a pair consisting of the particle index as iterated during `CellGrid`
    /// construction and the particle data itself.
    // TODO: should probably rather impl IntoIterator to match consuming/copy behaviour of neighbors()/point_pairs()?
    pub fn iter(self) -> Iter<'g, (usize, P)> {
        self.grid
            .cell_lists
            .cell_slice(
                self.grid
                    .cells
                    .get(&self.index)
                    // FIXME: cf. `CellStorage::cell_slice()` to see why we have to clone here
                    // CellSliceMeta::default represents an empty slice, so that's exactly what we want
                    // keep in mind that this only works for unique flat cell indices
                    // (ie. anything outside of the `CellGrid` might produce garbage/helical boundaries)
                    .map_or_else(CellSliceMeta::default, |meta| meta.clone()),
            )
            .iter()
    }

    /// Returns an iterator over non-empty neighboring cells.
    ///
    /// <div class="warning">
    ///
    /// Methods such as [`particle_pairs()`](GridCell::particle_pairs()) only access
    /// **half** of the neighboring cells per grid cell
    /// (so-called _half-space_) in order to iterate over _unique_ pairs.
    /// In contrast, [`CellGrid::query_neighbors()`] queries **all** neighboring cells (_full-space_).
    ///
    /// See [`neighborhood`] for details.
    ///
    /// </div>
    ///
    /// # Examples
    /// ```
    /// # use zelll::CellGrid;
    /// use zelll::cellgrid::neighborhood::{Full, Half};
    /// # let points = vec![[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
    /// # let cg = CellGrid::new(points.iter().copied(), 1.0);
    /// # let cell = cg.query([0.5, 1.0, 0.1]).unwrap();
    /// // half-space
    /// cell.neighbors::<Half>().for_each(|_| {});
    /// // full-space
    /// cell.neighbors::<Full>().for_each(|_| {});
    /// ```
    // TODO: currently only aperiodic boundaries
    // TODO: (helical would be simple enough, periodic requires a bit more work)
    pub fn neighbors<S: SpaceConfig>(
        self,
    ) -> impl FusedIterator<Item = GridCell<'g, P, N, F>> + Clone {
        S::neighbors_as_slice(&self.grid.index.neighbor_indices, private::Token)
            .iter()
            .filter_map(move |rel| {
                let index = rel + self.index();
                // TODO: a bit wasteful to throw away the `CellSliceMeta` after the hash map lookup
                // TODO: either make `GridCell` contain this as well or `CellSliceMeta` the cell index too
                // TODO: yet another reason for `cells` to store an enum?
                // TODO: now that iter() supports "empty" cells, filter_map is not even needed anymore
                // TODO: at least internally, have to be careful about this
                self.grid.cells.get(&index).map(|_| GridCell {
                    grid: self.grid,
                    index,
                })
            })
    }

    /// Returns an iterator over all unique pairs of points in this `GridCell`.
    #[inline]
    fn intra_cell_pairs(self) -> impl FusedIterator<Item = ((usize, P), (usize, P))> + Clone {
        // this is equivalent to
        // self.iter().copied().tuple_combinations::<((usize, P), (usize, P))>()
        // but faster for our specific case (pairs from slice of `Copy` values)
        self.iter()
            .copied()
            .enumerate()
            .flat_map(move |(n, i)| self.iter().copied().skip(n + 1).map(move |j| (i, j)))
    }

    /// Returns an iterator over all unique pairs of points in this `GridCell` with points of the neighboring cells.
    #[inline]
    fn inter_cell_pairs(self) -> impl FusedIterator<Item = ((usize, P), (usize, P))> + Clone {
        self.iter().copied().cartesian_product(
            self.neighbors::<Half>()
                .flat_map(|cell| cell.iter().copied()),
        )
    }

    /// Returns an iterator over all _relevant_ pairs of particles within in the neighborhood of this `GridCell`.
    ///
    /// _Relevant_ means the distance between paired particles might be less than the `cutoff` but
    /// this cannot be guaranteed.\
    /// This method consumes `self` but `GridCell` implements [`Copy`].
    //TODO: handle full-space as well
    //TODO: document that we're relying on GridCell impl'ing Copy here (so we can safely consume `self`)
    pub fn particle_pairs(
        self,
    ) -> impl FusedIterator<Item = ((usize, P), (usize, P))> + Clone + Send + Sync {
        self.intra_cell_pairs().chain(self.inter_cell_pairs())
    }
}

impl<P, const N: usize, F> CellGrid<P, N, F>
where
    F: Float + NumAssignOps + ConstOne + AsPrimitive<i32> + Send + Sync + std::fmt::Debug,
    P: Particle<[F; N]>,
{
    /// Returns an iterator over all [`GridCell`]s in this `CellGrid`, excluding empty cells.
    ///
    /// <div class="warning">A particular iteration order is not guaranteed.</div>
    ///
    /// # Examples
    /// ```
    /// # use zelll::CellGrid;
    /// # let points = vec![[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
    /// # let cg = CellGrid::new(points.iter().copied(), 1.0);
    /// assert_eq!(points.len(), cg.iter().flat_map(|cell| cell.iter()).count());
    /// ```
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn iter(&self) -> impl FusedIterator<Item = GridCell<'_, P, N, F>> + Clone {
        // note that ::keys() does not keep a stable iteration order!
        self.cells
            .keys()
            .map(|&index| GridCell { grid: self, index })
    }

    /// Returns a parallel iterator over all [`GridCell`]s in this `CellGrid`, excluding empty cells.
    ///
    /// <div class="warning">A particular iteration order is not guaranteed.</div>
    ///
    /// # Examples
    /// ```
    /// # use zelll::{CellGrid, rayon::ParallelIterator};
    /// # let points = vec![[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
    /// # let cg = CellGrid::new(points.iter().copied(), 1.0);
    /// // The number of non-empty cells in this cell grid does, in fact, not change
    /// // when counted in parallel instead of sequentially.
    /// assert_eq!(cg.iter().count(), cg.par_iter().count());
    /// ```
    #[cfg(feature = "rayon")]
    pub fn par_iter(&self) -> impl ParallelIterator<Item = GridCell<'_, P, N, F>>
    where
        P: Send + Sync,
    {
        self.cells
            .par_keys()
            .map(|&index| GridCell { grid: self, index })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cellgrid::util::generate_pointcloud;

    #[test]
    fn test_cellgrid_iter() {
        // Using 0-origin to avoid floating point errors
        let points = generate_pointcloud([3, 3, 3], 1.0, [0.0, 0.0, 0.0]);
        let cell_grid = CellGrid::new(points.iter().copied(), 1.0);

        assert_eq!(cell_grid.iter().count(), 14, "testing iter()");

        #[cfg(feature = "rayon")]
        assert_eq!(cell_grid.par_iter().count(), 14, "testing par_iter()");
    }

    #[test]
    fn test_gridcell_iter() {
        // Using 0-origin to avoid floating point errors
        let points = generate_pointcloud([3, 3, 3], 1.0, [0.0, 0.0, 0.0]);
        let cell_grid = CellGrid::new(points.iter().copied(), 1.0);

        assert_eq!(
            cell_grid.iter().flat_map(|cell| cell.iter()).count(),
            points.len(),
            "testing iter()"
        );

        #[cfg(feature = "rayon")]
        assert_eq!(
            cell_grid
                .par_iter()
                .flat_map_iter(|cell| cell.iter())
                .count(),
            points.len(),
            "testing par_iter()"
        );
    }

    #[test]
    fn test_neighborcell_particle_pairs() {
        // Using 0-origin to avoid floating point errors
        let points = generate_pointcloud([2, 2, 2], 1.0, [0.0, 0.0, 0.0]);
        let cell_grid = CellGrid::new(points.iter().copied(), 1.0);

        assert_eq!(
            cell_grid
                .iter()
                .map(|cell| cell.intra_cell_pairs().count())
                .sum::<usize>(),
            4,
            "testing intra_cell_pairs()"
        );

        assert_eq!(
            cell_grid
                .iter()
                .map(|cell| cell.inter_cell_pairs().count())
                .sum::<usize>(),
            24,
            "testing inter_cell_pairs()"
        );
    }
}
