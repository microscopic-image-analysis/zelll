//TODO iterate over all neighboured cells (full/half space), pairs of particles
//TODO: perhaps move parallel iteration into separate submodule
use super::{CellGrid, CellNeighbors};
#[cfg(feature = "rayon")]
use ndarray::parallel::prelude::*;

#[derive(Debug)]
pub struct GridCell<'g, const N: usize> {
    //TODO: maybe provide proper accessors to these fields for neighbors.rs to use?
    pub(crate) grid: &'g CellGrid<N>,
    //TODO: I guess I could just store Option<usize> directly as well
    pub(crate) head: &'g Option<usize>,
}

impl<const N: usize> GridCell<'_, N> {
    pub fn is_empty(&self) -> bool {
        self.head.is_none()
    }

    /// Return the (multi-)index of this (non-empty) `GridCell`.
    /// Returns `None` if the cell is empty.
    //TODO: However, the public API should not provide a way to address empty cells
    pub(crate) fn index(&self) -> Option<[usize; N]> {
        let idx = (*self.head)?;
        Some(self.grid.index.index[idx])
    }

    pub fn iter(&self) -> GridCellIterator<N> {
        GridCellIterator {
            grid: self.grid,
            state: self.head,
        }
    }

    /// Check whether this `GridCell` is on the boundary of the [`CellGrid`].
    /// Returns None if the cell is empty.
    //TODO: Again, empty cells shouldn't be accessible via the public API
    //TODO: I don't think I need it but let's keep it anyway
    pub fn on_boundary(&self) -> Option<bool> {
        let idx = self.index()?;
        let shape = self.grid.index.grid_info.shape;

        for (i, dim) in shape.iter().enumerate() {
            if idx[i] == 0 || idx[i] + 1 == *dim {
                return Some(true);
            }
        }
        Some(false)
    }

    /// Return [`CellNeighbors`], an iterator over all (currently half-space) non-empty neighboring cells.
    //TODO: currently only half-space and aperiodic boundaries
    //TODO: handle half-/full-space  and (a-)periodic boundary conditions
    pub fn neighbors(&self) -> CellNeighbors<N> {
        CellNeighbors::half_space(self)
    }
}
/// Iterates over all points (or rather their indices) in the [`GridCell`] this `GridCellIterator` was created from.
#[derive(Debug)]
#[must_use = "iterators are lazy and do nothing unless consumed"]
pub struct GridCellIterator<'g, const N: usize> {
    grid: &'g CellGrid<N>,
    //TODO: I guess I could just store Option<usize> directly as well
    state: &'g Option<usize>,
}

impl<const N: usize> Iterator for GridCellIterator<'_, N> {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(index) = self.state {
            self.state = &self.grid.cell_lists[*index];
            Some(*index)
        } else {
            None
        }
    }
}

impl<const N: usize> CellGrid<N> {
    /// Returns an iterator over all [`GridCell`]s in this `CellGrid`, excluding empty cells.
    ///
    /// # Examples
    ///
    /// ```
    /// # use zelll::cellgrid::CellGrid;
    /// # let points = [[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
    /// # let cell_grid = CellGrid::new(&points, 1.0);
    /// cell_grid.iter().filter(|cell| !cell.is_empty());
    /// ```
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn iter(&self) -> impl Iterator<Item = GridCell<N>> {
        self.cells
            .iter()
            // It seems a bit weird but I'm just moving a reference to self (if I'm not mistaken).
            .filter_map(move |head| {
                if head.is_some() {
                    Some(GridCell { grid: self, head })
                } else {
                    None
                }
            })
    }

    #[cfg(feature = "rayon")]
    #[must_use = "iterators are lazy and do nothing unless consumed"]
    pub fn par_iter(&self) -> impl ParallelIterator<Item = GridCell<N>> {
        self.cells
            .par_iter()
            // It seems a bit weird but I'm just moving a reference to self (if I'm not mistaken).
            .filter_map(move |head| {
                if head.is_some() {
                    Some(GridCell { grid: self, head })
                } else {
                    None
                }
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const points: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.1, 0.2]];

    #[test]
    fn test_cellgrid_iter() {
        let cell_grid: CellGrid<3> = CellGrid::new(points, 1.0);

        for cell in cell_grid.iter() {
            println!("{:?}", cell);
        }

        // Doing it twice to check if there were issues with moving &self
        for cell in cell_grid.iter() {
            println!("{:?}", cell);
        }
    }

    #[test]
    fn test_gridcell_iter() {
        let cell_grid: CellGrid<3> = CellGrid::new(points, 1.0);

        for cell in cell_grid.iter() {
            for point in cell.iter() {
                println!("{}", point);
            }
        }
    }

    #[cfg(feature = "rayon")]
    #[test]
    fn test_cellgrid_par_iter() {
        let cell_grid: CellGrid = CellGrid::new(points, 1.0);

        for cell in cell_grid.iter() {
            println!("{:?}", cell);
        }
    }
}

