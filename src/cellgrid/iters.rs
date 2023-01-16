//TODO iterate over all neighboured cells (full/half space), pairs of particles
//TODO: perhaps move parallel iteration into separate submodule
use super::CellGrid;
#[cfg(feature = "rayon")]
use ndarray::{parallel::prelude::*, Zip};
use ndarray::{Dim, IntoDimension};

#[derive(Debug)]
pub struct GridCell<'g> {
    grid: &'g CellGrid,
    index: Dim<[usize; 3]>, //<Dim<[usize; 3]> as Dimension>::Pattern,
    //TODO: I guess I could just store Option<usize> directly as well
    head: &'g Option<usize>,
}

impl GridCell<'_> {
    pub fn is_empty(&self) -> bool {
        self.head.is_none()
    }

    pub fn iter(&self) -> GridCellIterator {
        GridCellIterator {
            grid: self.grid,
            state: self.head,
        }
    }
}
/// Iterates over all points (or rather their indices) in the [`GridCell`] this `GridCellIterator` was created from.
#[derive(Debug)]
pub struct GridCellIterator<'g> {
    grid: &'g CellGrid,
    //TODO: I guess I could just store Option<usize> directly as well
    state: &'g Option<usize>,
}

impl Iterator for GridCellIterator<'_> {
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

impl CellGrid {
    /// Returns an iterator over all [`GridCell`]s in this `CellGrid`, including empty cells.
    ///
    /// # Examples
    ///
    /// ```
    /// # use zelll::cellgrid::CellGrid;
    /// # let points = [[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
    /// # let cell_grid = CellGrid::new(&points, 1.0);
    /// cell_grid.iter().filter(|cell| !cell.is_empty());
    /// ```
    pub fn iter(&self) -> impl Iterator<Item = GridCell> {
        //TODO: could also filter_map() to remove empty cells here...
        self.cells
            .indexed_iter()
            // It seems a bit weird but I'm just moving a reference to self (if I'm not mistaken).
            .map(move |(index, head)| GridCell {
                grid: self,
                //TODO: could use into_dimensions() in order to not having to store Dimension::Pattern in GridCell?
                index: index.into_dimension(),
                head,
            })
    }

    #[cfg(feature = "rayon")]
    //TODO: using ndarray's Zip like that has overhead but it seems to be the sensible approach
    pub fn par_iter(&self) -> impl ParallelIterator<Item = GridCell> {
        Zip::indexed(&self.cells)
            .into_par_iter()
            .map(move |(index, head)| GridCell {
                grid: self,
                index: index.into_dimension(),
                head,
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const points: &[[f64; 3]] = &[[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.1, 0.2]];

    #[test]
    fn test_cellgrid_iter() {
        let cell_grid: CellGrid = CellGrid::new(points, 1.0);

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
        let cell_grid: CellGrid = CellGrid::new(points, 1.0);

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

