//TODO iterate over all neighboured cells (full/half space), pairs of particles
use super::CellGrid;
use ndarray::{Dim, Dimension};

#[derive(Debug)]
pub struct GridCell<'g> {
    grid: &'g CellGrid,
    index: <Dim<[usize; 3]> as Dimension>::Pattern,
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
                index,
                head,
            })
    }
}

//TODO: iterate over neighbored grid cells
//TODO: see https://docs.rs/ndarray/latest/ndarray/struct.ArrayBase.html#method.windows
//TODO: but this doesn't solve the issues with boundary conditions (i.e. wrapping around or partial windows)
//TODO: could introduce empty auxiliary cells on the boundary

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
}

