//TODO: iterate over individual cell lists, all neighboured cells (full/half space), pairs of particles
use super::CellGrid;
use ndarray::{iter::IndexedIter, Dim, Dimension};

/// Iterates over all cells in a [`CellGrid`], including empty cells.
//TODO: I didn't wrap Filter in my CellGridIterator because I didn't know how to (problematic closure whose type I can't name)
/// Use [`Filter`] to omit e.g. empty cells.
/// The item type is `([usize; 3], &Option<usize>)`.
///
/// # Examples
///
/// ```
/// # use zelll::cellgrid::CellGrid;
/// # let points = [[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
/// # let cell_grid = CellGrid::new(&points, 1.0);
/// cell_grid.cells_iter().filter(|(index, head)| head.is_some());
/// ```
//TODO: I'm only doing this because I don't want to erase the type by returning impl Iterator in cells_iter().
pub struct CellGridIterator<'a> {
    state: IndexedIter<'a, Option<usize>, Dim<[usize; 3]>>,
}

impl CellGrid {
    pub fn cells_iter<'a>(&'a self) -> CellGridIterator<'a> {
        CellGridIterator {
            state: self.cells.indexed_iter(),
        }
    }
}

impl<'a> Iterator for CellGridIterator<'a> {
    type Item = (<Dim<[usize; 3]> as Dimension>::Pattern, &'a Option<usize>);

    fn next(&mut self) -> Option<Self::Item> {
        self.state.next()
    }
}

