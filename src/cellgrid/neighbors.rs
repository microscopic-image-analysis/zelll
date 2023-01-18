// In order to enumerate all neighboring cells of an individual GridCell in N dimensions
// we iterate through a 3^N box centered around that GridCell.
// Starting at a corner of that box we can use ternary numbers where each trit corresponds to a coordinate.
// (This could be generalized to "hierarchical" neighbors, i.e.  5-ary numbers if we want to include neighbors of neighbors)
// (could also use ndarray::indices((3,3)) for this and then map to balanced trits... would be more sensible)
// However, we'd prefer neighbor coordinates relative to the center.
// To do that, we map trits (0, 1, 2) to balanced trits (-1, 0, 1).
// Enumerating the neighboring cells' indices  is now simply done by incrementing a ternary number iteratively by 1,
// starting from [0; N] for half-space neighbors and [-1; N] for full-space neighbors until the upper boundary for our trit length is reached.
//
// This is equivalent to the cartesian product (-1, 0, 1)^N for which we could use itertools.
// Long story short: Using balanced ternary numbers saves us a (single) dependency.
// Might fall back to itertools if we end up needing it anyway
// (or ndarray::indices(3,3), which has split_at() which might parallel iteration simpler?)
// (or just ArrayBase.windows() + empty (for half-space)/opposite (for full-space) margin around grid or manual edge case handling).
// Balanced ternary numbers are probably not more performant anyway. It's just a bit convoluted and arguably more fun.
//TODO: would be nice to have a minimal example of these approaches to have a look at on godbolt.org

use super::GridCell;
use core::ops::{Add, AddAssign};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default)]
#[repr(i8)]
enum BalancedTrit {
    Negative = -1,
    #[default]
    Zero = 0,
    Positive = 1,
}

impl Add for BalancedTrit {
    // Carrying trit, actual trit
    type Output = (Self, Self);

    fn add(self, rhs: Self) -> Self::Output {
        use BalancedTrit::{Negative, Positive, Zero};
        match (self, rhs) {
            (Positive, Negative) | (Negative, Positive) => (Zero, Zero),
            (Positive, Positive) => (Positive, Negative),
            (Negative, Negative) => (Negative, Positive),
            (Zero, val) | (val, Zero) => (Zero, val),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct BalancedTernary<const N: usize>([BalancedTrit; N]);

impl<const N: usize> BalancedTernary<N> {
    const ZERO: BalancedTernary<N> = BalancedTernary([BalancedTrit::Zero; N]);
    const MAX: BalancedTernary<N> = BalancedTernary([BalancedTrit::Positive; N]);
    const MIN: BalancedTernary<N> = BalancedTernary([BalancedTrit::Negative; N]);

    pub const fn trit_len(&self) -> usize {
        self.0.len()
    }
}

impl<const N: usize> Default for BalancedTernary<N> {
    #[inline]
    fn default() -> Self {
        Self::ZERO
    }
}

//TODO: just implementing this for isize here
//TODO: could also do it for other integral types using num or num_traits
//TODO: This is infallible but will wrap around if BalancedTernary<N> is too large
impl<const N: usize> From<BalancedTernary<N>> for isize {
    fn from(value: BalancedTernary<N>) -> Self {
        let mut place_value = 1isize;

        value.0.as_slice().iter().fold(0isize, |acc, trit| {
            let new = acc.wrapping_add(*trit as isize).wrapping_mul(place_value);
            place_value *= 3;
            new
        })
    }
}
//TODO: note that this does not over-/underflow but caps at MAX and MIN respectively
impl<const N: usize> AddAssign<BalancedTrit> for BalancedTernary<N> {
    fn add_assign(&mut self, rhs: BalancedTrit) {
        let mut carry = rhs;
        self.0.iter_mut().for_each(|trit| {
            (carry, *trit) = carry + *trit;
        });
    }
}

//TODO: note that this does not over-/underflow but caps (saturates) at MAX and MIN respectively
//TODO: see also https://doc.rust-lang.org/std/primitive.usize.html#method.overflowing_add
//TODO: and https://doc.rust-lang.org/std/primitive.usize.html#method.carrying_add
impl<const N: usize> Add<BalancedTrit> for BalancedTernary<N> {
    type Output = Self;
    fn add(self, rhs: BalancedTrit) -> Self::Output {
        let mut carry = rhs;

        BalancedTernary(self.0.map(|mut trit| {
            (carry, trit) = carry + trit;
            trit
        }))
    }
}

//TODO: To make this work properly in arbitrary dimensions, I need to address generic CellGrids, GridCells, MultiIndex etc. first
//TODO: and/or reintroduce ndarray::Dim instead of using const arrays
//TODO: The problem being that I keep references to 3D grid in each GridCell
//TODO: next step: introduce const N: usize constraints
#[derive(Debug)]
#[must_use = "iterators are lazy and do nothing unless consumed"]
//TODO: we already have a const parameter here.
//TODO: At some point GridCell will have this too (or sth. like that)
pub struct CellNeighbors<'c> {
    center: &'c GridCell<'c>,
    state: BalancedTernary<3>,
}

impl<'c> CellNeighbors<'c> {
    pub(crate) fn half_space(center: &'c GridCell<'c>) -> Self {
        Self {
            center,
            state: BalancedTernary::ZERO,
        }
    }

    pub(crate) fn full_space(center: &'c GridCell<'c>) -> Self {
        Self {
            center,
            state: BalancedTernary::MIN,
        }
    }
}

impl<'c> Iterator for CellNeighbors<'c> {
    type Item = GridCell<'c>;
    //TODO: this definitively getting rewritten...
    fn next(&mut self) -> Option<Self::Item> {
        // We shouldn't be able to access empty cells via the public API
        // If that happens this iterator is stopping because it's not really useful anyway
        let center_index = self.center.index()?;
        if self.state == BalancedTernary::MAX {
            None
        } else if self.state == BalancedTernary::ZERO {
            // skip the center cell
            self.state += BalancedTrit::Positive;
            self.next()
        } else {
            let mut index: [usize; 3] = core::array::from_fn(|i| i);
            for i in index {
                //check for underflow (i.e. negative index coordinates)
                if let Some(coord) = center_index[i].checked_add_signed(self.state.0[i] as isize) {
                    // check for overflow (out of grid dimensions)
                    if coord < self.center.grid.index.grid_info.shape[i] {
                        index[i] = coord;
                    } else {
                        //TODO: also aperiodic boundary
                        //TODO: refactor this, it's getting ugly and repetitive
                        self.state += BalancedTrit::Positive;
                        return self.next();
                    }
                } else {
                    //TODO: currently only aperiodic boundary: increment self.state and call self.next() again
                    //TODO: for periodic boundary: wrap out-of-bounds coordinates around and return new GridCell
                    self.state += BalancedTrit::Positive;
                    return self.next();
                }
            }

            // check whether the neighboring cell at the newly constructed index is empty
            if self.center.grid.cells[index].is_some() {
                Some(GridCell {
                    grid: self.center.grid,
                    head: &self.center.grid.cells[index],
                })
            } else {
                //TODO: again, refactoring needed
                // we skip empty neighboring cells
                self.state += BalancedTrit::Positive;
                self.next()
            }
        }
    }
}

#[cfg(test)]
mod tests {
    //use super::*;

    #[test]
    fn test_skeleton() {
        todo!();
    }
}

