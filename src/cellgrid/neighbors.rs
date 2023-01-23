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
    /// Return the default value for `BalancedTernary`.
    /// For practical reasons, this is not `BalancedTernary::ZERO` but
    /// the balanced ternary representation of `1`.
    #[inline]
    fn default() -> Self {
        // The reason for returning the "1" array is that `from_fn()` is not `const`
        // (actually it's that [expr; N] requires a const expression),
        // i.e. we can't define `BalancedTernary::ONE` as `const` for arbitrary `const N`.
        let mut one = BalancedTernary::ZERO;
        one.0[0] = BalancedTrit::Positive;
        one
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

#[derive(Debug)]
#[must_use = "iterators are lazy and do nothing unless consumed"]
pub struct CellNeighbors<'c, const N: usize> {
    center: &'c GridCell<'c, N>,
    state: BalancedTernary<N>,
}

// This probably the only slightly relevant advantage of using balanced ternary numbers;
// The distinction between half and full space is simply given by the initial state of the iterator
impl<'c, const N: usize> CellNeighbors<'c, N> {
    pub(crate) fn half_space(center: &'c GridCell<'c, N>) -> Self {
        Self {
            center,
            state: BalancedTernary::default(),
        }
    }

    pub(crate) fn full_space(center: &'c GridCell<'c, N>) -> Self {
        Self {
            center,
            state: BalancedTernary::MIN,
        }
    }
}

//TODO: sum type/bool for boundary conditions would evaluated during runtime (pro: could change behavior during runtime)
//TODO: type states for boundary conditions would be during compiletime (pro: no runtime cost, con: slightly less flexible)
impl<'c, const N: usize> Iterator for CellNeighbors<'c, N> {
    type Item = GridCell<'c, N>;

    fn next(&mut self) -> Option<Self::Item> {
        // stop the iterator (by returning None) if the center cell is empty
        // Semantically, empty cells could have non-empty neighbors but we don't need them in that case
        // and empty cells should not be accessible using the public API anyway
        let center = self.center.index()?;

        //TODO: Note the recursive calls to self.next().
        //TODO: I'm okay with that (for now) since in practice the number of neighbor cells is pretty limited.
        //TODO: I'm sure the call stack won't mind too much.
        // See https://doc.rust-lang.org/book/ch18-03-pattern-syntax.html#matching-named-variables
        // and https://doc.rust-lang.org/book/ch18-03-pattern-syntax.html#extra-conditionals-with-match-guards
        match self.state {
            // Iterator stops if the current state is already the last (which is always BalancedTernary::MAX)
            max if max == BalancedTernary::MAX => None,
            // Skip the center cell
            zero if zero == BalancedTernary::ZERO => {
                self.state += BalancedTrit::Positive;
                self.next()
            }
            current_state => {
                // `BalancedTernary` is `Copy` i.e. we can already increment `self.state` and use `current_state`
                self.state += BalancedTrit::Positive;
                //TODO: waiting for stabilization of
                //TODO: https://doc.rust-lang.org/std/primitive.array.html#method.try_map or
                //TODO: https://doc.rust-lang.org/std/iter/trait.Iterator.html#method.try_collect
                //TODO: I'm working around this using mutable state and try_for_each (try_fold would be awkward in this case)
                let mut new_index = [0usize; N];

                // For convenience
                let shape = self.center.grid.index.grid_info.shape;

                new_index
                    .iter_mut()
                    .zip(center.iter())
                    .zip(current_state.0.iter())
                    .zip(shape.iter())
                    .try_for_each(|(((i, c), s), sh)| {
                        if let Some(new_i) = c.checked_add_signed(*s as isize) {
                            if new_i < *sh {
                                *i = new_i;
                                Some(())
                            } else {
                                // upper grid bound violated, skip this cell
                                //TODO: periodic boundaries: wrap around instead
                                None
                            }
                        } else {
                            // Lower grid bound violated, skip this cell
                            //TODO: periodic boundaries: see above
                            None
                        }
                    })
                    .map(|()| GridCell {
                        grid: self.center.grid,
                        head: &self.center.grid.cells[new_index.as_slice()],
                    })
                    .or_else(|| self.next())
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

