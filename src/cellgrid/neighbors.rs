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
//TODO: turns out I'll use itertools anyway but actually I started liking the manual BalancedTernary approach better
//TODO: Itertools::multi_cartesian_product() uses Vec<> and BalancedTernary is const generic
// (or ndarray::indices(3,3), which has split_at() which might make parallel iteration simpler?)
// (or just ArrayBase.windows() + empty (for half-space)/opposite (for full-space) margin around grid or manual edge case handling).
// Balanced ternary numbers are probably not more performant anyway. It's just a bit convoluted and arguably more fun.
//TODO: would be nice to have a minimal example of these approaches to have a look at on godbolt.org

use super::GridCell;
use core::iter::FusedIterator;
use core::ops::{Add, AddAssign};
use nalgebra::wrap;

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

#[derive(Debug, Clone, Copy)]
#[must_use = "iterators are lazy and do nothing unless consumed"]
pub struct CellNeighbors<'c, const N: usize> {
    center: &'c GridCell<'c, N>,
    state: BalancedTernary<N>,
}

// This is probably the only slightly relevant advantage of using balanced ternary numbers;
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

    //TODO: Ideally, I'd use something like Either<L, R>?
    //TODO: we can always compute a new index (i.e. using periodic boundaries)
    //TODO: but would need to indicate if it was wrapped around the boundary
    //TODO: and decide later to omit this if we have aperiodic boundaries?
    fn absolute_index(&self) -> Option<[i32; N]> {
        let shape = self.center.grid.index.grid_info.shape;
        let mut index = [0; N];
        //TODO: this is really ugly...
        //TODO: would like to use nalgebra::clamp() or PartialOrd::clamp()
        //TODO: but I don't really clamp but check if it is in the clamping interval
        self.center
            .index()
            .iter()
            .zip(self.state.0.iter())
            .map(|(coord, relative)| coord + *relative as i32)
            .zip(shape.iter())
            .zip(index.iter_mut())
            //TODO: waiting for stabilization of
            //TODO: https://doc.rust-lang.org/std/primitive.array.html#method.try_map or
            //TODO: https://doc.rust-lang.org/std/iter/trait.Iterator.html#method.try_collect
            //TODO: I'm working around this using mutable state and try_for_each (try_fold would be awkward in this case)
            .try_for_each(|((coord, dim), new_coord)| {
                if coord < *dim && coord >= 0 {
                    *new_coord = coord;
                    Some(())
                } else {
                    None
                }
            })
            .map(|()| index)
    }

    fn absolute_index_periodic(&self) -> [i32; N] {
        let shape = self.center.grid.index.grid_info.shape;
        let mut index = [0; N];
        //TODO: this is really ugly...
        self.center
            .index()
            .iter()
            .zip(self.state.0.iter())
            .map(|(coord, relative)| coord + *relative as i32)
            .zip(shape.iter())
            .map(|(coord, dim)| wrap(coord, 0, *dim - 1))
            .zip(index.iter_mut())
            .for_each(|(coord, new_coord)| *new_coord = coord);

        index
    }
}

//TODO: sum type/bool for boundary conditions would evaluated during runtime (pro: could change behavior during runtime)
//TODO: type states for boundary conditions would be during compiletime (pro: no runtime cost, con: slightly less flexible)
//TODO: could override size_hint() method relatively easily
//TODO: could implement DoubleEndedIterator but that would need extra state and currently there's no reason
impl<'c, const N: usize> Iterator for CellNeighbors<'c, N> {
    type Item = GridCell<'c, N>;
    //TODO: cachegrind attributes many D1mr misses here (related to copy_nonoverlapping, AddAssign, intrinsics.rs ptr alignment?, mixed_integer_ops, checked add?, equality.rs)
    fn next(&mut self) -> Option<Self::Item> {
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
            _ => {
                let new_index = self.absolute_index();

                self.state += BalancedTrit::Positive;

                new_index
                    .and_then(|idx| self.center.grid.cells.get(&idx).copied())
                    .map(|head| GridCell {
                        grid: self.center.grid,
                        head,
                    })
                    .or_else(|| self.next())
            }
        }
    }
}

impl<const N: usize> FusedIterator for CellNeighbors<'_, N> {}

#[cfg(test)]
mod tests {
    //use super::*;

    #[test]
    fn test_neighbors() {
        todo!();
    }
}

