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

use super::{GridCell, GridInfo};
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

    //TODO: document that gridinfo contains strides
    pub fn flatten_with_info(&self, info: &GridInfo<N>) -> i32 {
        info.flatten_index(self.0.map(|trit| trit as i32))
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
pub struct RelativeNeighborIndices<const N: usize> {
    state: BalancedTernary<N>,
}

// This is probably the only slightly relevant advantage of using balanced ternary numbers;
// The distinction between half and full space is simply given by the initial state of the iterator
impl<const N: usize> RelativeNeighborIndices<N> {
    pub(crate) fn half_space() -> Self {
        Self {
            state: BalancedTernary::default(),
        }
    }

    pub(crate) fn full_space() -> Self {
        Self {
            state: BalancedTernary::MIN,
        }
    }
}

impl<const N: usize> Iterator for RelativeNeighborIndices<N> {
    type Item = [i32; N];

    fn next(&mut self) -> Option<Self::Item> {
        match self.state {
            // Iterator stops if the current state is already the last (which is always BalancedTernary::MAX)
            max if max == BalancedTernary::MAX => None,
            // Skip the center cell
            zero if zero == BalancedTernary::ZERO => {
                self.state += BalancedTrit::Positive;
                self.next()
            }
            _ => {
                let new_index = self.state.0.map(|trit| trit as i32);
                self.state += BalancedTrit::Positive;
                Some(new_index)
            }
        }
    }
}

impl<const N: usize> FusedIterator for RelativeNeighborIndices<N> {}

#[derive(Debug, Clone, Copy)]
#[must_use = "iterators are lazy and do nothing unless consumed"]
#[deprecated(note = "Please use GridCell::neighbors() directly.")]
//TODO: make sure neighbors() handles half-/full-space and (a)periodic boundaries
pub struct CellNeighbors<'c, const N: usize> {
    center: &'c GridCell<'c, N>,
    //state: BalancedTernary<N>,
    state: usize,
}

// This is probably the only slightly relevant advantage of using balanced ternary numbers;
// The distinction between half and full space is simply given by the initial state of the iterator
//TODO: The distinction between half and full-space now happens in RelativeNeighborIndices<N>
//TODO: Which is done in FlatIndex computation, although this is a bit annoying too
//TODO: could always compute full space neighbors in FlatIndex construction and then start from (number-of-neighbors)/2 in half_space() here
impl<'c, const N: usize> CellNeighbors<'c, N> {
    pub(crate) fn half_space(center: &'c GridCell<'c, N>) -> Self {
        Self {
            center,
            //state: BalancedTernary::default(),
            state: 0,
        }
    }

    pub(crate) fn full_space(center: &'c GridCell<'c, N>) -> Self {
        Self {
            center,
            //state: BalancedTernary::MIN,
            state: 0,
        }
    }

    //TODO: Ideally, I'd use something like Either<L, R>?
    //TODO: we can always compute a new index (i.e. using periodic boundaries)
    //TODO: but would need to indicate if it was wrapped around the boundary
    //TODO: and decide later to omit this if we have aperiodic boundaries?
    /// TODO: Returns None if self.state is out of bounds
    fn absolute_index(&self) -> Option<i32> {
        self.center
            .grid
            .index
            .neighbor_indices
            .get(self.state)
            .map(|rel| *rel + self.center.index())
    }

    //TODO: This is probably a bit tricky with flat indices (+ the strides+1 hack?)
    /*fn absolute_index_periodic(&self) -> [i32; N] {
        let shape = self.center.grid.index.grid_info.shape;
        let mut index = [0; N];
        //TODO: this is still a bit ugly...
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
    }*/
}

//TODO: sum type/bool for boundary conditions would evaluated during runtime (pro: could change behavior during runtime)
//TODO: type states for boundary conditions would be during compiletime (pro: no runtime cost, con: slightly less flexible)
//TODO: could override size_hint() method relatively easily
//TODO: could implement DoubleEndedIterator but that would need extra state and currently there's no reason
impl<'c, const N: usize> Iterator for CellNeighbors<'c, N> {
    type Item = GridCell<'c, N>;
    //TODO: cachegrind attributes most D1mr misses here (likely related to matching the Option, or failed lookup in hashbrown)
    //TODO: Now that we just have a simple Vec containing the relative indices, we probably should just iterate over that
    //TODO: to compute absolute indices, instead of having CellNeighbors manually incrementing self.state
    fn next(&mut self) -> Option<Self::Item> {
        //TODO: Note the recursive calls to self.next().
        //TODO: I'm okay with that (for now) since in practice the number of neighbor cells is pretty limited.
        //TODO: I'm sure the call stack won't mind too much.
        if let Some(new_index) = self.absolute_index() {
            self.state += 1;
            //TODO: maybe make this more clear, we're just checking if cells[&new_index] is occupied
            self.center
                .grid
                .cells
                .get(&new_index)
                .map(|_| GridCell {
                    grid: self.center.grid,
                    index: new_index,
                })
                .or_else(|| self.next())
        } else {
            None
        }
    }
}

impl<const N: usize> FusedIterator for CellNeighbors<'_, N> {}

