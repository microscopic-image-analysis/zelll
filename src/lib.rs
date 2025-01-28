//! `zelll`[^etymology] provides a Rust implementation of the __cell lists__ algorithm.
//!
//! Particle simulations usually require to compute interactions between those particles.
//! Considering all _pairwise_ interactions of _`n`_ particles would be of time complexity _`O(n²)`_.\
//! Cell lists facilitate _linear-time_ enumeration of particle pairs closer than a certain
//! cutoff distance by dividing the enclosing bounding box into cuboid grid cells.
//!
//! # Caveats
//!
//! `zelll` is motivated by _coarse-grained_ (bio-)molecular simulations but is not restricted to that.\
//! This is reflected by a few things:
//!
//! - internally, the simulation box is represented by a (sparse) hash map only storing non-empty grid cells,
//!   which gives an upper bound for memory usage given by _`n`_
//! - bounding boxes are assumed to change and are computed from particle data
//!   (future APIs may be added to set a fixed bounding box)
//! - instead of cell _lists_, slices into a contiguous storage buffer are used
//! - periodic boundary conditions are currently not supported
//! - parts of this implementation are more cache-aware than others which becomes noticeable with\
//!   larger data sets (at `10⁶` -- `10⁷` particles, mostly depending on L2 cache size)\
//!   but is less pronounced with sequential (or otherwise structured) data, such as biomolecules
//!
//! # Usage
//!
//! The general pattern in which this crate is intended to be used is roughly:
//!
//! 1. construct `CellGrid` from particle positions
//! 2. enumerate pairs in order to compute particle interactions
//! 3. simulate particle motion
//! 4. rebuild `CellGrid` from updated particle positions
//!
//! This crate only provides iteration over particle pairs.
//! It is left to the user to filter (eg. by distance) and compute interaction potentials.
//! The `rayon` feature enables parallel iteration. Performance gains depend on data size and
//! computational cost per pair though. Benchmarks are encouraged.
//!
//! While the main struct [`CellGrid`] is generic over dimension `N`,
//! it is intended to be used with `N = 2` or `N = 3`.
//! Particle data represented as fixed-size arrays is supported without additional work.\
//! Additionally, implementing [`Particle`] allows usage of custom types as particle data.
//! This can be used to encode different kinds of particles or enable interior mutability if required.
//!
//! # Examples
//! ```
//! use zelll::CellGrid;
//!
//! let points = vec![[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
//! let mut cg = CellGrid::new(points.iter().copied(), 1.0);
//!
//! for ((i, p), (j, q)) in cg.particle_pairs() {
//!     /* do some work */
//! }
//!
//! cg.rebuild_mut(points.iter().copied(), Some(0.5));
//! ```
//!
//! [^etymology]: abbrv. from German _Zelllisten_ /ˈʦɛlɪstən/, for cell lists.
#[allow(dead_code)]
pub mod cellgrid;

#[cfg(feature = "rayon")]
pub mod rayon {
    //! Re-export of the [`ParallelIterator`] trait.
    pub use rayon::prelude::ParallelIterator;
}

// inlined re-exports
#[doc(inline)]
pub use crate::cellgrid::CellGrid;

/// Particle data trait.
///
/// This trait is required for types used with [`CellGrid`] which needs how to get coordinate data.\
/// Only [`Copy`] types can be used.
/// In general, the smaller the type, the better (for the CPU cache).
///
/// A blanket implementation for `Into<T> + Copy` types is provided.\
/// [`CellGrid`] is slightly more specific and requires impl'ing `Particle<[{float}; N]>`.
/// Therefore, fixed-size float arrays, [`nalgebra::SVector`], or types that can be `Deref`-coerced
/// into the former or [`mint`](https://docs.rs/mint/latest/mint/) types, can be directly used.
///
/// Having custom types implement this trait allows for patterns like interior mutability,
/// referencing separate storage (eg. with ECS, or concurrent storage types),
/// or particle data having differend kinds.\
///
/// # Examples
/// ```
/// # use zelll::Particle;
/// # #[derive(Clone, Copy)]
/// # enum AtomKind {
/// #    Hydrogen, // no associated coordinate data since it's the same for all atom kinds
/// #    Oxygen,
/// #    // ...
/// # }
/// #
/// // Typically, we would associate data with `AtomKind` variants for concise code
/// // but here, every variant would carry the same type of data.
/// #[derive(Clone, Copy)]
/// struct Atom {
///     kind: AtomKind,
///     coords: [f64; 3],
/// }
/// impl Particle for Atom {
///     #[inline]
///     fn coords(&self) -> [f64; 3] {
///         self.coords // no pattern matching required here
///     }
/// }
/// ```
pub trait Particle<T = [f64; 3]>: Copy {
    /// Return a copy of this particle's coordinates
    fn coords(&self) -> T;
}

// TODO: Might consider restricting this impl.
// TODO: While we can be this generic, this might help articulating our intentions better:
// TODO: impl<P, T, const N: usize> Particle<[T; N]> for P where P: Into<[T; N]> + Copy {
impl<P, T> Particle<T> for P
where
    P: Into<T> + Copy,
{
    #[inline]
    fn coords(&self) -> T /* [T; N] */ {
        <P as Into<T>>::into(*self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::SVector;

    // TODO: API experiments for GridStorage trait in CellGrid
    trait Grid<T = ()>: Default {}

    #[derive(Default)]
    struct SparseGrid;
    #[derive(Default)]
    struct DenseGrid;

    impl<T> Grid<T> for SparseGrid {}
    impl<T> Grid<T> for DenseGrid {}

    struct PStorage<P, G: Grid<()> = SparseGrid> {
        buffer: Vec<P>,
        grid: G,
    }

    impl<P> PStorage<P, SparseGrid> {
        #[inline]
        fn new_sparse<I>(points: I) -> Self
        where
            I: IntoIterator<Item = P>,
            P: Copy,
        {
            PStorage::new(points)
        }
    }

    impl<P> PStorage<P, DenseGrid> {
        #[inline]
        fn new_dense<I>(points: I) -> Self
        where
            I: IntoIterator<Item = P>,
            P: Copy,
        {
            PStorage::new(points)
        }
    }

    impl<P, G: Grid<()>> PStorage<P, G> {
        fn new<I>(points: I) -> Self
        where
            I: IntoIterator<Item = P>,
            P: Copy,
        {
            Self {
                buffer: points.into_iter().collect(),
                grid: <G as Default>::default(),
            }
        }

        fn convert<T>(&self) -> Vec<T>
        where
            P: Particle<T>,
        {
            self.buffer.iter().map(|p| p.coords()).collect()
        }
    }

    #[test]
    fn test_impl_particle() {
        let points = vec![[0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3]];

        let ps = PStorage::<_, SparseGrid>::new(points.iter().copied());
        let ps: PStorage<_> = PStorage::new(points.clone().into_iter());

        let points: Vec<_> = points.into_iter().map(|p| SVector::from(p)).collect();
        let ps = PStorage::new_sparse(points.iter().copied());

        let _: Vec<[_; 3]> = ps.convert();
    }

    #[test]
    fn test_impl_particle_ref() {
        #[derive(Clone, Copy)]
        struct ParticleRef<'p>(&'p [f64; 3]);

        impl Particle<[f64; 3]> for ParticleRef<'_> {
            #[inline]
            fn coords(&self) -> [f64; 3] {
                (*self.0).coords() // equivalent to *self.0
            }
        }

        let points = vec![[0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3]];

        let ps: PStorage<ParticleRef> =
            PStorage::new(points.iter().map(|p| ParticleRef(p)).clone());
        let _: Vec<[_; 3]> = ps.convert();
    }
}
