// FIXME: derived serde traits do not yet show up as feature-gated
// FIXME: cf. https://github.com/rust-lang/rust/issues/103300
#![cfg_attr(docsrs, feature(doc_cfg))]
//! `zelll`[^etymology] provides a Rust implementation of the __cell lists__ algorithm.
//!
//! Particle simulations usually require to compute interactions between those particles.
//! Considering all _pairwise_ interactions of _`n`_ particles would be of time complexity _`O(n²)`_.\
//! Cell lists facilitate _linear-time_ enumeration of particle pairs closer than a certain
//! cutoff distance by dividing the enclosing bounding box into (cuboid) grid cells.
//!
//! # Caveats
//!
//! `zelll` is motivated by _coarse-grained_ (bio-)molecular simulations but is not restricted to that.\
//! This is reflected by a few points:
//!
//! - internally, the simulation box is represented by a (sparse) hash map only storing non-empty grid cells,
//!   which gives an upper bound for memory usage of _`n`_
//! - bounding boxes are assumed to change and are computed from particle data\
//!   (future APIs may be added to set a fixed bounding box)
//! - instead of cell _lists_, slices into a contiguous storage buffer are used
//! - periodic boundary conditions are currently not supported
//! - parts of this implementation are more cache-aware than others, which becomes noticeable with
//!   larger data sets\
//!   (at `10⁶` -- `10⁷` particles, mostly depending on last-level cache size)
//!   but is less pronounced with structured data[^structureddata]
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
//! It is left to the user to filter (e.g. by distance) and compute interaction potentials.
//! The `rayon` feature enables parallel iteration. Performance gains depend on data size and
//! computational cost per pair though. Benchmarks are encouraged.
//!
//! While the main struct [`CellGrid`] is generic over dimension `N`,
//! it is intended to be used with `N = 2` or `N = 3`.
//! Data represented as fixed-size arrays is supported by wrapping in [`Particle`].\
//! Additionally, implementing [`ParticleLike`] allows usage of custom types as particle data.
//! This can be used to encode different kinds of particles or enable interior mutability if required.
//!
//! # Examples
//! ```
//! use zelll::CellGrid;
//!
//! let data = vec![[0.0, 0.0, 0.0], [1.0,2.0,0.0], [0.0, 0.1, 0.2]];
//! let mut cg = CellGrid::new(data.iter().copied(), 1.0);
//!
//! for ((i, p), (j, q)) in cg.particle_pairs() {
//!     /* do some work */
//! }
//!
//! cg.rebuild_mut(data.iter().copied(), Some(0.5));
//! ```
//!
//! [^etymology]: abbrv. from German _Zelllisten_ /ˈʦɛlɪstən/, for cell lists.
//! [^structureddata]: Usually, (bio-)molecular data files are not completely unordered
//! even though they could be.
//! In practice, it may be a reasonable assumption that sequentially proximate
//! particles often have spatially clustered coordinates as well.
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
/// This trait is required for types used with [`CellGrid`]
/// which needs to know how to get coordinate data.\
/// Only [`Copy`] types can be used.
/// In general, the smaller the type, the better (for the CPU cache).
///
/// Note that [`CellGrid`] is more specific than this trait and requires implementing `ParticleLike<[{float}; N]>`.
///
/// We do not provide a blanket implementation for types implementing `Into<[T; N]> + Copy` but
/// a wrapper type instead.
/// Therefore, [`nalgebra::SVector`](https://docs.rs/nalgebra/latest/nalgebra/base/type.SVector.html), or types that can be `Deref`-coerced
/// into the former or [`mint`](https://docs.rs/mint/latest/mint/) types, can be used
/// by wrapping them in [`Particle`].
///
/// For convenience, this trait is implemented for `[{float}; N]` and key-value tuples.
///
/// Having custom types implement this trait allows for patterns like interior mutability,
/// referencing separate storage (e.g. with ECS, or concurrent storage types),
/// or particle data being of different kinds.
///
/// # Examples
/// ```
/// # use zelll::ParticleLike;
/// # #[derive(Clone, Copy)]
/// # enum Element {
/// #    Hydrogen, // no associated coordinate data since it's the same for all variants
/// #    Oxygen,
/// #    // ...
/// # }
/// #
/// // Typically, we would associate data with `Element` variants for concise code
/// // but here, every variant would carry the same type of data.
/// #[derive(Clone, Copy)]
/// struct Atom {
///     kind: Element,
///     coords: [f64; 3],
/// }
/// impl ParticleLike for Atom {
///     #[inline]
///     fn coords(&self) -> [f64; 3] {
///         self.coords // no pattern matching required here
///     }
/// }
/// ```
pub trait ParticleLike<T = [f64; 3]>: Copy {
    /// Returns a copy of this particle's coordinates.
    fn coords(&self) -> T;
}

/// Wrapper type that implements [`ParticleLike`] for types that are `Into<[T; N]> + Copy`.
///
/// Notable types that can be used with `Particle` include `({float}, ...)`,
/// [`nalgebra::SVector`](https://docs.rs/nalgebra/latest/nalgebra/base/type.SVector.html), or types
/// that can be `Deref`-coerced into the former or [`mint`](https://docs.rs/mint/latest/mint/) types,
/// respectively.
///
/// `Particle` can be `Deref`-coerced into the wrapped type.
///
/// # Examples
/// ```
/// use zelll::{Particle, ParticleLike};
///
/// let x: Particle<_> = (0.0f64, 0.0f64, 0.0f64).into();
/// x.coords();
/// ```
#[derive(Copy, Clone, Debug, Default)]
pub struct Particle<P> {
    inner: P,
}

impl<P, T: Copy, const N: usize> ParticleLike<[T; N]> for Particle<P>
where
    P: Into<[T; N]> + Copy,
{
    #[inline]
    fn coords(&self) -> [T; N] {
        <P as Into<[T; N]>>::into(self.inner)
    }
}

/// # Examples
/// Wrapping an array is redundant but illustrates the `Deref` behavior:
/// ```
/// # use zelll::{Particle, ParticleLike};
/// let x = Particle::from([0.0f64; 3]);
/// x.coords(); // calls impl ParticleLike<T> for Particle<P>
/// (*x).coords(); // calls impl ParticleLike<T> for P where P: ParticleLike<T>
/// ```
impl<P> std::ops::Deref for Particle<P> {
    type Target = P;

    fn deref(&self) -> &Self::Target {
        &self.inner
    }
}

impl<P> From<P> for Particle<P> {
    fn from(value: P) -> Self {
        Self { inner: value }
    }
}

/// Sometimes, it is useful to store indices alongside particle data.
/// This is easily facilitated by enumerating the iterator used to construct a `CellGrid`.
///
/// Other types such as [`std::collections::HashMap`] can also be used
/// although it is less straight-forward
///
/// # Examples
/// Here, `ip` is a tuple `(usize, P)`.
/// ```
/// # use zelll::ParticleLike;
/// # let x = [0.0f64; 3];
/// # let points: Vec<_> = std::iter::repeat_n(x, 10).collect();
/// points.iter()
///     .copied()
///     .enumerate()
///     .map(|ip| ip.coords());
/// ```
/// ```
/// # use zelll::ParticleLike;
/// # use std::collections::HashMap;
/// # let data = [("a", [0.0f64; 3]); 10];
/// let points = HashMap::from(data);
/// // this iterator is consuming
/// points.into_iter()
///     .map(|kv| kv.coords());
/// // after constructing a CellGrid
/// // and iterating over pairs, the particles have to be inserted
/// // into a hash map again to do any meaningful work
/// ```
impl<L, P, T, const N: usize> ParticleLike<[T; N]> for (L, P)
where
    L: Copy,
    P: ParticleLike<[T; N]>,
{
    #[inline]
    fn coords(&self) -> [T; N] {
        self.1.coords()
    }
}

impl<T, const N: usize> ParticleLike<[T; N]> for [T; N]
where
    T: Copy,
{
    #[inline]
    fn coords(&self) -> [T; N] {
        *self
    }
}

#[allow(dead_code)]
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
            P: ParticleLike<T>,
        {
            self.buffer.iter().map(|p| p.coords()).collect()
        }
    }

    #[test]
    fn test_wrapped_particle() {
        // let x: Particle<_> = [0.0f64; 3].into();
        let x = Particle::from([0.0f64; 3]);
        x.coords(); // calls impl ParticleLike<T> for Particle<P>
        (*x).coords(); // calls impl ParticleLike<T> for P where P: ParticleLike<T>

        let points = std::iter::repeat_n(x, 10);

        let _ = points.clone().enumerate().map(|ix| ix.coords());
        let _ = points.enumerate().map(|(_i, x)| x.coords());

        let points = std::iter::repeat_n([0.0f64; 3], 10);

        let _ = points.clone().enumerate().map(|ix| ix.coords());
        let _ = points.enumerate().map(|(_i, x)| x.coords());

        let points = std::iter::repeat_n(SVector::from([0.0f64; 3]), 10)
            .map(Particle::from)
            .enumerate();

        let _ = points.clone().map(|ix| -> [f64; 3] { ix.coords() });

        let _ = points.map(|(_i, x)| -> [f64; 3] { x.coords() });
    }

    #[test]
    fn test_impl_particle() {
        let points = vec![[0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3]];

        let _ps = PStorage::<_, SparseGrid>::new(points.iter().copied());
        let _ps: PStorage<_> = PStorage::new(points.clone().into_iter());

        let points: Vec<_> = points
            .into_iter()
            .map(SVector::from)
            .map(Particle::from)
            .collect();
        let ps = PStorage::new_sparse(points.iter().copied());

        let _: Vec<[_; 3]> = ps.convert();
    }

    #[test]
    fn test_impl_particle_ref() {
        #[derive(Clone, Copy)]
        struct ParticleRef<'p>(&'p [f64; 3]);

        impl ParticleLike<[f64; 3]> for ParticleRef<'_> {
            #[inline]
            fn coords(&self) -> [f64; 3] {
                (*self.0).coords() // equivalent to *self.0
            }
        }

        let points = vec![[0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3]];

        let ps: PStorage<ParticleRef> = PStorage::new(points.iter().map(ParticleRef).clone());
        let _: Vec<[_; 3]> = ps.convert();
    }
}
