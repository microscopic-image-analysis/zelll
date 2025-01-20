/// TODO

///cell grid stuff
#[allow(dead_code)]
pub mod cellgrid;

// TODO: could reduce boilerplate in current trait bounds
// TODO: implementing this trait on a custom type allows to store metadata for each point
// TODO: e.g. having an enum WaterAtom::{Hydrogen([T, N]), Oxygen([T, N])} impl Particle
// TODO: allows to filter easily by those variants
// TODO: another example would be enum Nucleotide::{A(...), C(...), G(...), U(...)}
// TODO: or enum DiploidChromosome::{A(...), B(...)}
// TODO: maybe make it a subtrait of Copy? (check.)
// TODO: implementors could choose to make their type only be a reference/slice (check. cf. tests below)
// TODO: which means we would only store a reference in our cached CellStorage
// TODO: (almost as before when we didn't store anything besides the index)
// TODO: which would make sense for atomically/Mutex protected access
// TODO: but usually implementors would keep their type small so it can be cached directly in CellStorage

pub trait Particle<T = [f64; 3]>: Copy {
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
