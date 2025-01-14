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
        (*self).into()
    }
}

struct PStorage<P> {
    buffer: Vec<P>,
}

impl<P> PStorage<P> {
    fn new<Q, I>(points: I) -> Self
    where
        I: IntoIterator<Item = Q>,
        P: Copy,
        Q: std::borrow::Borrow<P>,
    {
        Self {
            buffer: points.into_iter().map(|p| *p.borrow()).collect(),
        }
    }

    fn convert<T>(&self) -> Vec<T>
    where
        P: Particle<T>,
    {
        let out = self.buffer.iter().map(|p| p.coords()).collect();

        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::SVector;

    #[test]
    fn test_impl_particle() {
        let p0 = [0.0f64; 3];
        let p1 = &p0;
        let p2 = SVector::from(p0);

        fn test<P>(point: P) -> [f64; 3]
        where
            P: Particle<[f64; 3]>,
            // Q: std::borrow::Borrow<P>,
        {
            // (*point.borrow()).coords()
            point.coords()
        }

        let points = vec![[0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3]];

        points.iter().for_each(|&p| assert_eq!(test(p), p));

        assert_eq!(test(p0), p0);
        assert_eq!(test(*p1), *p1);
        // assert_eq!(test(p1), *p1);
        assert_eq!(test(p2), p0);
        // assert_eq!(test(&p2), p0);
        assert_eq!(test(*&p2), p0);

        let _: [_; 3] = p0.coords();
        let _: [_; 3] = p2.coords();
    }

    #[test]
    fn test_impl_particle2() {
        let points = vec![[0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3]];

        let ps: PStorage<[f32; 3]> = PStorage::new(points.iter().clone());
        let ps: PStorage<[f32; 3]> = PStorage::new(points.into_iter().clone());

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
