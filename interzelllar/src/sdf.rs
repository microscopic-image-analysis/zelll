use crate::Angstrom;
use crate::atom::Atom;
use nalgebra::{ComplexField, Const, Point, SVector};
use num_dual::*;
use zelll::{CellGrid, Particle};

type DsVec<T> = DualVec<T, T, Const<3>>;

pub struct SmoothDistanceField {
    inner: CellGrid<Atom, 3, Angstrom>,
}

impl SmoothDistanceField {
    /// Returns distance and gradient at `pos`.
    pub fn evaluate(&self, pos: [Angstrom; 3]) -> Option<(Angstrom, [Angstrom; 3])> {
        // FIXME: when calling log() on this later, we need to handle this still being zero
        let mut sum_exp_dists = 0.0;
        let mut grad: SVector<Angstrom, 3> = self
            .inner
            .query_neighbors(pos)?
            .filter_map(|(_, atom)| {
                // FIXME: need to incorporate element radius into SDF and gradient
                // FIXME: which might mean we need some more mutable accumulators
                // let radius = atom.element.radius();
                let other: [Angstrom; 3] = atom.coords();
                let mut diff: SVector<Angstrom, 3> = Point::from(pos) - Point::from(other);
                let dist = diff.normalize_mut();

                if dist <= self.inner.info().cutoff() {
                    let exp_dist = (-dist).exp();
                    sum_exp_dists += exp_dist;
                    diff.scale_mut(exp_dist);
                    Some(diff)
                } else {
                    None
                }
            })
            .sum();

        grad /= sum_exp_dists;

        Some((-sum_exp_dists.ln(), grad.into()))
    }

    /// Computes the gradient and its norm at each "support" point, ie. each inner particle location.
    fn normals(&self) -> Vec<(Angstrom, [Angstrom; 3])> {
        vec![]
    }

    fn sdf<F: DualNumFloat, D: DualNum<F> + ComplexField<RealField = D> + Copy>(
        &self,
        x: SVector<D, 3>,
        neighbors: impl Iterator<Item = (F, SVector<D, 3>)>,
    ) -> D {
        let sum_exp_dists: D = neighbors
            .filter_map(|(_radius, other)| {
                // FIXME: need to incorporate element radius into SDF and gradient
                // FIXME: which might mean we need some more mutable accumulators
                // let radius = atom.element.radius();

                let diff: SVector<D, 3> = x - other;
                let dist: D = diff.norm();

                if let Some(cutoff) = F::from(self.inner.info().cutoff()) {
                    if dist.re() <= cutoff {
                        // this is correct in the sense
                        // that it approximates lim grad(·) for dist->0
                        // FIXME: maybe use Float::epsilon() here
                        if dist.re() != F::zero() {
                            Some((-dist).exp())
                        } else {
                            Some(F::one().into())
                            // FIXME: actually not sure what's more suitable here
                            // FIXME: above makes it actually continuous
                            // FIXME: below fits softmax -> hard max better for single neighbors?
                            // None
                        }
                    } else {
                        None
                    }
                } else {
                    None
                }
            })
            .sum();

        -sum_exp_dists.ln()
    }

    // FIXME: might have to change the sign somewhere to make it work with HMC
    fn harmonic_potential<D: DualNum<Angstrom> + ComplexField<RealField = D> + Copy>(
        &self,
        x: D,
        radius: D,
    ) -> D {
        const SPRING: Angstrom = 2.0;
        <Angstrom as Into<D>>::into(0.5 * SPRING) * (x - radius).powi(2)
    }

    pub fn gradient(&self, pos: [Angstrom; 3]) -> Option<(Angstrom, [Angstrom; 3])> {
        let neighbors = self.inner.query_neighbors(pos)?.map(|(_, atom)| {
            let coords: [Angstrom; 3] = atom.coords();
            let coords = coords.map(DsVec::from_re);
            (atom.element.radius(), SVector::from(coords))
        });
        let pos = SVector::from(pos);

        let (val, grad) = gradient(|x| self.sdf(x, neighbors), pos);

        Some((val.re(), grad.into()))
    }

    pub fn harmonic_gradient(&self, pos: [Angstrom; 3]) -> Option<(Angstrom, [Angstrom; 3])> {
        let neighbors = self.inner.query_neighbors(pos)?.map(|(_, atom)| {
            let coords: [Angstrom; 3] = atom.coords();
            let coords = coords.map(DsVec::from_re);
            (atom.element.radius(), SVector::from(coords))
        });
        let pos = SVector::from(pos);

        let (val, grad) = gradient(
            |x| self.harmonic_potential(self.sdf(x, neighbors), 1.05.into()),
            pos,
        );

        Some((val.re(), grad.into()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::Element;

    // TODO: adjust reference values for SDF with atom radii
    #[test]
    fn test_sdf_autodiff() {
        let points = vec![
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.5, 0.5, 0.5],
            [1.5, 1.5, 1.5],
        ];

        let reference_values = vec![
            -0.9259472,
            -0.9259472,
            -0.9259472,
            -0.9259472,
            -0.9259472,
            -0.9259472,
            -0.9259472,
            -1.0800674,
            -1.473609,
            -0.35109338,
        ];

        let reference_grads = vec![
            [-0.24194217, -0.24194217, -0.24194217],
            [-0.24194217, -0.24194217, 0.24194217],
            [-0.24194217, 0.24194217, -0.24194217],
            [0.24194217, -0.24194217, -0.24194217],
            [0.24194217, 0.24194217, -0.24194217],
            [-0.24194217, 0.24194217, 0.24194217],
            [0.24194217, -0.24194217, 0.24194217],
            [0.1249218, 0.1249218, 0.1249218],
            [6.8276282e-9, -0.0, -0.0],
            [0.17094302, 0.17094302, 0.17094302],
        ];

        // let reference_values = vec![
        //     -0.42150798,
        //     -0.42150798,
        //     -0.42150798,
        //     -0.42150798,
        //     -0.42150798,
        //     -0.42150798,
        //     -0.42150798,
        //     -0.6651994,
        //     -1.2134161,
        //     0.8660254,
        // ];

        // let reference_grads = vec![
        //     [-0.40066993, -0.40066993, -0.40066993],
        //     [-0.40066993, -0.40066993, 0.40066993],
        //     [-0.40066993, 0.40066993, -0.40066993],
        //     [0.40066993, -0.40066993, -0.40066993],
        //     [0.40066993, 0.40066993, -0.40066993],
        //     [-0.40066993, 0.40066993, 0.40066993],
        //     [0.40066993, -0.40066993, 0.40066993],
        //     [0.18915294, 0.18915294, 0.18915294],
        //     [8.856665e-9, -0.0, -0.0],
        //     [0.5773502, 0.5773502, 0.5773502],
        // ];

        let sdf = SmoothDistanceField {
            inner: CellGrid::new(
                points.iter().map(|&coords| Atom {
                    element: Element::default(),
                    coords,
                }),
                1.0,
            ),
        };

        let (sdf_values, sdf_grads): (Vec<Angstrom>, Vec<[Angstrom; 3]>) = points
            .iter()
            .copied()
            // .chain(std::iter::once([1.501; 3]))
            // .chain(std::iter::once([0.001; 3]))
            .filter_map(|atom| sdf.gradient(atom))
            .unzip();

        assert_eq!(reference_values, sdf_values);
        assert_eq!(reference_grads, sdf_grads);
    }

    #[test]
    fn test_harmonic_1d() {
        let positions = vec![
            -10.0f32, -8.0, -6.0, -4.0, -2.0, -1.0, 0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0,
        ];
        let radius = 0.0f32;

        let sdf = SmoothDistanceField {
            inner: CellGrid::new(std::iter::empty(), 1.0),
        };

        for pos in positions {
            let (_y, _dy) = first_derivative(|x| sdf.harmonic_potential(x, radius.into()), pos);
        }
    }
}
