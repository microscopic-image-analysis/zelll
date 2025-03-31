use crate::Angstrom;
use crate::atom::Atom;
use nalgebra::{ComplexField, Const, SVector};
use num_dual::*;
use zelll::{CellGrid, Particle};

type DsVec<T> = DualVec<T, T, Const<3>>;

pub struct SmoothDistanceField {
    inner: CellGrid<Atom, 3, Angstrom>,
}

impl SmoothDistanceField {
    /// Computes the gradient and its norm at each "support" point, ie. each inner particle location.
    fn normals(&self) -> Vec<(Angstrom, [Angstrom; 3])> {
        todo!()
    }

    fn sdf<F: DualNumFloat, D: DualNum<F> + ComplexField<RealField = D> + Copy>(
        &self,
        x: SVector<D, 3>,
        neighbors: impl Iterator<Item = (F, SVector<D, 3>)>,
    ) -> D {
        let (scaled_exp_dists, atom_radii, total_exp_dists): (D, D, D) = neighbors
            .filter_map(|(radius, other)| {
                let cutoff = F::from(self.inner.info().cutoff())?;

                let diff = x - other;
                let dist = diff.norm();

                if dist.re() <= cutoff {
                    // L2-norm is not continuous at zero
                    // so at zero, we handle it manually with this conditional
                    // FIXME: maybe use Float::epsilon() here
                    if dist.re() != F::zero() {
                        Some(((-dist / radius).exp(), -dist.exp() * radius, -dist.exp()))
                    } else {
                        let one = F::one().into();
                        // without this, num-dual would produce NaN gradients
                        Some((one, radius.into(), one))
                    }
                } else {
                    None
                }
            })
            .fold(
                (F::zero().into(), F::zero().into(), F::zero().into()),
                |acc, curr| {
                    let (scaled_exp_dists, atom_radii, total_exp_dists) = acc;
                    let (scaled_exp_dist, atom_radius, exp_dist) = curr;
                    (
                        scaled_exp_dists + scaled_exp_dist,
                        atom_radii + atom_radius,
                        total_exp_dists + exp_dist,
                    )
                },
            );

        // average atom radius in neighborhood
        let sigma = atom_radii / total_exp_dists;

        -sigma * scaled_exp_dists.ln()
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

    /// Returns the (approximate) smooth distance of `pos` to the internal point cloud and its gradient.
    ///
    /// If `pos` is too far away from the point cloud (ie. its neighborhood cannot be queried),
    /// `None` is returned.
    pub fn evaluate(&self, pos: [Angstrom; 3]) -> Option<(Angstrom, [Angstrom; 3])> {
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
            -2.0124574, -2.0124574, -2.0124574, -2.0124574, -2.0124574, -2.0124574, -2.0124574,
            -2.2994778, -2.990327, -0.7998983,
        ];

        let reference_grads = vec![
            [-0.27617636, -0.27617636, -0.27617636],
            [-0.27617636, -0.27617636, 0.27617636],
            [-0.27617636, 0.27617636, -0.27617636],
            [0.27617636, -0.27617636, -0.27617636],
            [0.27617636, 0.27617636, -0.27617636],
            [-0.27617636, 0.27617636, 0.27617636],
            [0.27617636, -0.27617636, 0.27617636],
            [0.14357911, 0.14357911, 0.14357911],
            [1.5707174e-8, -0.0, -0.0],
            [0.21669577, 0.21669577, 0.21669577],
        ];

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
            .filter_map(|atom| sdf.evaluate(atom))
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
