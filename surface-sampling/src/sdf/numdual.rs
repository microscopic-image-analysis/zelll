use crate::Angstrom;
use crate::sdf::SmoothDistanceField;
use nalgebra::{ComplexField, SVector};
use num_dual::*;
use zelll::ParticleLike;

impl SmoothDistanceField {
    // TODO: this actually looks cleaner with a for loop...
    // FIXME: this should probably be done similar to LJTS?
    // FIXME: because having the hard cutoff affects the actual sdf values
    fn sdf<D: DualNum<Angstrom> + ComplexField<RealField = D> + Copy>(
        &self,
        x: SVector<D, 3>,
    ) -> Option<D> {
        let at: [D; 3] = x.into();
        let at = at.map(|coord| coord.re());
        let neighbors = self.inner.query_neighbors(at)?.map(|(_, atom)| {
            let coords: [Angstrom; 3] = atom.coords();
            let coords = coords.map(|c| c.into());
            (atom.element.radius(), SVector::from(coords))
        });

        let (scaled_exp_dists, atom_radii, total_exp_dists): (D, D, D) = neighbors
            .filter_map(|(radius, other)| {
                let cutoff = self.inner.info().cutoff();

                let diff = x - other;
                let dist = diff.norm();

                if dist.re() <= cutoff.into() {
                    // L2-norm is not continuous at zero
                    // so at zero, we handle it manually with this conditional
                    // FIXME: maybe use Float::epsilon() here
                    if dist.re() != 0.0 {
                        Some((
                            (-dist / radius).exp(),
                            (-dist).exp() * radius,
                            (-dist).exp(),
                        ))
                    } else {
                        // without this, num-dual would produce NaN gradients
                        Some((D::one(), radius.into(), D::one()))
                    }
                } else {
                    None
                }
            })
            .fold((D::zero(), D::zero(), D::zero()), |acc, curr| {
                let (scaled_exp_dists, atom_radii, total_exp_dists) = acc;
                let (scaled_exp_dist, atom_radius, exp_dist) = curr;
                (
                    scaled_exp_dists + scaled_exp_dist,
                    atom_radii + atom_radius,
                    total_exp_dists + exp_dist,
                )
            });

        // average atom radius in neighborhood
        let sigma = atom_radii / total_exp_dists;
        Some(-sigma * scaled_exp_dists.ln())
    }

    /// Returns the (approximate) smooth distance of `pos` to the internal point cloud and its gradient.
    ///
    /// If `pos` is too far away from the point cloud (ie. its neighborhood cannot be queried),
    /// `None` is returned.
    pub fn evaluate(&self, pos: [Angstrom; 3]) -> Option<(Angstrom, [Angstrom; 3])> {
        let (val, grad) = gradient(|x| self.sdf(x), &pos.into())?;
        Some((val.re(), grad.into()))
    }

    pub fn hmc_gradient(
        &self,
        pos: [Angstrom; 3],
        isoradius: Angstrom,
    ) -> Option<(Angstrom, [Angstrom; 3])> {
        let (val, grad) = gradient(
            |x| {
                self.sdf(x)
                    .map(|val| self.harmonic_potential(val, isoradius.into()))
            },
            &pos.into(),
        )?;

        Some((val.re(), grad.into()))
    }

    fn poly_potential<D: DualNum<Angstrom> + ComplexField<RealField = D> + Copy>(
        &self,
        x: D,
        radius: D,
    ) -> D {
        let offset_diff = x - radius + D::one();
        // linear term not strictly necessary since we have sdf as a potential as well
        D::from(self.k_force) * (offset_diff + offset_diff.powi(3) - offset_diff.powi(4))
    }

    fn harmonic_potential<D: DualNum<Angstrom> + ComplexField<RealField = D> + Copy>(
        &self,
        x: D,
        radius: D,
    ) -> D {
        -D::from(self.k_force) * (x - radius).powi(2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom::{Atom, Element};
    use zelll::CellGrid;

    // TODO: should use `approx` for this
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
            -2.012457244274712,
            -2.012457244274712,
            -2.012457244274712,
            -2.012457244274712,
            -2.012457244274712,
            -2.012457244274712,
            -2.012457244274712,
            -2.2994776285300675,
            -2.990326826730122,
            -0.7998983683589523,
        ];

        let reference_grads = vec![
            [
                -0.2761763132292168,
                -0.2761763132292168,
                -0.2761763132292168,
            ],
            [-0.2761763132292168, -0.2761763132292168, 0.2761763132292168],
            [-0.2761763132292168, 0.2761763132292168, -0.2761763132292168],
            [0.2761763132292168, -0.2761763132292168, -0.2761763132292168],
            [0.2761763132292168, 0.2761763132292168, -0.2761763132292168],
            [-0.2761763132292168, 0.2761763132292168, 0.2761763132292168],
            [0.2761763132292168, -0.2761763132292168, 0.2761763132292168],
            [
                0.14357909754235015,
                0.14357909754235015,
                0.14357909754235015,
            ],
            [6.651802279961878e-17, -0.0, -0.0],
            [
                0.21669568034989597,
                0.21669568034989597,
                0.21669568034989597,
            ],
        ];

        let sdf = SmoothDistanceField {
            inner: CellGrid::new(
                points.iter().map(|&coords| Atom {
                    element: Element::default(),
                    coords,
                }),
                1.0,
            ),
            surface_radius: 1.05,
            k_force: 10.0,
        };

        let (sdf_values, sdf_grads): (Vec<Angstrom>, Vec<[Angstrom; 3]>) = points
            .iter()
            .copied()
            // .chain(std::iter::once([1.501; 3]))
            // .chain(std::iter::once([0.001; 3]))
            .filter_map(|atom| sdf.evaluate(atom))
            .unzip();

        assert_eq!(&reference_values, &sdf_values);
        assert_eq!(&reference_grads, &sdf_grads);
    }
}
