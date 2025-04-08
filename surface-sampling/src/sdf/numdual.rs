use crate::Angstrom;
use crate::sdf::SmoothDistanceField;
use nalgebra::{ComplexField, Const, SVector};
use num_dual::*;
use zelll::Particle;

type DsVec<T> = DualVec<T, T, Const<3>>;
type Ds2Vec<T> = Dual2Vec<T, T, Const<3>>;

impl SmoothDistanceField {
    // TODO: this actually looks cleaner with a for loop...
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

    /// Returns the (approximate) smooth distance of `pos` to the internal point cloud and its gradient.
    ///
    /// If `pos` is too far away from the point cloud (ie. its neighborhood cannot be queried),
    /// `None` is returned.
    pub fn evaluate(&self, pos: [Angstrom; 3]) -> Option<(Angstrom, [Angstrom; 3])> {
        let neighbors = self.inner.query_neighbors(pos)?.map(|(_, atom)| {
            let coords: [Angstrom; 3] = atom.coords();
            let coords = coords.map(DsVec::from_re);
            // let coords = coords.map(Ds2Vec::from_re);
            (atom.element.radius(), SVector::from(coords))
        });
        let pos = SVector::from(pos);

        let (val, grad) = gradient(|x| self.sdf(x, neighbors), pos);
        // let (val, grad, hess) = hessian(|x| self.sdf(x, neighbors), pos);
        // FIXME: this is only the curvature at critical points
        // FIXME: to make it always work, I'd have to normalize the gradient and compute its Jacobian?
        // FIXME: cf. gauss map
        // let _det: Angstrom = hess.determinant();
        // FIXME: ie. like this?
        // FIXME: question is whether `grad` still contains the necessary information
        // let _det: Angstrom = jacobian(|x| x.normalize(), grad).1.determinant();
        // FIXME: normalizing makes determinant vanish (which I checked to make sense)
        // dbg!(_det);
        // FIXME: probably should impl global function (try_)curvature() based on (try_)gradient()
        // FIXME: which normalizes gradient, then computes jacobian of it before discarding
        // FIXME: dual derivative and then finally returns determinant

        Some((val.re(), grad.into()))
    }

    pub fn hmc_gradient(
        &self,
        pos: [Angstrom; 3],
        isoradius: Angstrom,
    ) -> Option<(Angstrom, [Angstrom; 3])> {
        let neighbors = self.inner.query_neighbors(pos)?.map(|(_, atom)| {
            let coords: [Angstrom; 3] = atom.coords();
            let coords = coords.map(DsVec::from_re);
            (atom.element.radius(), SVector::from(coords))
        });
        let pos = SVector::from(pos);

        let (val, grad) = gradient(
            |x| self.poly_potential(self.sdf(x, neighbors), isoradius.into()),
            pos,
        );

        Some((val.re(), grad.into()))
    }

    fn poly_potential<D: DualNum<Angstrom> + ComplexField<RealField = D> + Copy>(
        &self,
        x: D,
        radius: D,
    ) -> D {
        const KFORCE: Angstrom = 10.0;
        let offset_diff = x - radius + D::one();
        // linear term not strictly necessary since we have sdf as a potential as well
        D::from(KFORCE) * (offset_diff + offset_diff.powi(3) - offset_diff.powi(4))
        // D::from(FORCE) * (offset_diff - offset_diff.powi(3) - offset_diff.powi(4))
        // FIXME: actually using offset +1.0 with polynomial below works better but is arguably incorrect?
        // TODO: this works similarly well but is arguably not as nice
        // let offset_diff = x - radius + D::from(0.5);
        // D::from(FORCE) * (offset_diff - offset_diff.powi(3) - offset_diff.powi(4))
    }

    // fn harmonic_potential<D: DualNum<Angstrom> + ComplexField<RealField = D> + Copy>(
    //     &self,
    //     x: D,
    //     radius: D,
    // ) -> D {
    //     const KFORCE: Angstrom = 10.0;
    //     const INTERCEPT: Angstrom = 1.0;
    //     D::from(INTERCEPT) - D::from(KFORCE) * (x - radius).powi(2)
    // }
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
            -0.7998983683589524,
        ];

        let reference_grads = vec![
            [-0.276176313229217, -0.276176313229217, -0.276176313229217],
            [-0.276176313229217, -0.276176313229217, 0.276176313229217],
            [-0.276176313229217, 0.276176313229217, -0.276176313229217],
            [0.276176313229217, -0.276176313229217, -0.276176313229217],
            [0.276176313229217, 0.276176313229217, -0.276176313229217],
            [-0.276176313229217, 0.276176313229217, 0.276176313229217],
            [0.276176313229217, -0.276176313229217, 0.276176313229217],
            [
                0.14357909754235013,
                0.14357909754235013,
                0.14357909754235013,
            ],
            [-2.9256894966310793e-17, -0.0, -0.0],
            [
                0.21669568034989586,
                0.21669568034989586,
                0.21669568034989586,
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
