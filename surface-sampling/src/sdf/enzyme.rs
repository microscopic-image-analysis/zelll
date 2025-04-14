// RUSTFLAGS="-C target-cpu=native -C target-feature=+avx2 -Zautodiff=Enable" cargo +enzyme run --example psssh --release --package psssh --features enzyme
// RUSTFLAGS="-Zautodiff=Enable" cargo +enzyme test --package psssh --features enzyme -- --show-output
use crate::Angstrom;
use crate::sdf::SmoothDistanceField;
use nalgebra::SVector;
use std::autodiff::autodiff;
use zelll::Particle;

impl SmoothDistanceField {
    /// Returns the (approximate) smooth distance of `pos` to the internal point cloud and its gradient.
    ///
    /// If `pos` is too far away from the point cloud (ie. its neighborhood cannot be queried),
    /// `None` is returned.
    pub fn evaluate(&self, pos: [Angstrom; 3]) -> Option<(Angstrom, [Angstrom; 3])> {
        let neighbors = self.inner.query_neighbors(pos)?.map(|(_, atom)| {
            let coords: [Angstrom; 3] = atom.coords();
            (atom.element.radius(), coords)
        });

        let mut grad = [0.0; 3];

        let val = d_sdf(&pos, &mut grad, &neighbors, self.inner.info().cutoff(), 1.0);

        Some((val, grad))
    }

    pub fn hmc_gradient(
        &self,
        pos: [Angstrom; 3],
        isoradius: Angstrom,
    ) -> Option<(Angstrom, [Angstrom; 3])> {
        let neighbors = self.inner.query_neighbors(pos)?.map(|(_, atom)| {
            let coords: [Angstrom; 3] = atom.coords();
            (atom.element.radius(), coords)
        });

        // TODO: for nuts_rs it would make sense to have hmc_gradient take &mut[f64] as gradient input
        // TODO: can use SVectorViewMut::from_slice() for this
        let mut grad = [0.0; 3];

        let val = d_hmc(
            &pos,
            &mut grad,
            &neighbors,
            self.inner.info().cutoff(),
            isoradius,
            1.0,
        );

        Some((val, grad.into()))
    }
}

#[autodiff(d_sdf, Reverse, Duplicated, Const, Const, Active)]
fn sdf(
    x: &[Angstrom; 3],
    neighbors: &(impl Iterator<Item = (Angstrom, [Angstrom; 3])> + Clone),
    cutoff: Angstrom,
) -> Angstrom {
    let mut scaled_exp_dists: Angstrom = 0.0;
    let mut atom_radii: Angstrom = 0.0;
    let mut total_exp_dists: Angstrom = 0.0;

    for (radius, atom) in neighbors.clone() {
        let this = SVector::from(*x);
        let other = SVector::from(atom);
        let dist = (this - other).norm();

        if dist <= cutoff {
            scaled_exp_dists += (-dist / radius).exp();
            atom_radii += (-dist).exp() * radius;
            total_exp_dists += (-dist).exp();
        }
    }

    // average atom radius in neighborhood
    let sigma = atom_radii / total_exp_dists;

    -sigma * scaled_exp_dists.ln()
}

#[autodiff(d_hmc, Reverse, Duplicated, Const, Const, Const, Active)]
fn hmc_potential(
    x: &[Angstrom; 3],
    neighbors: &(impl Iterator<Item = (Angstrom, [Angstrom; 3])> + Clone),
    cutoff: Angstrom,
    isoradius: Angstrom,
) -> Angstrom {
    let dist = sdf(x, neighbors, cutoff);
    poly_potential(dist, isoradius)
}

fn poly_potential(x: Angstrom, radius: Angstrom) -> Angstrom {
    // const KFORCE: Angstrom = 10.0;
    let offset_diff = x - radius + 1.0;
    self.k_force * (offset_diff + offset_diff.powi(3) - offset_diff.powi(4))
    // let offset_diff = x - radius + 0.5;
    // KFORCE * (offset_diff - offset_diff.powi(3) - offset_diff.powi(4))
    // FIXME: cf. numdual.rs
    // FIXME: actually using offset +1.0 with second polynomial works better?
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
            -2.9903268267301217,
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
            .filter_map(|atom| sdf.evaluate(atom))
            .unzip();

        assert_eq!(&reference_values, &sdf_values);
        assert_eq!(&reference_grads, &sdf_grads);
    }
}
