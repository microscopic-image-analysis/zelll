use dashmap::DashMap;
use itertools::Itertools;
use kiss3d::camera::ArcBall;
use kiss3d::light::Light;
use kiss3d::nalgebra::{Point3, Vector3};
use kiss3d::window::Window;
use rand::distributions::Standard;
use rand::prelude::*;
use soa_derive::StructOfArray;
use zelll::cellgrid::*;

const NBEADS: usize = 2000;
const TDELTA: f64 = 0.015;
const MAX_SPEED: f64 = 10.0;
//TODO: this is for the cellgrid, should LJTS use this or rather 2.3 * SIGMA?
const CUTOFF: f64 = 5.0;
// Harmonic potential constants
const LINK_LENGTH: f64 = 1.0;
const KSPRING: f64 = 1000.0;
const BEADMASS: f64 = 1.0;
// Lennard-Jones potential constants
const EPSILON: f64 = 0.1;
const SIGMA: f64 = 0.01;
//TODO: This is a little hack since most float operations (e.g. f64::powi()) are not (yet) `const`
//TODO: and I don't see why I should use OnceCell for this:
const _SIGMA_POW6: f64 = SIGMA * SIGMA * SIGMA * SIGMA * SIGMA * SIGMA;
const _SIGMA_POW12: f64 = _SIGMA_POW6 * _SIGMA_POW6;
const LJA: f64 = 4.0 * EPSILON * _SIGMA_POW12;
const LJB: f64 = 4.0 * EPSILON * _SIGMA_POW6;

// TODO: SoAIndex::get() conflicts with Itertools::get() (>= "0.13")
// TODO: there's nobody really at fault but I think
// TODO: the easiest fix would be to use fully qualified syntax for SoAIndex::get()
// TODO: where it's used inside of that proc_macro
// TODO: probably should file an issue? or fix it myself, it's relatively simple
// TODO: see https://github.com/lumol-org/soa-derive/blob/master/soa-derive-internal/src/index.rs
// TODO: the reason is: some Range* types impl Iterator
#[derive(Debug, StructOfArray)]
#[soa_derive(Debug)]
pub struct Bead {
    position: Point3<f64>,
    velocity: Vector3<f64>,
}

impl Bead {
    fn new_random() -> Self {
        Self {
            position: ((Vector3::from_iterator(thread_rng().sample_iter(Standard))
                - Vector3::new(0.5, 0.5, 0.5))
                * 100.0)
                .into(),
            velocity: Vector3::new(0.0, 0.0, 0.0),
        }
    }
}

// This implements 3 out of 4 steps of velocity-verlet integration *with* half-step velocites
// (equivalent to kick-drift-kick leap frog integration)
// Using half-step velocities because this way I don't need to save the acceleration of the previous time step
// The 3rd step (deriving new acceleration from force-field interactions) needs to be done separately
impl BeadRefMut<'_> {
    // 1st and 4th step of velocity-verlet
    // needs to be called twice with intermittent update of the derived accelaration using the force field
    fn update_velocity_halfstep(&mut self, acceleration: Vector3<f64>) {
        *self.velocity += 0.5 * acceleration * TDELTA;
        *self.velocity = self.velocity.cap_magnitude(MAX_SPEED);
    }
    // 2nd step of velocity-verlet
    fn update_position(&mut self) {
        *self.position += *self.velocity * TDELTA;
    }
}

/// potentials are using global constants for now
mod pairpotentials {
    /// Lennard-Jones potential
    #[inline]
    pub fn lennard_jones(radius: f64) -> f64 {
        12.0 * super::LJA / radius.powi(13) - 6.0 * super::LJB / radius.powi(7)
    }

    /// Lennard Jones, truncated and shifted
    /// EPSILON, SIGMA are global constants currently
    #[inline]
    pub fn lennard_jones_ts(radius: f64, cutoff: f64) -> f64 {
        if radius < cutoff {
            lennard_jones(radius) - lennard_jones(cutoff)
        } else {
            0.0
        }
    }
    /// Harmonic Oscillator with Bond Stretching
    /// spring constant is a global constant currently
    #[inline]
    pub fn harmonic_oscillator(radius: f64, bond_length: f64) -> f64 {
        super::KSPRING * (radius - bond_length)
    }
}

fn main() {
    let mut window = Window::new_with_size("[zelll] beads-on-a-chain", 1920, 1080);
    window.set_light(Light::StickToCamera);
    window.set_framerate_limit(Some(60));

    let mut cam = ArcBall::new_with_frustrum(
        std::f32::consts::PI / 4.0,
        0.1,
        3072.0,
        Point3::new(40.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 0.0),
    );

    let white = Point3::new(1.0, 1.0, 1.0);
    let red = Point3::new(1.0, 0.0, 0.0);

    let mut beads = BeadVec::with_capacity(NBEADS);

    for _ in 0..NBEADS {
        beads.push(Bead::new_random());
    }

    let accelerations: DashMap<usize, Vector3<f64>> = DashMap::with_capacity(NBEADS);

    let mut cell_grid = CellGrid::new(beads.position.iter().map(|p| p.coords.as_ref()), CUTOFF);

    while window.render_with_camera(&mut cam) {
        //TODO: unfortunately, soa_derive doesn't support rayon directly
        for (i, mut bead) in beads.iter_mut().enumerate() {
            if let Some(acceleration) = accelerations.insert(i, Vector3::new(0.0, 0.0, 0.0)) {
                bead.update_velocity_halfstep(acceleration);
            }
            bead.update_position();
        }

        //TODO: this is the 3rd verlet step
        // 3a. bonded interactions:
        for i in 0..NBEADS - 1 {
            let mut acc = beads.position[i + 1] - beads.position[i];
            let radius = acc.norm();
            let magnitude = pairpotentials::harmonic_oscillator(radius, LINK_LENGTH);
            acc.set_magnitude(magnitude);
            accelerations.alter(&i, |_, a| a + acc);
            accelerations.alter(&(i + 1), |_, a| a - acc);
        }
        // 3b. non-bonded interactions:
        #[cfg(not(feature = "rayon"))]
        cell_grid.for_each_point_pair(|bead, other| {
            let mut acc = beads.position[other] - beads.position[bead];
            let radius = acc.norm();

            if radius <= CUTOFF {
                let magnitude = pairpotentials::lennard_jones_ts(radius, CUTOFF);

                acc.set_magnitude(magnitude);
                accelerations.alter(&bead, |_, a| a + acc);
                accelerations.alter(&other, |_, a| a - acc);
            }
        });
        #[cfg(feature = "rayon")]
        cell_grid.par_for_each_point_pair(|bead, other| {
            let mut acc = beads.position[other] - beads.position[bead];
            let radius = acc.norm();

            if radius <= CUTOFF {
                let magnitude = pairpotentials::lennard_jones_ts(radius, CUTOFF);

                acc.set_magnitude(magnitude);
                accelerations.alter(&bead, |_, a| a + acc);
                accelerations.alter(&other, |_, a| a - acc);
            }
        });
        //TODO: 4th verlet step: again half-step velocity update
        for (i, mut bead) in beads.iter_mut().enumerate() {
            if let Some(acceleration) = accelerations.get(&i) {
                bead.update_velocity_halfstep(*acceleration);
            }
            window.draw_point(&bead.position.cast::<f32>(), &white);
        }

        for (a, b) in beads.position.iter().tuple_windows() {
            window.draw_line(&a.cast::<f32>(), &b.cast::<f32>(), &red);
        }

        cell_grid.rebuild_mut(beads.position.iter().map(|p| p.coords.as_ref()), None);
    }
}
