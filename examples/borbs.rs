use zelll::cellgrid::*;

use kiss3d::camera::ArcBall;
use kiss3d::light::Light;
use kiss3d::nalgebra::{distance, Point3, UnitVector3, Vector3};
use kiss3d::window::Window;
use soa_derive::StructOfArray;

const NBORBS: usize = 10000;
const SPEED: f64 = 0.2;
const OUTER_RADIUS: f64 = 5.0;
const INNER_RADIUS: f64 = 2.0;
const ALIGNMENT: f64 = 0.05;
const SEPARATION: f64 = 0.7;

fn main() {
    let mut borbs = BorbVec::with_capacity(NBORBS);

    for _ in 0..NBORBS {
        borbs.push(Borb::new_random());
    }

    let mut window = Window::new_with_size("zelll borbs", 1920, 1080);
    window.set_light(Light::StickToCamera);
    // Doesn't work with my driver though
    window.set_point_size(5.0);

    let mut cam = ArcBall::new_with_frustrum(
        std::f32::consts::PI / 4.0,
        0.1,
        3072.0,
        Point3::new(500.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 0.0),
    );

    let white = Point3::new(1.0, 1.0, 1.0);

    let mut cell_grid = CellGrid::new(&borbs.position, OUTER_RADIUS);

    // average neighbor position (center of mass)
    // with weighted contributions for separation and cohesion
    let mut neighbor_com: Vec<Point3<f64>> = vec![Point3::default(); NBORBS];
    // neighbor direction sum (i.e. not normalized)
    let mut neighbor_dir: Vec<Vector3<f64>> = borbs
        .direction
        .iter()
        .map(|uvec| (*uvec).into_inner())
        .collect();

    while window.render_with_camera(&mut cam) {
        //TODO: this is definetly not correct yet
        //TODO: and it's ugly
        for (borb, other) in cell_grid.point_pairs() {
            let dist = distance(&borbs.position[borb], &borbs.position[other]);
            if dist <= OUTER_RADIUS {
                if dist < INNER_RADIUS {
                    neighbor_com[borb] -= borbs.position[other].coords * SEPARATION;
                    neighbor_dir[borb] -= borbs.direction[other].into_inner() * SEPARATION;
                } else {
                    neighbor_com[borb] += borbs.position[other].coords * (1.0 - SEPARATION);
                    neighbor_dir[borb] += borbs.direction[other].into_inner() * (1.0 - SEPARATION);
                }
            }
        }

        for (i, borb) in borbs.iter_mut().enumerate() {
            //TODO: proper error handling
            let new_direction = (1.0 - ALIGNMENT)
                * (neighbor_com[i] - *borb.position)
                    .try_normalize(f64::EPSILON)
                    .unwrap_or(**borb.direction)
                + ALIGNMENT
                    * neighbor_dir[i]
                        .try_normalize(f64::EPSILON)
                        .unwrap_or(**borb.direction);

            *borb.direction = UnitVector3::new_normalize(new_direction);
        }

        for mut borb in borbs.iter_mut() {
            borb.fly(SPEED);
            window.draw_point(&borb.position.cast::<f32>(), &white);
        }

        cell_grid = cell_grid.rebuild_mut(&borbs.position, None);

        neighbor_dir.fill_with(Default::default);
        neighbor_com.fill_with(Default::default);
    }
}

#[derive(Debug, StructOfArray)]
#[soa_derive(Debug)]
pub struct Borb {
    position: Point3<f64>,
    direction: UnitVector3<f64>,
}

impl Borb {
    fn new_random() -> Self {
        Self {
            position: ((Vector3::new_random() - Vector3::new(0.5, 0.5, 0.5)) * 100.0).into(),
            direction: UnitVector3::new_normalize(Vector3::new_random()),
        }
    }
}

impl BorbRefMut<'_> {
    fn fly(&mut self, speed: f64) {
        self.direction.renormalize_fast();
        *self.position += (*self.direction).into_inner() * speed;
    }
}

