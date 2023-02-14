use kiss3d::camera::ArcBall;
use kiss3d::light::Light;
use kiss3d::nalgebra::{distance, Point3, UnitVector3, Vector3};
use kiss3d::window::Window;
use soa_derive::StructOfArray;
use std::sync::{Arc, RwLock};
use zelll::cellgrid::*;

const NBORBS: usize = 25000;
const SPEED: f64 = 0.25;
const FOV: f64 = 1.7;
const OUTER_RADIUS: f64 = 6.0;
const ALIGNMENT: f64 = 0.85;
//TODO: mixed up some signs, behaves more like COHESION
const SEPARATION: f64 = 0.2;
const STUBBORNESS: f64 = 0.4;
const HOME: f64 = 0.04;

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
        Point3::new(350.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 0.0),
    );

    let white = Point3::new(1.0, 1.0, 1.0);

    let mut cell_grid = CellGrid::new(&borbs.position, OUTER_RADIUS);

    while window.render_with_camera(&mut cam) {
        // average neighbor position (center of mass)
        // with weighted contributions for separation and cohesion
        let mut neighbor_com: Vec<Arc<RwLock<Point3<f64>>>> = Vec::with_capacity(NBORBS);
        (0..NBORBS).for_each(|_| neighbor_com.push(Arc::new(RwLock::new(Point3::default()))));
        // neighbor direction sum (i.e. not normalized)
        let neighbor_dir: Vec<Arc<RwLock<Vector3<f64>>>> = borbs
            .direction
            .iter()
            .map(|uvec| Arc::new(RwLock::new((*uvec).into_inner())))
            .collect();
        //TODO: this is definetly not correct yet
        //TODO: and it's ugly
        cell_grid.par_point_pairs().for_each(|(borb, other)| {
            let dist = distance(&borbs.position[borb], &borbs.position[other]);
            let is_behind = borbs.direction[borb]
                .angle(&(&borbs.position[borb] - &borbs.position[other]))
                > FOV;
            if dist <= OUTER_RADIUS && !is_behind {
                let mut com_lock = neighbor_com[borb].write().unwrap();
                let mut dir_lock = neighbor_dir[borb].write().unwrap();
                if dist < OUTER_RADIUS * SEPARATION {
                    *com_lock -= borbs.position[other].coords;
                } else {
                    *com_lock += borbs.position[other].coords;
                }
                *dir_lock += borbs.direction[other].into_inner();
            }
        });
        for (i, borb) in borbs.iter_mut().enumerate() {
            //TODO: proper error handling
            let cohesion_separation = (*neighbor_com[i].read().unwrap() - *borb.position)
                .try_normalize(f64::EPSILON)
                .unwrap_or(Vector3::default());
            let alignment = neighbor_dir[i]
                .read()
                .unwrap()
                .try_normalize(f64::EPSILON)
                .unwrap_or(**borb.direction);

            let home = borb
                .position
                .coords
                .try_normalize(f64::EPSILON)
                .unwrap_or(Vector3::default());

            *borb.direction = UnitVector3::new_normalize(*borb.direction.slerp(
                &UnitVector3::new_normalize(
                    (1.0 - ALIGNMENT) * cohesion_separation + ALIGNMENT * alignment - HOME * home,
                ),
                1.0 - STUBBORNESS,
            ));
        }

        for mut borb in borbs.iter_mut() {
            borb.fly(SPEED);
            window.draw_point(&borb.position.cast::<f32>(), &white);
        }

        cell_grid = cell_grid.rebuild_mut(&borbs.position, None);

        //neighbor_dir.fill_with(Default::default);
        //neighbor_com.fill_with(Default::default);
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

