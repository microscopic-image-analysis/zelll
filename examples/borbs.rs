use kiss3d::camera::ArcBall;
use kiss3d::light::Light;
use kiss3d::nalgebra::{Point3, Vector3};
use kiss3d::window::Window;
use soa_derive::StructOfArray;
use zelll::cellgrid::*;

const NBORBS: usize = 20000;
const DELTA: f64 = 0.015;
const OUTER_RADIUS: f64 = 1.0;
const ALIGNMENT: f64 = 0.4;
const SEPARATION: f64 = 0.4;
const COHESION: f64 = 0.4;
//const STUBBORNNESS: f64 = 1.0;

//TODO: lessons from flamegraph:
//TODO: flattening in point_pairs() is expensive due to allocating Vecs for nested iterators
//TODO: par_point_pairs() in this naive implementation is inefficient

fn main() {
    let mut borbs = BorbVec::with_capacity(NBORBS);

    for _ in 0..NBORBS {
        borbs.push(Borb::new_random());
    }

    let mut window = Window::new_with_size("zelll borbs", 1920, 1080);
    window.set_light(Light::StickToCamera);
    window.set_framerate_limit(Some(60));
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
    let red = Point3::new(1.0, 0.0, 0.0);

    let mut cell_grid = CellGrid::new(&borbs.position, OUTER_RADIUS);

    //TODO: scale by number of neighbors for each boid
    let mut cohesion: Vec<Point3<f64>> = vec![Point3::default(); NBORBS];
    let mut separation: Vec<Vector3<f64>> = vec![Vector3::default(); NBORBS];
    let mut alignment: Vec<Vector3<f64>> = vec![Vector3::default(); NBORBS];
    let mut neighborhood: Vec<usize> = vec![0; NBORBS];

    while window.render_with_camera(&mut cam) {
        cell_grid.point_pairs().for_each(|(borb, other)| {
            let relative_pos = borbs.position[borb] - borbs.position[other];

            if relative_pos.norm() <= OUTER_RADIUS
            //&& borbs.direction[borb].angle(&relative_pos) < 1.7
            {
                neighborhood[borb] += 1;

                separation[borb] += relative_pos;
                cohesion[borb] += borbs.position[other].coords;
                alignment[borb] += borbs.direction[other];
            }
        });
        for (i, borb) in borbs.iter_mut().enumerate() {
            //TODO: proper error handling
            let borb_cohesion = cohesion[i] / neighborhood[i].max(1) as f64 - *borb.position;
            let borb_separation = separation[i];
            let borb_alignment = alignment[i] / neighborhood[i].max(1) as f64;

            *borb.direction += COHESION * borb_cohesion
                + SEPARATION * borb_separation
                + ALIGNMENT * borb_alignment;
        }

        for mut borb in borbs.iter_mut() {
            borb.fly(DELTA);
            window.draw_point(&borb.position.cast::<f32>(), &white);
            /*window.draw_line(
                &borb.position.cast::<f32>(),
                &(borb.position.cast::<f32>() + borb.direction.normalize().cast::<f32>()),
                &red,
            );*/
        }

        cell_grid = cell_grid.rebuild_mut(&borbs.position, None);

        cohesion.fill_with(Default::default);
        separation.fill_with(Default::default);
        alignment.fill_with(Default::default);
        neighborhood.fill_with(Default::default);
    }
}

#[derive(Debug, StructOfArray)]
#[soa_derive(Debug)]
pub struct Borb {
    position: Point3<f64>,
    direction: Vector3<f64>,
}

impl Borb {
    fn new_random() -> Self {
        Self {
            position: ((Vector3::new_random() - Vector3::new(0.5, 0.5, 0.5)) * 100.0).into(),
            direction: Vector3::new_random(),
        }
    }
}

impl BorbRefMut<'_> {
    fn fly(&mut self, delta: f64) {
        //self.direction.renormalize_fast();
        let magnitude = self.direction.norm();
        self.direction.normalize_mut();
        *self.direction *= magnitude.min(10.0);
        *self.position += *self.direction * delta;
    }
}

