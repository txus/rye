use crate::color::Color;
use crate::linear::{Point, Vector};
use crate::world::World;
use crate::jitter::RandomJitter;
use crate::rays::Precomputation;
use crate::materials::Material;

use std::cell::RefCell;
use std::rc::Rc;
use std::ops::DerefMut;

pub trait Light {
    fn intensity(&self) -> &Color;
    fn intensity_at(&self, color: &Color, comps: &Precomputation, mat: &Material, w: &World) -> Color;
}

#[derive(Debug)]
pub struct PointLight {
    pub position: Point,
    pub intensity: Color,
}

impl Light for PointLight {
    fn intensity(&self) -> &Color { &self.intensity }
    fn intensity_at(&self, color: &Color, comps: &Precomputation, mat: &Material, w: &World) -> Color {
        if w.is_shadowed(&self.position, &comps.over_point) {
            Color::black()
        } else {
            let lightv = (self.position - comps.over_point).normalize();
            mat.lighten_hit(&color, &lightv, &self.intensity, &comps)
        }
    }
}

impl PointLight {
    pub fn new(position: Point, intensity: Color) -> PointLight {
        PointLight {
            position,
            intensity,
        }
    }
}

#[derive(Debug)]
pub struct AreaLight {
    pub corner: Point,
    pub uvec: Vector,
    pub usteps: usize,
    pub vvec: Vector,
    pub vsteps: usize,
    pub samples: usize,
    pub intensity: Color,
    jitter: Rc<RefCell<RandomJitter>>,
}

impl Light for AreaLight {
    fn intensity(&self) -> &Color { &self.intensity }
    fn intensity_at(&self, color: &Color, comps: &Precomputation, mat: &Material, w: &World) -> Color {
        let jitter = self.jitter.clone();

        let mut result = Color::black();

        for v in 0..self.vsteps {
            for u in 0..self.usteps {
                let mut the_jitter = jitter.borrow_mut();
                let mut the = the_jitter.deref_mut();
                let light_position = self.point_on_light(&mut the, u, v);
                if !w.is_shadowed(&light_position, &comps.over_point) {
                    let lightv = (light_position - comps.point).normalize();
                    result = result + mat.lighten_hit(&color, &lightv, &self.intensity, &comps);
                }
            }
        }

        result / (self.samples as f32)
    }
}

impl AreaLight {
    pub fn new(corner: Point, full_uvec: Vector, usteps: usize, full_vvec: Vector, vsteps: usize, intensity: Color) -> AreaLight {
        let jitter = Rc::new(RefCell::new(RandomJitter::new()));

        AreaLight {
            jitter,
            corner,
            uvec: full_uvec / (usteps as f32),
            usteps,
            vvec: full_vvec / (vsteps as f32),
            vsteps,
            intensity,
            samples: usteps * vsteps,
        }
    }

    fn point_on_light<I: Iterator<Item=f32>>(&self, jitter: &mut I, u: usize, v: usize) -> Point {
        self.corner + self.uvec * (u as f32 + jitter.next().unwrap()) + self.vvec * (v as f32 + jitter.next().unwrap())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::jitter;

    #[test]
    fn initialize() {
        let pos = Point::origin();
        let intensity = Color::white();
        let light = PointLight::new(pos, intensity);
        assert_eq!(light.position, pos);
        assert_eq!(light.intensity, intensity);
    }

    #[test]
    fn point_lights_evaluate_light_intensity_at_a_given_point() {
        let w = World::default();
        let light = &w.light_source;
        let examples: Vec<(Point, f32)> = vec![
        (Point::new(0.0, 1.0001, 0.0)   , 1.0 ),
        (Point::new(-1.0001, 0.0, 0.0)  , 1.0 ),
        (Point::new(0.0, 0.0, -1.0001), 1.0 ),
        (Point::new(0.0, 0.0, 1.0001)   , 0.0 ),
        (Point::new(1.0001, 0.0, 0.0)   , 0.0 ),
        (Point::new(0.0, -1.0001, 0.0)  , 0.0 ),
        (Point::new(0.0, 0.0, 0.0)      , 0.0 ),
        ];

        for (p, intensity) in &examples {
            assert_eq!(light.intensity_at(p, &w), *intensity);
        }
    }

    #[test]
    fn finding_single_point_on_area_light() {
        let corner = Point::origin();
        let v1 = Vector::new(2.0, 0.0, 0.0);
        let v2 = Vector::new(0.0, 0.0, 1.0);
        let light = AreaLight::new(corner, v1, 4, v2, 2, Color::white());

        let examples: Vec<(usize, usize, Point)> = vec![
            (0, 0, Point::new(0.25, 0.0, 0.25)),
            (1, 0, Point::new(0.75, 0.0, 0.25)),
            (0, 1, Point::new(0.25, 0.0, 0.75)),
            (2, 0, Point::new(1.25, 0.0, 0.25)),
            (3, 1, Point::new(1.75, 0.0, 0.75)),
        ];

        let jitter = jitter::ConstantJitter::new(0.5);

        for (u, v, result) in &examples {
            assert_eq!(light.point_on_light(jitter, *u, *v), *result);
        }
    }

    #[test]
    fn area_light_intensity_function() {
        let w = World::default();
        let corner = Point::new(-0.5, -0.5, -5.0);
        let v1 = Vector::new(1.0, 0.0, 0.0);
        let v2 = Vector::new(0.0, 1.0, 0.0);
        let light = AreaLight::new(corner, v1, 2, v2, 2, Color::white());

        let examples: Vec<(Point, f32)> = vec![
            (Point::new(0.0, 0.0, 2.0), 0.0),
            (Point::new(1.0, -1.0, 2.0), 0.25),
            (Point::new(1.5, 0.0, 2.0), 0.5),
            (Point::new(1.25, 1.25, 3.0), 0.75),
            (Point::new(0.0, 0.0, -2.0), 1.0),
        ];

        for (point, result) in &examples {
            assert_eq!(light.intensity_at(point, &w), *result);
        }
    }
}
