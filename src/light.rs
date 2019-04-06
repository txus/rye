use crate::color::Color;
use crate::linear::Point;
use crate::world::World;

pub trait Light {
    fn intensity_at(&self, p: &Point, w: &World) -> f32;
}

#[derive(Debug)]
pub struct PointLight {
    pub position: Point,
    pub intensity: Color,
}

impl Light for PointLight {
    fn intensity_at(&self, p: &Point, w: &World) -> f32 {
        if w.is_shadowed(&self.position, p) {
            0.0
        } else {
            1.0
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

#[cfg(test)]
mod tests {
    use super::*;

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
}
