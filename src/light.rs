use crate::primitives::{Color, Point};

pub trait Light
where
    Self: std::fmt::Debug + PartialEq,
{
}

#[derive(Debug, PartialEq)]
pub struct PointLight {
    pub position: Point,
    pub intensity: Color,
}

impl Light for PointLight {}

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
        let pos = Point::new(0.0, 0.0, 0.0);
        let intensity = Color::new(1.0, 1.0, 1.0);
        let light = PointLight::new(pos, intensity);
        assert_eq!(light.position, pos);
        assert_eq!(light.intensity, intensity);
    }
}
