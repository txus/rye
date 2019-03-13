use crate::light::PointLight;
use crate::primitives::{Color, Point, Vector};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Material {
    pub color: Color,
    pub ambient: f32,
    pub diffuse: f32,
    pub specular: f32,
    pub shininess: f32,
}

impl Material {
    pub fn new(
        color: Color,
        ambient: f32,
        diffuse: f32,
        specular: f32,
        shininess: f32,
    ) -> Material {
        Material {
            color,
            ambient,
            diffuse,
            specular,
            shininess,
        }
    }

    pub fn default() -> Material {
        Self::new(Color::new(1.0, 1.0, 1.0), 0.1, 0.9, 0.9, 200.0)
    }

    pub fn lighting(
        &self,
        light: &PointLight,
        position: &Point,
        eye: &Vector,
        normal: &Vector,
        in_shadow: bool,
    ) -> Color {
        let effective_color = self.color * light.intensity;
        let ambient = effective_color * self.ambient;

        if in_shadow {
            return ambient;
        }

        let lightv = (light.position - *position).normalize();
        let light_dot_normal = lightv.dot(normal);
        let diffuse: Color;
        let specular: Color;

        if light_dot_normal < 0.0 {
            diffuse = Color::black();
            specular = Color::black();
        } else {
            diffuse = effective_color * self.diffuse * light_dot_normal;
            let reflectv = (-lightv).reflect(normal);
            let reflect_dot_eye = reflectv.dot(eye);
            if reflect_dot_eye <= 0.0 {
                specular = Color::black();
            } else {
                let factor = reflect_dot_eye.powf(self.shininess);
                specular = light.intensity * self.specular * factor;
            }
        }
        ambient + diffuse + specular
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{Color, Point, Vector};

    #[test]
    fn lighting_eye_between_light_and_surface() {
        let mat = Material::default();
        let position = Point::new(0.0, 0.0, 0.0);
        let eye = Vector::new(0.0, 0.0, -1.0);
        let normal = Vector::new(0.0, 0.0, -1.0);
        let light = PointLight::new(Point::new(0.0, 0.0, -10.0), Color::new(1.0, 1.0, 1.0));
        let result = mat.lighting(&light, &position, &eye, &normal, false);
        assert_eq!(result, Color::new(1.9, 1.9, 1.9));
    }

    #[test]
    fn lighting_eye_between_light_and_surface_offset_45_deg() {
        let mat = Material::default();
        let position = Point::new(0.0, 0.0, 0.0);
        let eye = Vector::new(0.0, 2_f32.sqrt() / 2.0, -2_f32.sqrt() / 2.0);
        let normal = Vector::new(0.0, 0.0, -1.0);
        let light = PointLight::new(Point::new(0.0, 0.0, -10.0), Color::new(1.0, 1.0, 1.0));
        let result = mat.lighting(&light, &position, &eye, &normal, false);
        assert_eq!(result, Color::new(1.0, 1.0, 1.0));
    }

    #[test]
    fn lighting_eye_opposite_surface_light_offset_45_deg() {
        let mat = Material::default();
        let position = Point::new(0.0, 0.0, 0.0);
        let eye = Vector::new(0.0, 0.0, -1.0);
        let normal = Vector::new(0.0, 0.0, -1.0);
        let light = PointLight::new(Point::new(0.0, 10.0, -10.0), Color::new(1.0, 1.0, 1.0));
        let result = mat.lighting(&light, &position, &eye, &normal, false);
        assert_eq!(result, Color::new(0.7364, 0.7364, 0.7364));
    }

    #[test]
    fn lighting_eye_path_of_reflection() {
        let mat = Material::default();
        let position = Point::new(0.0, 0.0, 0.0);
        let eye = Vector::new(0.0, -2_f32.sqrt() / 2.0, -2_f32.sqrt() / 2.0);
        let normal = Vector::new(0.0, 0.0, -1.0);
        let light = PointLight::new(Point::new(0.0, 10.0, -10.0), Color::new(1.0, 1.0, 1.0));
        let result = mat.lighting(&light, &position, &eye, &normal, false);
        assert_eq!(result, Color::new(1.6364, 1.6364, 1.6364));
    }

    #[test]
    fn lighting_behind_surface() {
        let mat = Material::default();
        let position = Point::new(0.0, 0.0, 0.0);
        let eye = Vector::new(0.0, 0.0, -1.0);
        let normal = Vector::new(0.0, 0.0, -1.0);
        let light = PointLight::new(Point::new(0.0, 0.0, 10.0), Color::new(1.0, 1.0, 1.0));
        let result = mat.lighting(&light, &position, &eye, &normal, false);
        assert_eq!(result, Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn lighting_with_surface_in_shadow() {
        let mat = Material::default();
        let position = Point::new(0.0, 0.0, 0.0);
        let eye = Vector::new(0.0, 0.0, -1.0);
        let normal = Vector::new(0.0, 0.0, -1.0);
        let light = PointLight::new(Point::new(0.0, 0.0, -10.0), Color::new(1.0, 1.0, 1.0));
        let in_shadow = true;
        let result = mat.lighting(&light, &position, &eye, &normal, in_shadow);
        assert_eq!(result, Color::new(0.1, 0.1, 0.1));
    }
}
