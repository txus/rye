use crate::color::Color;
use crate::light::PointLight;
use crate::linear::{Point, Vector};
use crate::patterns::Pattern;
use crate::shapes::Shape;

pub struct Material {
    pub color: Color,
    pub ambient: f32,
    pub diffuse: f32,
    pub specular: f32,
    pub shininess: f32,
    pub reflective: f32,
    pub transparency: f32,
    pub refractive_index: f32,
    pub pattern: Option<Box<Pattern>>,
}

const BASE: Material = Material {
    color: Color { red: 1.0, green: 1.0, blue: 1.0 },
    ambient: 0.1,
    diffuse: 0.9,
    specular: 0.9,
    shininess: 200.0,
    reflective: 0.0,
    transparency: 0.0,
    refractive_index: 1.0,
    pattern: None,
};

impl Material {
    pub fn base() -> &'static Material {
        &BASE
    }

    pub fn default() -> Material {
        Material {
            color: Color { red: 1.0, green: 1.0, blue: 1.0 },
            ambient: 0.1,
            diffuse: 0.9,
            specular: 0.9,
            shininess: 200.0,
            reflective: 0.0,
            transparency: 0.0,
            refractive_index: 1.0,
            pattern: None,
        }
    }

    pub fn lighting(
        &self,
        object: &Box<Shape>,
        light: &PointLight,
        position: &Point,
        eye: &Vector,
        normal: &Vector,
        intensity: f32,
    ) -> Color {
        let color = if let Some(pattern) = &self.pattern {
            pattern.color_at_object(&object, &position)
        } else {
            self.color
        };
        let effective_color = color * light.intensity;
        let ambient = effective_color * self.ambient;

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
        ambient + ((diffuse + specular) * intensity)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::patterns::StripePattern;
    use crate::shapes::Sphere;

    #[test]
    fn lighting_eye_between_light_and_surface() {
        let mat = Material::default();
        let s: Box<Shape> = Box::from(Sphere::new());
        let position = Point::origin();
        let eye = Vector::new(0.0, 0.0, -1.0);
        let normal = Vector::new(0.0, 0.0, -1.0);
        let light = PointLight::new(Point::new(0.0, 0.0, -10.0), Color::white());
        let result = mat.lighting(&s, &light, &position, &eye, &normal, 1.0);
        assert_eq!(result, Color::new(1.9, 1.9, 1.9));
    }

    #[test]
    fn lighting_eye_between_light_and_surface_offset_45_deg() {
        let mat = Material::default();
        let s: Box<Shape> = Box::from(Sphere::new());
        let position = Point::origin();
        let eye = Vector::new(0.0, 2_f32.sqrt() / 2.0, -2_f32.sqrt() / 2.0);
        let normal = Vector::new(0.0, 0.0, -1.0);
        let light = PointLight::new(Point::new(0.0, 0.0, -10.0), Color::white());
        let result = mat.lighting(&s, &light, &position, &eye, &normal, 1.0);
        assert_eq!(result, Color::white());
    }

    #[test]
    fn lighting_eye_opposite_surface_light_offset_45_deg() {
        let mat = Material::default();
        let s: Box<Shape> = Box::from(Sphere::new());
        let position = Point::origin();
        let eye = Vector::new(0.0, 0.0, -1.0);
        let normal = Vector::new(0.0, 0.0, -1.0);
        let light = PointLight::new(Point::new(0.0, 10.0, -10.0), Color::white());
        let result = mat.lighting(&s, &light, &position, &eye, &normal, 1.0);
        assert_eq!(result, Color::new(0.7364, 0.7364, 0.7364));
    }

    #[test]
    fn lighting_eye_path_of_reflection() {
        let mat = Material::default();
        let s: Box<Shape> = Box::from(Sphere::new());
        let position = Point::origin();
        let eye = Vector::new(0.0, -2_f32.sqrt() / 2.0, -2_f32.sqrt() / 2.0);
        let normal = Vector::new(0.0, 0.0, -1.0);
        let light = PointLight::new(Point::new(0.0, 10.0, -10.0), Color::white());
        let result = mat.lighting(&s, &light, &position, &eye, &normal, 1.0);
        assert_eq!(result, Color::new(1.6364, 1.6364, 1.6364));
    }

    #[test]
    fn lighting_behind_surface() {
        let mat = Material::default();
        let s: Box<Shape> = Box::from(Sphere::new());
        let position = Point::origin();
        let eye = Vector::new(0.0, 0.0, -1.0);
        let normal = Vector::new(0.0, 0.0, -1.0);
        let light = PointLight::new(Point::new(0.0, 0.0, 10.0), Color::white());
        let result = mat.lighting(&s, &light, &position, &eye, &normal, 1.0);
        assert_eq!(result, Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn lighting_with_surface_in_shadow() {
        let mat = Material::default();
        let s: Box<Shape> = Box::from(Sphere::new());
        let position = Point::origin();
        let eye = Vector::new(0.0, 0.0, -1.0);
        let normal = Vector::new(0.0, 0.0, -1.0);
        let light = PointLight::new(Point::new(0.0, 0.0, -10.0), Color::white());
        let result = mat.lighting(&s, &light, &position, &eye, &normal, 0.0);
        assert_eq!(result, Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn lighting_with_pattern_applied() {
        let s: Box<Shape> = Box::from(Sphere::new());
        let mut mat = Material::default();
        mat.pattern = Some(Box::from(StripePattern::new(
            Color::white(),
            Color::black(),
        )));
        mat.ambient = 1.0;
        mat.diffuse = 0.0;
        mat.specular = 0.0;
        let eye = Vector::new(0.0, 0.0, -1.0);
        let normal = Vector::new(0.0, 0.0, -1.0);
        let light = PointLight::new(Point::new(0.0, 0.0, 10.0), Color::white());
        let c1 = mat.lighting(&s, &light, &Point::new(0.9, 0.0, 0.0), &eye, &normal, 1.0);
        let c2 = mat.lighting(&s, &light, &Point::new(1.1, 0.0, 0.0), &eye, &normal, 1.0);
        assert_eq!(c1, Color::white());
        assert_eq!(c2, Color::black());
    }
}
