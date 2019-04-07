use crate::color::Color;
use crate::light::Light;
use crate::linear::{Vector, EPSILON, Point};
use crate::patterns::Pattern;
use crate::shapes::Shape;
use crate::world::World;
use crate::rays::Precomputation;

pub struct Material {
    pub color: Color,
    pub ambient: f32,
    pub diffuse: f32,
    pub specular: f32,
    pub shininess: f32,
    pub reflective: f32,
    pub transparency: f32,
    pub roughness: f32,
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
    roughness: 0.0,
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
            roughness: 0.0,
            refractive_index: 1.0,
            pattern: None,
        }
    }

    // OrenNayar implements the Van Ouwerkerks rewrite of Oren-Nayar model,
    // see: http://shaderjvo.blogspot.com/2011/08/van-ouwerkerks-rewrite-of-oren-nayar.html
    // In this model sigma represents the material roughness,
    // if sigma = 0 then the material is not rough at all and the formula simplifies to the Lambert model
    fn oren_nayar(eye: &Vector, light: &Vector, normal: &Vector, sigma: f32) -> f32 {
        if sigma == 0.0 {
            1.0
        } else {
            let sigma2 = sigma.powi(2);
            let a = 1.0 - (sigma2 / (2.0 * (sigma2 + 0.33)));
            let b = 0.45 * sigma2 / (sigma2 + 0.09);

            let l = normal.dot(&light).max(0.0); //cosThetaI
            let v = normal.dot(&eye).max(0.0); //cosTheta0

            if l >= (1.0 - EPSILON) || v >= (1.0 - EPSILON) {
                //cosPhi and sinTheta will be zero, exit now
                a
            } else {
                let n_l = *normal * l;
                let light_plane = (*light - n_l).normalize();
                let view_plane = (*eye - n_l).normalize();

                let p = view_plane.dot(&light_plane).max(0.0); //cosPhi

                let sin_theta = ((1.0 - l.powi(2)) * (1.0 - v.powi(2))).sqrt();
                let den = l.max(v);

                a + b*p*sin_theta/den
            }
        }
    }

    pub fn lighten_hit(&self, color: &Color, lightv: &Vector, intensity: &Color, precomputation: &Precomputation) -> Color {
        let n = precomputation.normal;
        let cos_theta = lightv.dot(&n);
        if cos_theta >= 0.0 {
            // The energy of the light hitting the surface depends on the cosine of the angle
		    // between the light direction and the surface normal (Lambert's cosine law) i.e. cosTheta
            let mut result = *color * self.diffuse * cos_theta * Self::oren_nayar(&precomputation.eye, &lightv, &n, self.roughness);

            // The Blinn-Phong model accounts for light that may be reflected directly towards the eye,
            // controlled by Specular (intensity of reflected light) and Shininess
            // (size of reflecting area, higher values yield a smaller area with harder reflection)
            let halfv = (*lightv + precomputation.eye).normalize();
            let n_dot_h = n.dot(&halfv);
            if n_dot_h > EPSILON {
                let f = n_dot_h.powf(4.0 * self.shininess); // Multiply shininess by 4 to keep "compatibility" with values tuned for the standard Phong model
                result = result + *intensity * self.specular * f; // add specular component
            }
            result * *intensity
        } else {
            Color::black()
        }
    }

    pub fn effective_color(&self, object: &Box<Shape>, point: &Point) -> Color {
        if let Some(pattern) = &self.pattern {
            pattern.color_at_object(&object, &point)
        } else {
            self.color
        }
    }

    pub fn lighting(
        &self,
        object: &Box<Shape>,
        lights: &[Box<Light>],
        comps: &Precomputation,
        world: &World,
    ) -> Color {
        let eff_color = self.effective_color(&object, &comps.over_point);
        let mut c = eff_color * self.ambient;

        for light in lights {
            c = c + light.intensity_at(&eff_color, &comps, &self, &world);
        }

        c
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::patterns::StripePattern;
    use crate::shapes::Sphere;
    use crate::linear::Point;

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
