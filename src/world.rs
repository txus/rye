use crate::color::Color;
use crate::light::PointLight;
use crate::linear::{Matrix4, Point};
use crate::materials::Material;
use crate::rays::{Intersection, Precomputation, Ray};
use crate::shapes::{Shape, Sphere};

pub struct World {
    pub objects: Vec<Box<Shape>>,
    pub light_source: PointLight,
}

impl World {
    pub fn default() -> Self {
        let light = PointLight::new(Point::new(-10.0, 10.0, -10.0), Color::white());
        let mut s1: Box<Shape> = Box::from(Sphere::new());
        s1.set_material(Material {
            color: Color::new(0.8, 1.0, 0.6),
            diffuse: 0.7,
            specular: 0.2,
            ..Material::default()
        });
        let mut s2: Box<Shape> = Box::from(Sphere::new());
        s2.set_transform(Matrix4::scaling(0.5, 0.5, 0.5));

        World {
            objects: vec![s1, s2],
            light_source: light,
        }
    }

    pub fn intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let mut out: Vec<Intersection> = vec![];
        for object in &self.objects {
            out.append(&mut object.intersect(&ray))
        }
        out.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
        out
    }

    pub fn shade(&self, c: &Precomputation) -> Color {
        c.object.material().lighting(
            c.object,
            &self.light_source,
            &c.over_point,
            &c.eye,
            &c.normal,
            self.is_shadowed(c.over_point),
        )
    }

    pub fn color_at(&self, r: &Ray) -> Color {
        if let Some(i) = self.intersect(&r).first() {
            self.shade(&i.precompute(&r))
        } else {
            Color::black()
        }
    }

    pub fn is_shadowed(&self, p: Point) -> bool {
        let v = self.light_source.position - p;
        let distance = v.magnitude();
        let direction = v.normalize();

        let r = Ray::new(p, direction);
        let mut intersections = self.intersect(&r);

        match Intersection::hit(&mut intersections) {
            Some(h) if h.t < distance => true,
            _ => false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::color::Color;
    use crate::light::PointLight;
    use crate::linear::{Matrix4, Point, Vector, EPSILON};
    use crate::materials::Material;

    #[test]
    fn default_world() {
        let light = PointLight::new(Point::new(-10.0, 10.0, -10.0), Color::white());
        let mut s1 = Sphere::new();
        s1.material = Material {
            color: Color::new(0.8, 1.0, 0.6),
            diffuse: 0.7,
            specular: 0.2,
            ..s1.material
        };
        let mut s2 = Sphere::new();
        s2.set_transform(Matrix4::scaling(0.5, 0.5, 0.5));

        let w = World::default();
        assert_eq!(w.light_source, light);
    }

    #[test]
    fn intersect_world() {
        let w = World::default();
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let hits = w.intersect(&r).iter().map(|x| x.t).collect::<Vec<f32>>();
        assert_eq!(hits, vec!(4.0, 4.5, 5.5, 6.0));
    }

    #[test]
    fn shading_intersection() {
        let w = World::default();
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let shape: &Box<Shape> = w.objects.first().unwrap();
        let i = Intersection {
            t: 4.0,
            object: &**shape,
        };
        let comps = i.precompute(&r);
        assert_eq!(w.shade(&comps), Color::new(0.38066, 0.47583, 0.2855));
    }

    #[test]
    fn shading_intersection_from_inside() {
        let mut w = World::default();
        w.light_source = PointLight::new(Point::new(0.0, 0.25, 0.0), Color::white());
        let r = Ray::new(Point::origin(), Vector::new(0.0, 0.0, 1.0));
        let shape: &Box<Shape> = w.objects.last().unwrap();
        let i = Intersection {
            t: 0.5,
            object: &**shape,
        };
        let comps = i.precompute(&r);
        assert_eq!(w.shade(&comps), Color::new(0.90498, 0.90498, 0.90498));
    }

    #[test]
    fn color_at_miss() {
        let w = World::default();
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 1.0, 0.0));
        let c = w.color_at(&r);
        assert_eq!(c, Color::black());
    }

    #[test]
    fn color_at_hit() {
        let w = World::default();
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let c = w.color_at(&r);
        assert_eq!(c, Color::new(0.38066, 0.47582, 0.2855));
    }

    #[test]
    fn no_shadow_when_nothing_is_collinear_with_point_and_light() {
        let w = World::default();
        let p = Point::new(0.0, 10.0, 0.0);
        assert_eq!(w.is_shadowed(p), false);
    }

    #[test]
    fn shadow_when_object_is_between_light_and_point() {
        let w = World::default();
        let p = Point::new(10.0, -10.0, 10.0);
        assert_eq!(w.is_shadowed(p), true);
    }

    #[test]
    fn no_shadow_when_object_is_behind_the_light() {
        let w = World::default();
        let p = Point::new(-20.0, 20.0, -20.0);
        assert_eq!(w.is_shadowed(p), false);
    }

    #[test]
    fn no_shadow_when_object_is_behind_the_point() {
        let w = World::default();
        let p = Point::new(-2.0, 2.0, -2.0);
        assert_eq!(w.is_shadowed(p), false);
    }

    #[test]
    fn hit_should_offset_the_point() {
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let mut s = Sphere::new();
        s.set_transform(Matrix4::translation(0.0, 0.0, 1.0));
        let i = Intersection { t: 5.0, object: &s };
        let comps = i.precompute(&r);
        assert!(comps.over_point.z < -EPSILON / 2.0);
        assert!(comps.point.z > comps.over_point.z);
    }
}
