use crate::color::Color;
use crate::light::{Light, PointLight};
use crate::linear::{Matrix4, Point};
use crate::materials::Material;
use crate::rays::{Intersection, Precomputation, Ray};
use crate::shapes::{Shape, Sphere};
use crate::registry::Registry;

use std::cell::RefCell;
use std::rc::Rc;

pub const MAX_REFLECTIONS: usize = 4;

pub struct World {
    pub registry: Rc<RefCell<Registry>>,
    pub lights: Vec<Box<Light>>,
}

impl World {
    pub fn default_sphere1() -> Sphere {
        let mut s1 = Sphere::new();
        s1.set_material(Material {
            color: Color::new(0.8, 1.0, 0.6),
            diffuse: 0.7,
            specular: 0.2,
            ..Material::default()
        });
        s1
    }

    pub fn default_sphere2() -> Sphere {
        let mut s2 = Sphere::new();
        s2.set_transform(Matrix4::scaling(0.5, 0.5, 0.5));
        s2
    }

    pub fn default() -> Self {
        let light = PointLight::new(Point::new(-10.0, 10.0, -10.0), Color::white());
        let registry = Rc::new(RefCell::new(Registry::new()));
        {
            let mut reg = registry.borrow_mut();
            let s1: Box<Shape> = Box::from(Self::default_sphere1());
            let s2: Box<Shape> = Box::from(Self::default_sphere2());
            let id1 = reg.register(s1);
            let id2 = reg.register(s2);
            println!("first id: {:?}", id1);
            println!("second id: {:?}", id2);
        }

        World {
            registry: registry,
            lights: vec![Box::from(light)],
        }
    }

    pub fn empty() -> Self {
        let light = PointLight::new(Point::new(-10.0, 10.0, -10.0), Color::white());
        World {
            registry: Rc::new(RefCell::new(Registry::new())),
            lights: vec![Box::from(light)],
        }
    }

    pub fn intersect<'a>(&'a self, ray: &Ray) -> Vec<Intersection> {
        let mut out: Vec<Intersection> = vec![];
        {
            let reg = self.registry.borrow();
            let ids = reg.all_ids();

            for id in ids {
                let object = reg.get(id.clone());
                for i in object.intersect(&reg, &ray) {
                    out.push(i);
                }
            }
        }
        out.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
        out
    }

    pub fn shade(&self, c: &Precomputation, remaining: usize) -> Color {
        let reg = self.registry.borrow();
        let object = reg.get(c.object);
        let material = object.material(&reg);

        let surface = material.lighting(
            &object,
            &self.lights,
            &c,
            &self,
        );
        let reflected = self.reflected_color(&c, remaining);
        let refracted = self.refracted_color(&c, remaining);

        if material.reflective > 0.0 && material.transparency > 0.0 {
            let reflectance = c.schlick();
            surface + reflected * reflectance + refracted * (1.0 - reflectance)
        } else {
            surface + reflected + refracted
        }
    }

    pub fn color_at(&self, r: &Ray, remaining: usize) -> Color {
        let intersections = self.intersect(&r);
        let maybe_hit = Intersection::hit(&intersections);
        let reg = self.registry.borrow();
        if let Some(hit) = maybe_hit {
            self.shade(&hit.precompute(&reg, &r, &intersections), remaining)
        } else {
            Color::black()
        }
    }

    pub fn reflected_color(&self, comps: &Precomputation, remaining: usize) -> Color {
        let reg = self.registry.borrow();
        let object = reg.get(comps.object);
        let mat = object.material(&reg);
        if remaining <= 0 || mat.reflective == 0.0 {
            Color::black()
        } else {
            let reflect_ray = Ray::new(comps.over_point, comps.reflect);
            let color = self.color_at(&reflect_ray, remaining - 1);
            color * mat.reflective
        }
    }

    pub fn refracted_color(&self, comps: &Precomputation, remaining: usize) -> Color {
        let reg = self.registry.borrow();
        let object = reg.get(comps.object);
        let mat = object.material(&reg);
        if remaining == 0 || mat.transparency == 0.0 {
            Color::black()
        } else {
            let n_ratio = comps.n1 / comps.n2;
            let cos_i = comps.eye.dot(&comps.normal);
            let sin2_t = n_ratio.powi(2) * (1.0 - cos_i.powi(2));
            if sin2_t > 1.0 {
                Color::black()
            } else {
                let cos_t = (1.0 - sin2_t).sqrt();
                let direction =
                    (comps.normal * ((n_ratio * cos_i) - cos_t)) - (comps.eye * n_ratio);
                let refract_ray = Ray::new(comps.under_point, direction);
                self.color_at(&refract_ray, remaining - 1) * mat.transparency
            }
        }
    }

    pub fn is_shadowed(&self, light_position: &Point, p: &Point) -> bool {
        let v = *light_position - *p;
        let distance = v.magnitude();
        let direction = v.normalize();

        let r = Ray::new(*p, direction);
        let intersections = self.intersect(&r);
        let maybe_hit = Intersection::hit(&intersections);

        match maybe_hit {
            Some(h) if h.t < distance => {
                let reg = self.registry.borrow();
                let object = reg.get(h.object);
                object.casts_shadows()
            },
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
    use crate::shapes::Plane;
    use crate::registry::{id, id_from};
    use crate::patterns::TestPattern;
    use std::cell::RefCell;
    use std::rc::Rc;

    #[test]
    fn intersect_world() {
        let w = World::default();
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let hits = w.intersect(&r).iter().map(|x| x.t).collect::<Vec<f32>>();
        assert_eq!(hits, vec!(4.0, 4.5, 5.5, 6.0));
    }

    #[test]
    fn precomputing_state_of_intersection() {
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let registry = Rc::new(RefCell::new(Registry::new()));
        let id;
        {
            let mut reg = registry.borrow_mut();
            id = reg.register(Box::from(Sphere::new()));
        }
        let reg = registry.borrow();
        let i = Intersection { uv: None, t: 4.0, object: id };
        let is = [i];
        let comps = i.precompute(&reg, &r, &is);
        assert_eq!(comps.point, Point::new(0.0, 0.0, -1.0));
        assert_eq!(comps.eye, Vector::new(0.0, 0.0, -1.0));
        assert_eq!(comps.normal, Vector::new(0.0, 0.0, -1.0));
    }

    #[test]
    fn shading_intersection() {
        let w = World::default();
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let reg = w.registry.borrow();
        let i = Intersection { uv: None,
            t: 4.0,
            object: id(),
        };
        let is = [i];
        let comps = i.precompute(&reg, &r, &is);
        assert_eq!(w.shade(&comps, 1), Color::new(0.38066, 0.47583, 0.2855));
    }

    #[test]
    fn shading_intersection_from_inside() {
        let mut w = World::default();
        w.light_source = PointLight::new(Point::new(0.0, 0.25, 0.0), Color::white());
        let r = Ray::new(Point::origin(), Vector::new(0.0, 0.0, 1.0));
        let i = Intersection { uv: None,
            t: 0.5,
            object: id_from(2),
        };
        let is = [i];
        let reg = w.registry.borrow();
        let comps = i.precompute(&reg, &r, &is);
        assert_eq!(w.shade(&comps, 1), Color::new(0.90498, 0.90498, 0.90498));
    }

    #[test]
    fn color_at_miss() {
        let w = World::default();
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 1.0, 0.0));
        let c = w.color_at(&r, 1);
        assert_eq!(c, Color::black());
    }

    #[test]
    fn color_at_hit() {
        let w = World::default();
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let c = w.color_at(&r, 1);
        assert_eq!(c, Color::new(0.38066, 0.47582, 0.2855));
    }

    #[test]
    fn hit_should_offset_the_point() {
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let registry = Rc::new(RefCell::new(Registry::new()));
        let id;
        {
            let mut reg = registry.borrow_mut();
            id = reg.register(Box::from(Sphere::new()));
            let s = reg.get_mut(id);
            s.set_transform(Matrix4::translation(0.0, 0.0, 1.0));
        }

        let reg = registry.borrow();
        let i = Intersection { uv: None, t: 5.0, object: id };
        let mut is = [i];
        let comps = i.precompute(&reg, &r, &mut is);
        assert!(comps.over_point.z < -EPSILON / 2.0);
        assert!(comps.point.z > comps.over_point.z);
    }

    #[test]
    fn reflected_color_for_nonreflective_material() {
        let w = World::default();
        let id = id_from(2);
        {
            let mut reg = w.registry.borrow_mut();
            let s = reg.get_mut(id);
            s.set_material(Material {
                ambient: 1.0,
                ..Material::default()
            });
        }

        let reg = w.registry.borrow();
        let i = Intersection { uv: None,
            t: 1.0,
            object: id,
        };
        let is = [i];
        let r = Ray::new(Point::origin(), Vector::new(0.0, 0.0, 1.0));
        let comps = i.precompute(&reg, &r, &is);
        let color = w.reflected_color(&comps, 1);
        assert_eq!(color, Color::black());
    }

    #[test]
    fn reflected_color_for_reflective_material() {
        let w = World::default();

        let plane_id = {
            let mut reg = w.registry.borrow_mut();
            let mut p = Plane::new();
            p.set_material(Material {
                reflective: 0.5,
                ..Material::default()
            });
            p.set_transform(Matrix4::translation(0.0, -1.0, 0.0));
            reg.register(Box::from(p))
        };

        let reg = w.registry.borrow();
        let r = Ray::new(
            Point::new(0.0, 0.0, -3.0),
            Vector::new(0.0, -2_f32.sqrt() / 2.0, 2_f32.sqrt() / 2.0),
        );
        let i = Intersection { uv: None,
            t: 2_f32.sqrt(),
            object: plane_id,
        };
        let is = [i];
        let comps = i.precompute(&reg, &r, &is);
        let color = w.reflected_color(&comps, 1);
        assert_eq!(color, Color::new(0.19032, 0.2379, 0.14274));
    }

    #[test]
    fn shade_with_a_reflective_material() {
        let w = World::default();
        let id = {
            let mut reg = w.registry.borrow_mut();
            let mut p = Plane::new();
            p.set_material(Material {
                reflective: 0.5,
                ..Material::default()
            });
            p.set_transform(Matrix4::translation(0.0, -1.0, 0.0));
            reg.register(Box::from(p))
        };

        let reg = w.registry.borrow();
        let r = Ray::new(
            Point::new(0.0, 0.0, -3.0),
            Vector::new(0.0, -2_f32.sqrt() / 2.0, 2_f32.sqrt() / 2.0),
        );
        let i = Intersection { uv: None,
            t: 2_f32.sqrt(),
            object: id,
        };
        let is = [i];
        let comps = i.precompute(&reg, &r, &is);
        let color = w.shade(&comps, 1);
        assert_eq!(color, Color::new(0.87677, 0.92436, 0.82918));
    }

    #[test]
    fn color_at_with_mutually_recursive_surfaces() {
        let w = Rc::new(RefCell::new(World::default()));
        {
            let mut world = w.borrow_mut();
            world.light_source = PointLight::new(Point::origin(), Color::white());
        }
        let world = w.borrow();
        {
            let mut reg = world.registry.borrow_mut();
            let mut lower = Plane::new();
            lower.set_material(Material {
                reflective: 1.0,
                ..Material::default()
            });
            lower.set_transform(Matrix4::translation(0.0, -1.0, 0.0));
            let mut upper = Plane::new();
            upper.set_material(Material {
                reflective: 1.0,
                ..Material::default()
            });
            upper.set_transform(Matrix4::translation(0.0, 1.0, 0.0));
            reg.register(Box::from(lower));
            reg.register(Box::from(upper));
        }
        let r = Ray::new(Point::origin(), Vector::new(0.0, 1.0, 0.0));
        let color = world.color_at(&r, MAX_REFLECTIONS);
        assert_eq!(color, Color::new(1.9, 1.9, 1.9));
    }

    #[test]
    fn refracted_color_with_opaque_surface() {
        let world = World::default();
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let is = vec![
            Intersection { uv: None,
                t: 4.0,
                object: id()
            },
            Intersection { uv: None,
                t: 6.0,
                object: id()
            },
        ];
        let i = is[0];
        let comps = i.precompute(&world.registry.borrow(), &r, &is);
        let color = world.refracted_color(&comps, MAX_REFLECTIONS);
        assert_eq!(color, Color::black());
    }

    #[test]
    fn refracted_color_at_maximum_recursive_depth() {
        let w = World::default();
        let id = {
            let mut reg = w.registry.borrow_mut();
            let id = id();
            let first = reg.get_mut(id);
            first.set_material(Material {
                transparency: 1.0,
                refractive_index: 1.5,
                ..Material::default()
            });
            id
        };

        let reg = w.registry.borrow();
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let is = vec![
            Intersection { uv: None,
                t: 4.0,
                object: id,
            },
            Intersection { uv: None,
                t: 6.0,
                object: id,
            },
        ];
        let i = is[0];
        let comps = i.precompute(&reg, &r, &is);
        let color = w.refracted_color(&comps, 0);
        assert_eq!(color, Color::black());
    }

    #[test]
    fn refracted_color_under_total_internal_reflection() {
        let w = World::default();
        let id = {
            let mut reg = w.registry.borrow_mut();
            let id = id();
            let first = reg.get_mut(id);
            first.set_material(Material {
                transparency: 1.0,
                refractive_index: 1.5,
                ..Material::default()
            });
            id
        };
        let reg = w.registry.borrow();

        let r = Ray::new(
            Point::new(0.0, 0.0, 2_f32.sqrt()),
            Vector::new(0.0, 1.0, 0.0),
        );
        let is = vec![
            Intersection { uv: None,
                t: -2_f32.sqrt(),
                object: id,
            },
            Intersection { uv: None,
                t: 2_f32.sqrt(),
                object: id,
            },
        ];
        let i = is[1];
        let comps = i.precompute(&reg, &r, &is);
        let color = w.refracted_color(&comps, MAX_REFLECTIONS);
        assert_eq!(color, Color::black());
    }

    #[test]
    fn refracted_color_with_refracted_ray() {
        let w = World::default();
        let id1 = {
            let mut reg = w.registry.borrow_mut();
            let first = reg.get_mut(id());
            first.set_material(Material {
                ambient: 1.0,
                pattern: Some(Box::from(TestPattern::new())),
                ..Material::default()
            });
            id()
        };
        let id2 = {
            let mut reg = w.registry.borrow_mut();
            let second = reg.get_mut(id_from(2));
            second.set_material(Material {
                transparency: 1.0,
                refractive_index: 1.5,
                ..Material::default()
            });
            id_from(2)
        };
        let reg = w.registry.borrow();

        let r = Ray::new(Point::new(0.0, 0.0, 0.1), Vector::new(0.0, 1.0, 0.0));
        let is = vec![
            Intersection { uv: None,
                t: -0.9899,
                object: id1,
            },
            Intersection { uv: None,
                t: -0.4899,
                object: id2,
            },
            Intersection { uv: None,
                t: 0.4899,
                object: id2,
            },
            Intersection { uv: None,
                t: 0.9899,
                object: id1,
            },
        ];
        let i = is[2];
        let comps = i.precompute(&reg, &r, &is);
        ();
        let color = w.refracted_color(&comps, MAX_REFLECTIONS);
        assert_eq!(color, Color::new(0.0, 0.99888, 0.04725));
    }

    #[test]
    fn shade_with_transparent_material() {
        let w = World::default();
        let (fid, _bid) = {
            let mut reg = w.registry.borrow_mut();
            let mut floor = Plane::new();
            floor.set_transform(Matrix4::translation(0.0, -1.0, 0.0));
            floor.set_material(Material {
                transparency: 0.5,
                refractive_index: 1.5,
                ..Material::default()
            });

            let mut ball = Sphere::new();
            ball.set_transform(Matrix4::translation(0.0, -3.5, -0.5));
            ball.set_material(Material {
                color: Color::new(1.0, 0.0, 0.0),
                ambient: 0.5,
                ..Material::default()
            });
            (reg.register(Box::from(floor)), reg.register(Box::from(ball)))
        };
        let reg = w.registry.borrow();

        let r = Ray::new(
            Point::new(0.0, 0.0, -3.0),
            Vector::new(0.0, -2_f32.sqrt() / 2.0, 2_f32.sqrt() / 2.0),
        );
        let is = vec![Intersection { uv: None,
            t: 2_f32.sqrt(),
            object: fid,
        }];
        let i = is[0];
        let comps = i.precompute(&reg, &r, &is);
        let color = w.shade(&comps, MAX_REFLECTIONS);
        assert_eq!(color, Color::new(0.93642, 0.68642, 0.68642));
    }

    #[test]
    fn shade_with_reflective_transparent_material() {
        let w = World::default();
        let (fid, _bid) = {
            let mut reg = w.registry.borrow_mut();
            let mut floor = Plane::new();
            floor.set_transform(Matrix4::translation(0.0, -1.0, 0.0));
            floor.set_material(Material {
                reflective: 0.5,
                transparency: 0.5,
                refractive_index: 1.5,
                ..Material::default()
            });

            let mut ball = Sphere::new();
            ball.set_transform(Matrix4::translation(0.0, -3.5, -0.5));
            ball.set_material(Material {
                color: Color::new(1.0, 0.0, 0.0),
                ambient: 0.5,
                ..Material::default()
            });
            (reg.register(Box::from(floor)), reg.register(Box::from(ball)))
        };
        let reg = w.registry.borrow();

        let r = Ray::new(
            Point::new(0.0, 0.0, -3.0),
            Vector::new(0.0, -2_f32.sqrt() / 2.0, 2_f32.sqrt() / 2.0),
        );
        let is = vec![Intersection { uv: None,
            t: 2_f32.sqrt(),
            object: fid,
        }];
        let i = is[0];
        let comps = i.precompute(&reg, &r, &is);
        let color = w.shade(&comps, MAX_REFLECTIONS);
        assert_eq!(color, Color::new(0.93391, 0.69643, 0.69243));
    }

    #[test]
    fn is_shadowed_test_for_occlusion_between_two_points() {
        let w = World::default();
        let light_pos = Point::new(-10.0, -10.0, -10.0);

        let examples: Vec<(Point, bool)> = vec![
            (Point::new(-10.0, -10.0, 10.0), false),
            (Point::new(10.0, 10.0, 10.0), true),
            (Point::new(-20.0, -20.0, -20.0), false),
            (Point::new(-5.0, -5.0, -5.0), false),
        ];

        for (p, result) in &examples {
            assert_eq!(w.is_shadowed(&light_pos, &p), *result);
        }
    }
}
