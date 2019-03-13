use crate::materials::Material;
use crate::primitives::{EPSILON, Matrix, Matrix4, Point, Vector};
use crate::rays::{Intersection, Ray};

pub trait Shape {
    fn transform(&self) -> Matrix4;
    fn set_transform(&mut self, t: Matrix4);
    fn material(&self) -> Material;
    fn set_material(&mut self, m: Material);
    fn normal(&self, world_point: Point) -> Vector {
        let inverse = self.transform().inverse();
        let object_point = inverse * world_point;
        let object_normal = self.local_normal_at(&object_point);
        let world_normal = inverse.transpose() * object_normal;
        world_normal.normalize()
    }
    fn local_normal_at(&self, p: &Point) -> Vector;
    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection>;
    fn intersect(&self, r: &Ray) -> Vec<Intersection> {
        let ray = r.transform(self.transform().inverse());
        self.local_intersect(&ray)
    }
}

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct Sphere {
    origin: Point,
    radius: f32,
    pub transform: Matrix4,
    pub material: Material,
}

impl Sphere {
    pub fn unit() -> Self {
        Sphere::new(1.0, Point::new(0.0, 0.0, 0.0))
    }

    pub fn new(radius: f32, origin: Point) -> Self {
        let material = Material::default();
        Sphere {
            radius,
            origin,
            material,
            transform: Matrix4::id(),
        }
    }
}

impl Shape for Sphere {
    fn transform(&self) -> Matrix4 {
        self.transform
    }
    fn set_transform(&mut self, t: Matrix4) {
        self.transform = t
    }
    fn material(&self) -> Material {
        self.material
    }
    fn set_material(&mut self, m: Material) {
        self.material = m;
    }

    fn local_normal_at(&self, p: &Point) -> Vector {
        *p - self.origin
    }

    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let sphere_to_ray = ray.origin - self.origin;
        let a = ray.direction.dot(&ray.direction);
        let b = ray.direction.dot(&sphere_to_ray) * 2.0;
        let c = sphere_to_ray.dot(&sphere_to_ray) - 1.0;
        let discriminant = b.powf(2.0) - 4.0 * a * c;

        let s: &Shape = self;

        if discriminant < 0.0 {
            vec![]
        } else {
            let t1 = (-b - discriminant.sqrt()) / (a * 2.0);
            let t2 = (-b + discriminant.sqrt()) / (a * 2.0);
            if t1 > t2 {
                vec![
                    Intersection { t: t2, object: s },
                    Intersection { t: t1, object: s },
                ]
            } else {
                vec![
                    Intersection { t: t1, object: s },
                    Intersection { t: t2, object: s },
                ]
            }
        }
    }
}

pub struct Plane {
    transform: Matrix4,
    material: Material
}

impl Plane {
    pub fn new() -> Plane {
        Plane { transform: Matrix4::id(), material: Material::default() }
    }
}

impl Shape for Plane {
    fn set_transform(&mut self, t: Matrix4) {
        self.transform = t
    }
    fn set_material(&mut self, m: Material) {
        self.material = m
    }
    fn material(&self) -> Material {
        self.material
    }
    fn transform(&self) -> Matrix4 {
        self.transform
    }
    fn local_normal_at(&self, _p: &Point) -> Vector {
        Vector::new(0.0, 1.0, 0.0)
    }
    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
        if ray.direction.y.abs() < EPSILON {
            vec!()
        } else {
            let t = -ray.origin.y / ray.direction.y;
            let shape: &Shape = self;
            vec!(Intersection { t, object: shape })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plane_normal_is_constant() {
        let p = Plane::new();
        assert_eq!(p.local_normal_at(&Point::new(0.0, 0.0, 0.0)), Vector::new(0.0, 1.0, 0.0));
        assert_eq!(p.local_normal_at(&Point::new(10.0, 0.0, -10.0)), Vector::new(0.0, 1.0, 0.0));
        assert_eq!(p.local_normal_at(&Point::new(-5.0, 0.0, 150.0)), Vector::new(0.0, 1.0, 0.0));
    }

    #[test]
    fn plane_intersect_with_ray_parallel_to_it() {
        let p = Plane::new();
        let r = Ray::new(Point::new(0.0, 10.0, 0.0), Vector::new(0.0, 0.0, 1.0));
        assert_eq!(p.local_intersect(&r).iter().map(|x| x.t).collect::<Vec<f32>>(), vec!());
    }
    #[test]
    fn plane_intersect_with_complanar_ray() {
        let p = Plane::new();
        let r = Ray::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0));
        assert_eq!(p.local_intersect(&r).iter().map(|x| x.t).collect::<Vec<f32>>(), vec!());
    }

    #[test]
    fn plane_intersect_with_ray_from_above() {
        let p = Plane::new();
        let r = Ray::new(Point::new(0.0, 1.0, 0.0), Vector::new(0.0, -1.0, 0.0));
        assert_eq!(
            p.local_intersect(&r).iter().map(|x| x.t).collect::<Vec<f32>>(),
            vec!(1.0)
            );
    }

    #[test]
    fn plane_intersect_with_ray_from_below() {
        let p = Plane::new();
        let r = Ray::new(Point::new(0.0, -1.0, 0.0), Vector::new(0.0, 1.0, 0.0));
        assert_eq!(
            p.local_intersect(&r).iter().map(|x| x.t).collect::<Vec<f32>>(),
            vec!(1.0)
            );
    }
}
