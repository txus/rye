use crate::linear::{Matrix, Matrix4, Point, Vector, EPSILON};
use crate::materials::Material;
use crate::rays::{Intersection, Ray};

pub trait Shape {
    fn transform(&self) -> &Matrix4;
    fn set_transform(&mut self, t: Matrix4);
    fn material(&self) -> &Material;
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

pub struct Sphere {
    pub transform: Matrix4,
    pub material: Material,
}

impl Sphere {
    pub fn new() -> Self {
        Sphere {
            material: Material::default(),
            transform: Matrix4::id(),
        }
    }
}

impl Shape for Sphere {
    fn transform(&self) -> &Matrix4 {
        &self.transform
    }
    fn set_transform(&mut self, t: Matrix4) {
        self.transform = t
    }
    fn material(&self) -> &Material {
        &self.material
    }
    fn set_material(&mut self, m: Material) {
        self.material = m;
    }

    fn local_normal_at(&self, p: &Point) -> Vector {
        *p - Point::origin()
    }

    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let sphere_to_ray = ray.origin - Point::origin();
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
    material: Material,
}

impl Plane {
    pub fn new() -> Plane {
        Plane {
            transform: Matrix4::id(),
            material: Material::default(),
        }
    }
}

impl Shape for Plane {
    fn set_transform(&mut self, t: Matrix4) {
        self.transform = t
    }
    fn set_material(&mut self, m: Material) {
        self.material = m
    }
    fn material(&self) -> &Material {
        &self.material
    }
    fn transform(&self) -> &Matrix4 {
        &self.transform
    }
    fn local_normal_at(&self, _p: &Point) -> Vector {
        Vector::new(0.0, 1.0, 0.0)
    }
    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
        if ray.direction.y.abs() < EPSILON {
            vec![]
        } else {
            let t = -ray.origin.y / ray.direction.y;
            let shape: &Shape = self;
            vec![Intersection { t, object: shape }]
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

    struct TestShape {
        transform: Matrix4,
        material: Material,
        expected_ray: Ray,
    }

    impl TestShape {
        pub fn new(expected_ray: Ray) -> TestShape {
            TestShape {
                transform: Matrix4::id(),
                material: Material::default(),
                expected_ray,
            }
        }
    }

    impl Shape for TestShape {
        fn set_material(&mut self, m: Material) {
            self.material = m;
        }
        fn material(&self) -> &Material {
            &self.material
        }

        fn set_transform(&mut self, t: Matrix4) {
            self.transform = t;
        }
        fn transform(&self) -> &Matrix4 {
            &self.transform
        }

        fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
            if self.expected_ray.origin != ray.origin {
                panic!("Ray origin is not like expected")
            }
            if self.expected_ray.direction != ray.direction {
                panic!("Ray direction is not like expected")
            }
            vec![]
        }

        fn local_normal_at(&self, p: &Point) -> Vector {
            Vector::new(p.x, p.y, p.z)
        }
    }

    #[test]
    fn intersecting_scaled_shape_with_ray() {
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let expected_ray = Ray::new(Point::new(0.0, 0.0, -2.5), Vector::new(0.0, 0.0, 0.5));
        let mut s = TestShape::new(expected_ray);
        s.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        s.intersect(&r);
    }

    #[test]
    fn intersecting_translated_shape_with_ray() {
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let expected_ray = Ray::new(Point::new(-5.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let mut s = TestShape::new(expected_ray);
        s.set_transform(Matrix4::translation(5.0, 0.0, 0.0));
        s.intersect(&r);
    }

    #[test]
    fn computing_normal_on_translated_shape() {
        let expected_ray = Ray::new(Point::new(-5.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let mut s = TestShape::new(expected_ray);
        s.set_transform(Matrix4::translation(0.0, 1.0, 0.0));
        assert_eq!(
            s.normal(Point::new(0.0, 1.70711, -0.70711)),
            Vector::new(0.0, 0.70711, -0.70711)
        );
    }

    #[test]
    fn computing_normal_on_transformed_shape() {
        let expected_ray = Ray::new(Point::new(-5.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let mut s = TestShape::new(expected_ray);
        s.set_transform(Matrix4::scaling(1.0, 0.5, 1.0) * Matrix4::rotation_z(PI / 5.0));
        assert_eq!(
            s.normal(Point::new(0.0, 2_f32.sqrt() / 2.0, -2_f32.sqrt() / 2.0)),
            Vector::new(0.0, 0.97014, -0.24254)
        );
    }

    mod plane {
        use super::*;

        #[test]
        fn normal_is_constant() {
            let p = Plane::new();
            assert_eq!(
                p.local_normal_at(&Point::origin()),
                Vector::new(0.0, 1.0, 0.0)
            );
            assert_eq!(
                p.local_normal_at(&Point::new(10.0, 0.0, -10.0)),
                Vector::new(0.0, 1.0, 0.0)
            );
            assert_eq!(
                p.local_normal_at(&Point::new(-5.0, 0.0, 150.0)),
                Vector::new(0.0, 1.0, 0.0)
            );
        }

        #[test]
        fn intersect_with_ray_parallel_to_it() {
            let p = Plane::new();
            let r = Ray::new(Point::new(0.0, 10.0, 0.0), Vector::new(0.0, 0.0, 1.0));
            assert_eq!(
                p.local_intersect(&r)
                    .iter()
                    .map(|x| x.t)
                    .collect::<Vec<f32>>(),
                vec!()
            );
        }
        #[test]
        fn intersect_with_coplanar_ray() {
            let p = Plane::new();
            let r = Ray::new(Point::origin(), Vector::new(0.0, 0.0, 1.0));
            assert_eq!(
                p.local_intersect(&r)
                    .iter()
                    .map(|x| x.t)
                    .collect::<Vec<f32>>(),
                vec!()
            );
        }

        #[test]
        fn intersect_with_ray_from_above() {
            let p = Plane::new();
            let r = Ray::new(Point::new(0.0, 1.0, 0.0), Vector::new(0.0, -1.0, 0.0));
            assert_eq!(
                p.local_intersect(&r)
                    .iter()
                    .map(|x| x.t)
                    .collect::<Vec<f32>>(),
                vec!(1.0)
            );
        }

        #[test]
        fn plane_intersect_with_ray_from_below() {
            let p = Plane::new();
            let r = Ray::new(Point::new(0.0, -1.0, 0.0), Vector::new(0.0, 1.0, 0.0));
            assert_eq!(
                p.local_intersect(&r)
                    .iter()
                    .map(|x| x.t)
                    .collect::<Vec<f32>>(),
                vec!(1.0)
            );
        }
    }

    mod sphere {
        use super::*;

        #[test]
        fn transform() {
            let mut s = Sphere::new();
            assert_eq!(*s.transform(), Matrix4::id());
            let t = Matrix4::translation(2.0, 3.0, 4.0);
            s.set_transform(t);
            assert_eq!(*s.transform(), t);
        }

        #[test]
        fn scaled_sphere_intersection() {
            let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
            let mut s = Sphere::new();
            s.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));

            assert_eq!(
                s.intersect(&r).iter().map(|x| x.t).collect::<Vec<f32>>(),
                vec!(3.0, 7.0)
            );

            s.set_transform(Matrix4::translation(5.0, 0.0, 0.0));
            assert_eq!(
                s.intersect(&r).iter().map(|x| x.t).collect::<Vec<f32>>(),
                vec!()
            );
        }

        #[test]
        fn normal() {
            let mut s = Sphere::new();
            assert_eq!(
                s.normal(Point::new(1.0, 0.0, 0.0)),
                Vector::new(1.0, 0.0, 0.0)
            );
            assert_eq!(
                s.normal(Point::new(0.0, 1.0, 0.0)),
                Vector::new(0.0, 1.0, 0.0)
            );
            assert_eq!(
                s.normal(Point::new(0.0, 0.0, 1.0)),
                Vector::new(0.0, 0.0, 1.0)
            );
            let v = s.normal(Point::new(
                3_f32.sqrt() / 3.0,
                3_f32.sqrt() / 3.0,
                3_f32.sqrt() / 3.0,
            ));
            assert_eq!(v, v.normalize());

            let translate = Matrix4::translation(0.0, 1.0, 0.0);
            let scale_and_rotate = Matrix4::scaling(1.0, 0.5, 1.0) * Matrix4::rotation_z(PI / 5.0);

            s.set_transform(translate);
            assert_eq!(
                s.normal(Point::new(0.0, 1.70711, -0.70711)),
                Vector::new(0.0, 0.70711, -0.70711)
            );

            s.set_transform(scale_and_rotate);
            assert_eq!(
                s.normal(Point::new(0.0, 2_f32.sqrt() / 2.0, -2_f32.sqrt() / 2.0)),
                Vector::new(0.0, 0.970140, -0.24254)
            );
        }

        #[test]
        fn intersection() {
            let mut r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
            let s = Sphere::new();
            let mut is = s.intersect(&r);
            assert_eq!(is.iter().map(|x| x.t).collect::<Vec<f32>>(), vec!(4.0, 6.0));

            r = Ray::new(Point::new(0.0, 1.0, -5.0), Vector::new(0.0, 0.0, 1.0));
            is = s.intersect(&r);
            assert_eq!(is.iter().map(|x| x.t).collect::<Vec<f32>>(), vec!(5.0, 5.0));

            r = Ray::new(Point::new(0.0, 2.0, -5.0), Vector::new(0.0, 0.0, 1.0));
            assert_eq!(
                s.intersect(&r).iter().map(|x| x.t).collect::<Vec<f32>>(),
                vec!()
            );

            r = Ray::new(Point::origin(), Vector::new(0.0, 0.0, 1.0));
            is = s.intersect(&r);
            assert_eq!(
                is.iter().map(|x| x.t).collect::<Vec<f32>>(),
                vec!(-1.0, 1.0)
            );

            r = Ray::new(Point::new(0.0, 0.0, 5.0), Vector::new(0.0, 0.0, 1.0));
            is = s.intersect(&r);
            assert_eq!(
                is.iter().map(|x| x.t).collect::<Vec<f32>>(),
                vec!(-6.0, -4.0)
            );
        }
    }
}
