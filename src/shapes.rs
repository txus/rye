use crate::linear::{Matrix, Matrix4, Point, Vector, EPSILON};
use crate::materials::Material;
use crate::rays::{Intersection, Ray};

use std::f32::INFINITY;

use rand::Rng;

fn gen_id() -> i32 {
    let mut rng = rand::thread_rng();
    rng.gen()
}

pub trait Shape: Send + Sync {
    fn id(&self) -> i32;
    fn casts_shadows(&self) -> bool;
    fn transform(&self) -> &Matrix4;
    fn inverse_transform(&self) -> &Matrix4;
    fn set_transform(&mut self, t: Matrix4);
    fn material(&self) -> &Material;
    fn set_material(&mut self, m: Material);
    fn normal(&self, world_point: Point) -> Vector {
        let inverse = self.inverse_transform();
        let object_point = *inverse * world_point;
        let object_normal = self.local_normal_at(&object_point);
        let world_normal = inverse.transpose() * object_normal;
        world_normal.normalize()
    }
    fn local_normal_at(&self, p: &Point) -> Vector;
    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection>;
    fn intersect(&self, r: &Ray) -> Vec<Intersection> {
        let ray = r.transform(self.inverse_transform());
        self.local_intersect(&ray)
    }
}

impl std::cmp::PartialEq for Shape {
    fn eq(&self, other: &Shape) -> bool {
        self.id() == other.id()
    }
}

pub struct Sphere {
    pub id: i32,
    pub transform: Matrix4,
    pub material: Material,
    pub casts_shadows: bool,
    pub inverse_transform: Matrix4,
}

impl Sphere {
    pub fn new() -> Self {
        Sphere {
            id: gen_id(),
            material: Material::default(),
            transform: Matrix4::id(),
            casts_shadows: true,
            inverse_transform: Matrix4::id().inverse(),
        }
    }

    pub fn glass() -> Self {
        Sphere {
            id: gen_id(),
            material: Material {
                transparency: 1.0,
                refractive_index: 1.5,
                ..Material::default()
            },
            transform: Matrix4::id(),
            casts_shadows: true,
            inverse_transform: Matrix4::id().inverse(),
        }
    }
}

impl Shape for Sphere {
    fn id(&self) -> i32 {
        self.id
    }
    fn casts_shadows(&self) -> bool {
        self.casts_shadows
    }
    fn inverse_transform(&self) -> &Matrix4 {
        &self.inverse_transform
    }
    fn transform(&self) -> &Matrix4 {
        &self.transform
    }
    fn set_transform(&mut self, t: Matrix4) {
        self.transform = t;
        self.inverse_transform = t.inverse();
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
    pub id: i32,
    pub transform: Matrix4,
    pub material: Material,
    pub casts_shadows: bool,
    inverse_transform: Matrix4,
}

impl Plane {
    pub fn new() -> Plane {
        Plane {
            id: gen_id(),
            transform: Matrix4::id(),
            material: Material::default(),
            casts_shadows: true,
            inverse_transform: Matrix4::id().inverse(),
        }
    }
}

impl Shape for Plane {
    fn id(&self) -> i32 {
        self.id
    }
    fn casts_shadows(&self) -> bool {
        self.casts_shadows
    }
    fn set_transform(&mut self, t: Matrix4) {
        self.transform = t;
        self.inverse_transform = t.inverse();
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
    fn inverse_transform(&self) -> &Matrix4 {
        &self.inverse_transform
    }
    fn local_normal_at(&self, _p: &Point) -> Vector {
        Vector::new(0.0, 1.0, 0.0)
    }
    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
        if ray.direction.y.abs() < EPSILON {
            vec![]
        } else {
            let t = -(ray.origin.y) / ray.direction.y;
            let shape: &Shape = self;
            vec![Intersection { t, object: shape }]
        }
    }
}

pub struct Cube {
    pub id: i32,
    pub transform: Matrix4,
    pub material: Material,
    pub casts_shadows: bool,
    inverse_transform: Matrix4,
}

impl Cube {
    pub fn new() -> Cube {
        Cube {
            id: gen_id(),
            transform: Matrix4::id(),
            material: Material::default(),
            casts_shadows: true,
            inverse_transform: Matrix4::id().inverse(),
        }
    }
}

fn check_axis(origin: f32, direction: f32) -> (f32, f32) {
    let tmin_numerator = -1.0 - origin;
    let tmax_numerator = 1.0 - origin;

    let (tmin, tmax) = if direction.abs() >= EPSILON {
        (tmin_numerator / direction, tmax_numerator / direction)
    } else {
        (tmin_numerator * INFINITY, tmax_numerator * INFINITY)
    };

    if tmin > tmax {
        (tmax, tmin)
    } else {
        (tmin, tmax)
    }
}

impl Shape for Cube {
    fn id(&self) -> i32 {
        self.id
    }
    fn casts_shadows(&self) -> bool {
        self.casts_shadows
    }
    fn set_transform(&mut self, t: Matrix4) {
        self.transform = t;
        self.inverse_transform = t.inverse();
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
    fn inverse_transform(&self) -> &Matrix4 {
        &self.inverse_transform
    }
    fn local_normal_at(&self, p: &Point) -> Vector {
        let (x, y, z) = (p.x.abs(), p.y.abs(), p.z.abs());
        let maxc = x.max(y).max(z);
        if maxc == x {
            Vector::new(p.x, 0.0, 0.0)
        } else if maxc == y {
            Vector::new(0.0, p.y, 0.0)
        } else {
            Vector::new(0.0, 0.0, p.z)
        }
    }
    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let (xtmin, xtmax) = check_axis(ray.origin.x, ray.direction.x);
        let (ytmin, ytmax) = check_axis(ray.origin.y, ray.direction.y);
        let (ztmin, ztmax) = check_axis(ray.origin.z, ray.direction.z);

        let tmin = xtmin.max(ytmin).max(ztmin);
        let tmax = xtmax.min(ytmax).min(ztmax);

        if tmin > tmax {
            vec![]
        } else {
            let shape: &Shape = self;
            vec![
                Intersection {
                    t: tmin,
                    object: shape,
                },
                Intersection {
                    t: tmax,
                    object: shape,
                },
            ]
        }
    }
}

pub struct Cylinder {
    pub id: i32,
    pub transform: Matrix4,
    pub material: Material,
    pub casts_shadows: bool,
    pub minimum: f32,
    pub maximum: f32,
    pub closed: bool,
    inverse_transform: Matrix4,
}

impl Cylinder {
    pub fn new() -> Cylinder {
        Cylinder {
            id: gen_id(),
            transform: Matrix4::id(),
            material: Material::default(),
            casts_shadows: true,
            minimum: -INFINITY,
            maximum: INFINITY,
            closed: false,
            inverse_transform: Matrix4::id().inverse(),
        }
    }

    pub fn open(minimum: f32, maximum: f32) -> Cylinder {
        let mut c = Cylinder::new();
        c.minimum = minimum;
        c.maximum = maximum;
        c
    }

    pub fn closed(minimum: f32, maximum: f32) -> Cylinder {
        let mut c = Cylinder::open(minimum, maximum);
        c.closed = true;
        c
    }

    fn check_cap(ray: &Ray, t: f32) -> bool {
        let x = ray.origin.x + t * ray.direction.x;
        let z = ray.origin.z + t * ray.direction.z;
        (x.powf(2.0) + z.powf(2.0)) <= 1.0
    }

    fn intersect_caps<'a>(&'a self, ray: &Ray) -> Vec<Intersection<'a>> {
        let mut is = vec![];
        if !self.closed || ray.direction.y.abs() < EPSILON {
            return is;
        }
        let shape: &Shape = self;
        let lower_t = (self.minimum - ray.origin.y) / ray.direction.y;
        if Self::check_cap(ray, lower_t) {
            is.push(Intersection {
                t: lower_t,
                object: shape,
            })
        }

        let upper_t = (self.maximum - ray.origin.y) / ray.direction.y;
        if Self::check_cap(ray, upper_t) {
            is.push(Intersection {
                t: upper_t,
                object: shape,
            })
        }
        is
    }
}

impl Shape for Cylinder {
    fn id(&self) -> i32 {
        self.id
    }
    fn casts_shadows(&self) -> bool {
        self.casts_shadows
    }
    fn set_transform(&mut self, t: Matrix4) {
        self.transform = t;
        self.inverse_transform = t.inverse();
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
    fn inverse_transform(&self) -> &Matrix4 {
        &self.inverse_transform
    }
    fn local_normal_at(&self, p: &Point) -> Vector {
        let dist = p.x.powf(2.0) + p.z.powf(2.0);
        if dist < 1.0 && p.y >= self.maximum - EPSILON {
            Vector::new(0.0, 1.0, 0.0)
        } else if dist < 1.0 && p.y <= self.minimum + EPSILON {
            Vector::new(0.0, -1.0, 0.0)
        } else {
            Vector::new(p.x, 0.0, p.z)
        }
    }
    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let mut intersections = vec![];
        let a = ray.direction.x.powf(2.0) + ray.direction.z.powf(2.0);
        if a.abs() > EPSILON {
            // if the body of the cylinder intersects
            let b = 2.0 * ray.origin.x * ray.direction.x + 2.0 * ray.origin.z * ray.direction.z;
            let c = ray.origin.x.powf(2.0) + ray.origin.z.powf(2.0) - 1.0;

            let disc = b.powf(2.0) - 4.0 * a * c;

            if disc >= 0.0 {
                let t0_ = (-b - disc.sqrt()) / (2.0 * a);
                let t1_ = (-b + disc.sqrt()) / (2.0 * a);
                let mut t0 = t0_;
                let mut t1 = t1_;
                if t0_ > t1_ {
                    t0 = t1_;
                    t1 = t0_;
                }
                let shape: &Shape = self;
                let y0 = ray.origin.y + t0 * ray.direction.y;
                if self.minimum < y0 && y0 < self.maximum {
                    intersections.push(Intersection {
                        t: t0,
                        object: shape,
                    });
                }
                let y1 = ray.origin.y + t1 * ray.direction.y;
                if self.minimum < y1 && y1 < self.maximum {
                    intersections.push(Intersection {
                        t: t1,
                        object: shape,
                    });
                }
            }
        }
        let mut cap_intersections = self.intersect_caps(ray);
        intersections.append(&mut cap_intersections);
        intersections
    }
}

pub struct Cone {
    pub id: i32,
    pub transform: Matrix4,
    pub material: Material,
    pub casts_shadows: bool,
    pub minimum: f32,
    pub maximum: f32,
    pub closed: bool,
    inverse_transform: Matrix4,
}

impl Cone {
    pub fn new() -> Cone {
        Cone {
            id: gen_id(),
            transform: Matrix4::id(),
            material: Material::default(),
            casts_shadows: true,
            minimum: -INFINITY,
            maximum: INFINITY,
            closed: false,
            inverse_transform: Matrix4::id().inverse(),
        }
    }

    pub fn open(minimum: f32, maximum: f32) -> Cone {
        let mut c = Cone::new();
        c.minimum = minimum;
        c.maximum = maximum;
        c
    }

    pub fn closed(minimum: f32, maximum: f32) -> Cone {
        let mut c = Cone::open(minimum, maximum);
        c.closed = true;
        c
    }

    fn check_cap(ray: &Ray, radius: f32, t: f32) -> bool {
        let x = ray.origin.x + t * ray.direction.x;
        let z = ray.origin.z + t * ray.direction.z;
        (x.powf(2.0) + z.powf(2.0)) <= radius
    }

    fn intersect_caps<'a>(&'a self, ray: &Ray) -> Vec<Intersection<'a>> {
        let mut is = vec![];
        if !self.closed || ray.direction.y.abs() < EPSILON {
            return is;
        }
        let shape: &Shape = self;
        let lower_t = (self.minimum - ray.origin.y) / ray.direction.y;
        if Self::check_cap(ray, lower_t, self.minimum) {
            is.push(Intersection {
                t: lower_t,
                object: shape,
            })
        }

        let upper_t = (self.maximum - ray.origin.y) / ray.direction.y;
        if Self::check_cap(ray, upper_t, self.maximum) {
            is.push(Intersection {
                t: upper_t,
                object: shape,
            })
        }
        is
    }
}

impl Shape for Cone {
    fn id(&self) -> i32 {
        self.id
    }
    fn casts_shadows(&self) -> bool {
        self.casts_shadows
    }
    fn set_transform(&mut self, t: Matrix4) {
        self.transform = t;
        self.inverse_transform = t.inverse();
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
    fn inverse_transform(&self) -> &Matrix4 {
        &self.inverse_transform
    }
    fn local_normal_at(&self, p: &Point) -> Vector {
        let dist = p.x.powf(2.0) + p.z.powf(2.0);
        if dist < 1.0 && p.y >= self.maximum - EPSILON {
            Vector::new(0.0, 1.0, 0.0)
        } else if dist < 1.0 && p.y <= self.minimum + EPSILON {
            Vector::new(0.0, -1.0, 0.0)
        } else {
            Vector::new(p.x, 0.0, p.z)
        }
    }
    fn local_intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let shape: &Shape = self;
        let mut intersections = vec![];
        let a = ray.direction.x.powf(2.0) - ray.direction.y.powf(2.0) + ray.direction.z.powf(2.0);
        let b = 2.0
            * (ray.origin.x * ray.direction.x - ray.origin.y * ray.direction.y
                + ray.origin.z * ray.direction.z);
        let c = ray.origin.x.powf(2.0) - ray.origin.y.powf(2.0) + ray.origin.z.powf(2.0);
        if a.abs() > EPSILON {
            // if the body of the cone intersects
            let disc = b.powf(2.0) - 4.0 * a * c;

            if disc >= 0.0 {
                let sq = disc.sqrt();
                let t0 = (-b - sq) / (2.0 * a);
                let t1 = (-b + sq) / (2.0 * a);

                let y0 = ray.origin.y + t0 * ray.direction.y;
                if self.minimum < y0 && y0 < self.maximum {
                    intersections.push(Intersection {
                        t: t0,
                        object: shape,
                    });
                }
                let y1 = ray.origin.y + t1 * ray.direction.y;
                if self.minimum < y1 && y1 < self.maximum {
                    intersections.push(Intersection {
                        t: t1,
                        object: shape,
                    });
                }
            }
        } else {
            if b.abs() > EPSILON {
                let t = -c / (2.0 * b);
                let y = ray.origin.y + t * ray.direction.y;
                if self.minimum < y && y < self.maximum {
                    intersections.push(Intersection {
                        t: (-c) / 2.0 * b,
                        object: shape,
                    });
                }
            }
        }
        if self.closed {
            let t0 = (self.minimum - ray.origin.y) / ray.direction.y;
            if Self::check_cap(ray, self.minimum.abs(), t0) {
                intersections.push(Intersection {
                    t: t0,
                    object: shape,
                });
            }

            let t1 = (self.maximum - ray.origin.y) / ray.direction.y;
            if Self::check_cap(ray, self.maximum.abs(), t1) {
                intersections.push(Intersection {
                    t: t1,
                    object: shape,
                });
            }
        }
        intersections
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;

    struct TestShape {
        transform: Matrix4,
        inverse_transform: Matrix4,
        material: Material,
        expected_ray: Ray,
    }

    impl TestShape {
        pub fn new(expected_ray: Ray) -> TestShape {
            TestShape {
                transform: Matrix4::id(),
                material: Material::default(),
                inverse_transform: Matrix4::id().inverse(),
                expected_ray,
            }
        }
    }

    impl Shape for TestShape {
        fn id(&self) -> i32 {
            0
        }
        fn casts_shadows(&self) -> bool {
            true
        }
        fn set_material(&mut self, m: Material) {
            self.material = m;
        }
        fn material(&self) -> &Material {
            &self.material
        }
        fn set_transform(&mut self, t: Matrix4) {
            self.transform = t;
            self.inverse_transform = t.inverse();
        }
        fn transform(&self) -> &Matrix4 {
            &self.transform
        }
        fn inverse_transform(&self) -> &Matrix4 {
            &self.inverse_transform
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

    mod cube {
        use super::*;

        #[test]
        fn intersection() {
            let c = Cube::new();
            let examples = [
                // origin, direction, t1, t2
                (
                    Point::new(5.0, 0.5, 0.0),
                    Vector::new(-1.0, 0.0, 0.0),
                    4.0,
                    6.0,
                ),
                (
                    Point::new(-5.0, 0.5, 0.0),
                    Vector::new(1.0, 0.0, 0.0),
                    4.0,
                    6.0,
                ),
                (
                    Point::new(0.5, 5.0, 0.0),
                    Vector::new(0.0, -1.0, 0.0),
                    4.0,
                    6.0,
                ),
                (
                    Point::new(0.5, -5.0, 0.0),
                    Vector::new(0.0, 1.0, 0.0),
                    4.0,
                    6.0,
                ),
                (
                    Point::new(0.5, 0.0, 5.0),
                    Vector::new(0.0, 0.0, -1.0),
                    4.0,
                    6.0,
                ),
                (
                    Point::new(0.5, 0.0, -5.0),
                    Vector::new(0.0, 0.0, 1.0),
                    4.0,
                    6.0,
                ),
                (
                    Point::new(0.0, 0.5, 0.0),
                    Vector::new(0.0, 0.0, 1.0),
                    -1.0,
                    1.0,
                ),
            ];

            for (origin, direction, t1, t2) in &examples {
                let r = Ray::new(*origin, *direction);
                let hits = c
                    .local_intersect(&r)
                    .iter()
                    .map(|x| x.t)
                    .collect::<Vec<f32>>();
                assert_eq!(hits, [*t1, *t2]);
            }
        }

        #[test]
        fn ray_misses_cube() {
            let c = Cube::new();
            let examples = [
                (
                    Point::new(-2.0, 0.0, 0.0),
                    Vector::new(0.2673, 0.5345, 0.8018),
                ),
                (
                    Point::new(0.0, -2.0, 0.0),
                    Vector::new(0.8018, 0.2673, 0.5345),
                ),
                (
                    Point::new(0.0, 0.0, -2.0),
                    Vector::new(0.5345, 0.8018, 0.2673),
                ),
                (Point::new(2.0, 0.0, 2.0), Vector::new(0.0, 0.0, -1.0)),
                (Point::new(0.0, 2.0, 2.0), Vector::new(0.0, -1.0, 0.0)),
                (Point::new(2.0, 2.0, 0.0), Vector::new(-1.0, 0.0, 0.0)),
            ];

            for (origin, direction) in &examples {
                let r = Ray::new(*origin, *direction);
                let hits = c
                    .local_intersect(&r)
                    .iter()
                    .map(|x| x.t)
                    .collect::<Vec<f32>>();
                assert_eq!(hits, []);
            }
        }

        #[test]
        fn normal_on_surface() {
            let c = Cube::new();
            let examples = [
                (Point::new(1.0, 0.5, -0.8), Vector::new(1.0, 0.0, 0.0)),
                (Point::new(-1.0, -0.2, 0.9), Vector::new(-1.0, 0.0, 0.0)),
                (Point::new(-0.4, 1.0, -0.1), Vector::new(0.0, 1.0, 0.0)),
                (Point::new(0.3, -1.0, -0.7), Vector::new(0.0, -1.0, 0.0)),
                (Point::new(-0.6, 0.3, 1.0), Vector::new(0.0, 0.0, 1.0)),
                (Point::new(0.4, 0.4, -1.0), Vector::new(0.0, 0.0, -1.0)),
                (Point::new(1.0, 1.0, 1.0), Vector::new(1.0, 0.0, 0.0)),
                (Point::new(-1.0, -1.0, -1.0), Vector::new(-1.0, 0.0, 0.0)),
            ];

            for (point, normal) in &examples {
                assert_eq!(c.local_normal_at(&point), *normal);
            }
        }
    }

    mod cylinder {
        use super::*;

        #[test]
        fn ray_strikes_cylinder() {
            let c = Cylinder::new();
            let examples = [
                // origin, direction, t1, t2
                (
                    Point::new(1.0, 0.0, -5.0),
                    Vector::new(0.0, 0.0, 1.0),
                    5.0,
                    5.0,
                ),
                (
                    Point::new(0.0, 0.0, -5.0),
                    Vector::new(0.0, 0.0, 1.0),
                    4.0,
                    6.0,
                ),
                (
                    Point::new(0.5, 0.0, -5.0),
                    Vector::new(0.1, 1.0, 1.0),
                    6.808006,
                    7.0886984,
                ),
            ];

            for (origin, direction, t1, t2) in &examples {
                let r = Ray::new(*origin, direction.normalize());
                let hits = c
                    .local_intersect(&r)
                    .iter()
                    .map(|x| x.t)
                    .collect::<Vec<f32>>();
                assert_eq!(hits, [*t1, *t2]);
            }
        }

        #[test]
        fn ray_misses_cylinder() {
            let c = Cylinder::new();
            let examples = [
                (Point::new(1.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
                (Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
                (Point::new(0.0, 0.0, -5.0), Vector::new(1.0, 1.0, 1.0)),
            ];

            for (origin, direction) in &examples {
                let r = Ray::new(*origin, direction.normalize());
                let hits = c
                    .local_intersect(&r)
                    .iter()
                    .map(|x| x.t)
                    .collect::<Vec<f32>>();
                assert_eq!(hits, []);
            }
        }

        #[test]
        fn normal_vector() {
            let c = Cylinder::new();
            let examples = [
                (Point::new(1.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)),
                (Point::new(0.0, 5.0, -1.0), Vector::new(0.0, 0.0, -1.0)),
                (Point::new(0.0, -2.0, 1.0), Vector::new(0.0, 0.0, 1.0)),
                (Point::new(-1.0, 1.0, 0.0), Vector::new(-1.0, 0.0, 0.0)),
            ];

            for (point, normal) in &examples {
                assert_eq!(c.local_normal_at(&point), *normal);
            }
        }

        #[test]
        fn default_minimum_maximum() {
            let c = Cylinder::new();
            assert_eq!(c.minimum, -INFINITY);
            assert_eq!(c.maximum, INFINITY);
        }

        #[test]
        fn intersecting_constrained_cylinder() {
            let c = Cylinder::open(1.0, 2.0);
            let examples: &[(Point, Vector, usize)] = &[
                // origin, direction, intersection count
                (Point::new(0.0, 1.5, 0.0), Vector::new(0.1, 1.0, 0.0), 0),
                (Point::new(0.0, 3.0, -5.0), Vector::new(0.0, 0.0, 1.0), 0),
                (Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0), 0),
                (Point::new(0.0, 2.0, -5.0), Vector::new(0.0, 0.0, 1.0), 0),
                (Point::new(0.0, 1.0, -5.0), Vector::new(0.0, 0.0, 1.0), 0),
                (Point::new(0.0, 1.5, -2.0), Vector::new(0.0, 0.0, 1.0), 2),
            ];

            for (origin, direction, hit_count) in examples {
                let r = Ray::new(*origin, direction.normalize());
                let hits = c
                    .local_intersect(&r)
                    .iter()
                    .map(|x| x.t)
                    .collect::<Vec<f32>>();
                assert_eq!(hits.len(), *hit_count);
            }
        }

        #[test]
        fn intersecting_caps_of_closed_cylinder() {
            let c = Cylinder::closed(1.0, 2.0);
            let examples: &[(Point, Vector, usize)] = &[
                // point, direction, intersection count
                (Point::new(0.0, 3.0, 0.0), Vector::new(0.0, -1.0, 0.0), 2),
                (Point::new(0.0, 3.0, -2.0), Vector::new(0.0, -1.0, 2.0), 2),
                (Point::new(0.0, 4.01, -2.0), Vector::new(0.0, -1.0, 1.0), 2), // edge case
                (Point::new(0.0, 0.0, -2.0), Vector::new(0.0, 1.0, 2.0), 2),
                (Point::new(0.0, -1.01, -2.0), Vector::new(0.0, 1.0, 1.0), 2), // edge case
            ];

            for (origin, direction, hit_count) in examples {
                let r = Ray::new(*origin, direction.normalize());
                let hits = c
                    .local_intersect(&r)
                    .iter()
                    .map(|x| x.t)
                    .collect::<Vec<f32>>();
                assert_eq!(hits.len(), *hit_count);
            }
        }

        #[test]
        fn normal_vector_on_caps() {
            let c = Cylinder::closed(1.0, 2.0);
            let examples = [
                // point, normal
                (Point::new(0.0, 1.0, 0.0), Vector::new(0.0, -1.0, 0.0)),
                (Point::new(0.5, 1.0, 0.0), Vector::new(0.0, -1.0, 0.0)),
                (Point::new(0.0, 1.0, 0.5), Vector::new(0.0, -1.0, 0.0)),
                (Point::new(0.0, 2.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
                (Point::new(0.5, 2.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
                (Point::new(0.0, 2.0, 0.5), Vector::new(0.0, 1.0, 0.0)),
            ];

            for (point, direction) in &examples {
                assert_eq!(c.local_normal_at(&point), *direction);
            }
        }
    }

    mod cone {
        use super::*;

        #[test]
        fn ray_strikes_cone() {
            let c = Cone::new();
            let examples = [
                // origin, direction, t1, t2
                (
                    Point::new(0.0, 0.0, -5.0),
                    Vector::new(0.0, 0.0, 1.0),
                    5.0,
                    5.0,
                ),
                (
                    Point::new(0.0, 0.0, -5.0),
                    Vector::new(1.0, 1.0, 1.0),
                    8.66025,
                    8.66025,
                ),
                (
                    Point::new(1.0, 1.0, -5.0),
                    Vector::new(-0.5, -1.0, 1.0),
                    4.55006,
                    49.44994,
                ),
            ];

            for (origin, direction, t1, t2) in &examples {
                let r = Ray::new(*origin, direction.normalize());
                let hits = c
                    .local_intersect(&r)
                    .iter()
                    .map(|x| x.t)
                    .collect::<Vec<f32>>();
                assert_eq!(hits, [*t1, *t2]);
            }
        }

        #[test]
        fn ray_strikes_cone_parallel_to_one_of_its_halves() {
            let c = Cone::new();
            let origin = Point::new(0.0, 0.0, -1.0);
            let direction = Vector::new(0.0, 1.0, 1.0).normalize();
            let r = Ray::new(origin, direction);
            let hits = c
                .local_intersect(&r)
                .iter()
                .map(|x| x.t)
                .collect::<Vec<f32>>();
            assert_eq!(hits, [0.35355]);
        }

        #[test]
        fn intersecting_cone_end_caps() {
            let c = Cone::closed(-0.5, 0.5);
            let examples = [
                // origin, direction, hit_count
                (Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 1.0, 0.0), 0),
                (Point::new(0.0, 0.0, -0.25), Vector::new(0.0, 1.0, 1.0), 2),
                (Point::new(0.0, 0.0, -0.25), Vector::new(0.0, 1.0, 0.0), 4),
            ];

            for (origin, direction, hit_count) in &examples {
                let r = Ray::new(*origin, direction.normalize());
                let hits = c
                    .local_intersect(&r)
                    .iter()
                    .map(|x| x.t)
                    .collect::<Vec<f32>>();
                assert_eq!(hits.len(), *hit_count);
            }
        }
    }
}
