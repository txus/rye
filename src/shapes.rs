use indextree::NodeId;
use crate::linear::{Matrix, Matrix4, Point, Vector, EPSILON};
use crate::materials::Material;
use crate::rays::{Intersection, Ray};
use crate::registry::Registry;

use std::f32::INFINITY;

use rand::Rng;

fn gen_id() -> i32 {
    let mut rng = rand::thread_rng();
    rng.gen()
}

pub trait Shape: Send + Sync {
    fn id(&self) -> i32;
    fn set_tag(&mut self, id: NodeId);
    fn tag(&self) -> NodeId;
    fn casts_shadows(&self) -> bool;
    fn transform(&self) -> &Matrix4;
    fn inverse_transform(&self) -> &Matrix4;
    fn set_transform(&mut self, t: Matrix4);
    fn material(&self) -> &Material;
    fn set_material(&mut self, m: Material);
    fn normal(&self, reg: &Registry, world_point: Point) -> Vector {
        println!("localpoint?");
        let local_point = self.world_to_object(&reg, &world_point);
        println!("localnormal?");
        let local_normal = self.local_normal_at(&local_point);
        println!("normaltoworld?");
        self.normal_to_world(&reg, &local_normal)
    }
    fn local_normal_at(&self, p: &Point) -> Vector;
    fn local_intersect<'a>(&'a self, reg: &'a Registry, ray: &Ray) -> Vec<Intersection<'a>>;
    fn intersect<'a>(&'a self, reg: &'a Registry, r: &Ray) -> Vec<Intersection<'a>> {
        let ray = r.transform(self.inverse_transform());
        self.local_intersect(&reg, &ray)
    }
    fn world_to_object(&self, reg: &Registry, p: &Point) -> Point {
        let point = if let Some(parent) = reg.parent(self.tag()) {
            parent.world_to_object(&reg, p)
        } else {
            p.clone()
        };
        *self.inverse_transform() * point
    }
    fn normal_to_world(&self, reg: &Registry, normal: &Vector) -> Vector {
        let n = (self.inverse_transform().transpose() * *normal).normalize();

        let norm = if let Some(parent) = reg.parent(self.tag()) {
            parent.normal_to_world(&reg, &n)
        } else {
            n
        };

        norm
    }
}

impl std::cmp::PartialEq for Shape {
    fn eq(&self, other: &Shape) -> bool {
        self.id() == other.id()
    }
}

pub struct Sphere {
    pub id: i32,
    pub tag: Option<NodeId>,
    pub transform: Matrix4,
    pub material: Material,
    pub casts_shadows: bool,
    pub inverse_transform: Matrix4,
}

impl Sphere {
    pub fn new() -> Self {
        Sphere {
            id: gen_id(),
            tag: None,
            material: Material::default(),
            transform: Matrix4::id(),
            casts_shadows: true,
            inverse_transform: Matrix4::id().inverse(),
        }
    }

    pub fn glass() -> Self {
        Sphere {
            id: gen_id(),
            tag: None,
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
    fn set_tag(&mut self, id: NodeId) {
        self.tag = Some(id);
    }
    fn tag(&self) -> NodeId {
        self.tag.unwrap_or_else(|| panic!("Object without tag"))
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

    fn local_intersect<'a>(&'a self, reg: &'a Registry, ray: &Ray) -> Vec<Intersection<'a>> {
        let sphere_to_ray = ray.origin - Point::origin();
        let a = ray.direction.dot(&ray.direction);
        let b = ray.direction.dot(&sphere_to_ray) * 2.0;
        let c = sphere_to_ray.dot(&sphere_to_ray) - 1.0;
        let discriminant = b.powi(2) - 4.0 * a * c;

        let s: &Box<Shape> = reg.get(self.tag()).unwrap();

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
    pub tag: Option<NodeId>,
    pub transform: Matrix4,
    pub material: Material,
    pub casts_shadows: bool,
    inverse_transform: Matrix4,
}

impl Plane {
    pub fn new() -> Plane {
        Plane {
            id: gen_id(),
            tag: None,
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
    fn set_tag(&mut self, id: NodeId) {
        self.tag = Some(id);
    }
    fn tag(&self) -> NodeId {
        self.tag.unwrap_or_else(|| panic!("Object without tag"))
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
    fn local_intersect<'a>(&'a self, reg: &'a Registry, ray: &Ray) -> Vec<Intersection<'a>> {
        if ray.direction.y.abs() < EPSILON {
            vec![]
        } else {
            let t = -(ray.origin.y) / ray.direction.y;
            let shape: &Box<Shape> = reg.get(self.tag()).unwrap();
            vec![Intersection { t, object: shape }]
        }
    }
}

pub struct Cube {
    pub id: i32,
    pub tag: Option<NodeId>,
    pub transform: Matrix4,
    pub material: Material,
    pub casts_shadows: bool,
    inverse_transform: Matrix4,
}

impl Cube {
    pub fn new() -> Cube {
        Cube {
            id: gen_id(),
            tag: None,
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
    fn set_tag(&mut self, id: NodeId) {
        self.tag = Some(id);
    }
    fn tag(&self) -> NodeId {
        self.tag.unwrap_or_else(|| panic!("Object without tag"))
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
    fn local_intersect<'a>(&'a self, reg: &'a Registry, ray: &Ray) -> Vec<Intersection<'a>> {
        let (xtmin, xtmax) = check_axis(ray.origin.x, ray.direction.x);
        let (ytmin, ytmax) = check_axis(ray.origin.y, ray.direction.y);
        let (ztmin, ztmax) = check_axis(ray.origin.z, ray.direction.z);

        let tmin = xtmin.max(ytmin).max(ztmin);
        let tmax = xtmax.min(ytmax).min(ztmax);

        if tmin > tmax {
            vec![]
        } else {
            let shape: &Box<Shape> = reg.get(self.tag()).unwrap();
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
    pub tag: Option<NodeId>,
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
            tag: None,
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
        (x.powi(2) + z.powi(2)) <= 1.0
    }

    fn intersect_caps<'a>(&'a self, reg: &'a Registry, ray: &Ray) -> Vec<Intersection<'a>> {
        let mut is = vec![];
        if !self.closed || ray.direction.y.abs() < EPSILON {
            return is;
        }
        let shape: &Box<Shape> = reg.get(self.tag()).unwrap();
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
    fn set_tag(&mut self, id: NodeId) {
        self.tag = Some(id);
    }
    fn tag(&self) -> NodeId {
        self.tag.unwrap_or_else(|| panic!("Object without tag"))
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
        let dist = p.x.powi(2) + p.z.powi(2);
        if dist < 1.0 && p.y >= self.maximum - EPSILON {
            Vector::new(0.0, 1.0, 0.0)
        } else if dist < 1.0 && p.y <= self.minimum + EPSILON {
            Vector::new(0.0, -1.0, 0.0)
        } else {
            Vector::new(p.x, 0.0, p.z)
        }
    }
    fn local_intersect<'a>(&'a self, reg: &'a Registry, ray: &Ray) -> Vec<Intersection<'a>> {
        let mut intersections = vec![];
        let a = ray.direction.x.powi(2) + ray.direction.z.powi(2);
        if a.abs() > EPSILON {
            // if the body of the cylinder intersects
            let b = 2.0 * ray.origin.x * ray.direction.x + 2.0 * ray.origin.z * ray.direction.z;
            let c = ray.origin.x.powi(2) + ray.origin.z.powi(2) - 1.0;

            let disc = b.powi(2) - 4.0 * a * c;

            if disc >= 0.0 {
                let t0_ = (-b - disc.sqrt()) / (2.0 * a);
                let t1_ = (-b + disc.sqrt()) / (2.0 * a);
                let mut t0 = t0_;
                let mut t1 = t1_;
                if t0_ > t1_ {
                    t0 = t1_;
                    t1 = t0_;
                }
                let shape: &Box<Shape> = reg.get(self.tag()).unwrap();
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
        let mut cap_intersections = self.intersect_caps(reg, ray);
        intersections.append(&mut cap_intersections);
        intersections
    }
}

pub struct Cone {
    pub id: i32,
    pub tag: Option<NodeId>,
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
            tag: None,
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
        (x.powi(2) + z.powi(2)) <= radius
    }
}

impl Shape for Cone {
    fn id(&self) -> i32 {
        self.id
    }
    fn set_tag(&mut self, id: NodeId) {
        self.tag = Some(id);
    }
    fn tag(&self) -> NodeId {
        self.tag.unwrap_or_else(|| panic!("Object without tag"))
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
        let dist = p.x.powi(2) + p.z.powi(2);
        if dist < 1.0 && p.y >= self.maximum - EPSILON {
            Vector::new(0.0, 1.0, 0.0)
        } else if dist < 1.0 && p.y <= self.minimum + EPSILON {
            Vector::new(0.0, -1.0, 0.0)
        } else {
            let mut y = (p.x.powi(2) + p.z.powi(2)).sqrt();
            if p.y > 0.0 {
                y = -y;
            }
            Vector::new(p.x, y, p.z)
        }
    }
    fn local_intersect<'a>(&'a self, reg: &'a Registry, ray: &Ray) -> Vec<Intersection<'a>> {
        let o = ray.origin;
        let d = ray.direction;
        let this: &Box<Shape> = reg.get(self.tag()).unwrap();
        let mut intersections = vec![];
        let a = d.x.powi(2) - d.y.powi(2) + d.z.powi(2);
        let b = 2.0 * o.x * d.x - 2.0 * o.y * d.y + 2.0 * o.z * d.z;
        let c = o.x.powi(2) - o.y.powi(2) + o.z.powi(2);
        if a.abs() > EPSILON {
            // if the body of the cone intersects
            let disc = b.powi(2) - 4.0 * a * c;

            if disc >= 0.0 {
                let sq = disc.sqrt();
                let t0 = (-b - sq) / (2.0 * a);
                let t1 = (-b + sq) / (2.0 * a);

                let y0 = ray.origin.y + t0 * ray.direction.y;
                if self.minimum < y0 && y0 < self.maximum {
                    intersections.push(Intersection {
                        t: t0,
                        object: this,
                    });
                }
                let y1 = ray.origin.y + t1 * ray.direction.y;
                if self.minimum < y1 && y1 < self.maximum {
                    intersections.push(Intersection {
                        t: t1,
                        object: this,
                    });
                }
            }
        } else {
            if b.abs() > EPSILON { // is the ray parallel to one of the cone halves?
                let t = -c / (2.0 * b);
                let y = ray.origin.y + t * ray.direction.y;
                if self.minimum < y && y < self.maximum {
                    intersections.push(Intersection {
                        t: t,
                        object: this,
                    });
                }
            }
        }
        // intersect with the caps
        if self.closed && ray.direction.y.abs() > EPSILON {
            let t0 = (self.minimum - ray.origin.y) / ray.direction.y;
            if Self::check_cap(ray, self.minimum.abs(), t0) {
                intersections.push(Intersection {
                    t: t0,
                    object: this,
                });
            }

            let t1 = (self.maximum - ray.origin.y) / ray.direction.y;
            if Self::check_cap(ray, self.maximum.abs(), t1) {
                intersections.push(Intersection {
                    t: t1,
                    object: this,
                });
            }
        }
        intersections
    }
}

pub struct Group {
    pub id: i32,
    pub tag: Option<NodeId>,
    pub transform: Matrix4,
    pub material: Material,
    pub casts_shadows: bool,
    inverse_transform: Matrix4,
}

impl Group {
    pub fn new() -> Group {
        let transform = Matrix4::id();
        Group {
            id: gen_id(),
            tag: None,
            transform: transform,
            inverse_transform: transform.inverse(),
            material: Material::default(),
            casts_shadows: true
        }
    }
}

impl Shape for Group {
    fn id(&self) -> i32 {
        self.id
    }
    fn set_tag(&mut self, id: NodeId) {
        self.tag = Some(id);
    }
    fn tag(&self) -> NodeId {
        self.tag.unwrap_or_else(|| panic!("Object without tag"))
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
        panic!("Calling local_normal_at in group");
    }
    fn local_intersect<'a>(&'a self, reg: &'a Registry, ray: &Ray) -> Vec<Intersection<'a>> {
        let children = reg.children(self.tag());
        let mut is = vec![];
        for child in children {
            is.append(&mut child.intersect(&reg, &ray));
        }
        is.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
        is
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;
    use std::cell::RefCell;
    use std::rc::Rc;

    struct TestShape {
        tag: Option<NodeId>,
        transform: Matrix4,
        inverse_transform: Matrix4,
        material: Material,
        expected_ray: Ray,
    }

    impl TestShape {
        pub fn new(expected_ray: Ray) -> TestShape {
            TestShape {
                tag: None,
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
        fn set_tag(&mut self, id: NodeId) {
            self.tag = Some(id);
        }
        fn tag(&self) -> NodeId {
            self.tag.unwrap()
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

        fn local_intersect(&self, reg: &Registry, ray: &Ray) -> Vec<Intersection> {
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
        let reg = Registry::new();
        s.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        s.intersect(&reg, &r);
    }

    #[test]
    fn intersecting_translated_shape_with_ray() {
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let expected_ray = Ray::new(Point::new(-5.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let mut s = TestShape::new(expected_ray);
        let reg = Registry::new();
        s.set_transform(Matrix4::translation(5.0, 0.0, 0.0));
        s.intersect(&reg, &r);
    }

    #[test]
    fn computing_normal_on_translated_shape() {
        let registry = Rc::new(RefCell::new(Registry::new()));
        let expected_ray = Ray::new(Point::new(-5.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let id;
        {
            let mut reg = registry.borrow_mut();
            let mut s = Box::from(TestShape::new(expected_ray));
            s.set_transform(Matrix4::translation(0.0, 1.0, 0.0));
            id = reg.register(s);
        }
        let reg = registry.borrow();
        let s = reg.get(id).unwrap();
        assert_eq!(
            s.normal(&reg, Point::new(0.0, 1.70711, -0.70711)),
            Vector::new(0.0, 0.70711, -0.70711)
        );
    }

    #[test]
    fn computing_normal_on_transformed_shape() {
        let registry = Rc::new(RefCell::new(Registry::new()));
        let expected_ray = Ray::new(Point::new(-5.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let id;
        {
            let mut reg = registry.borrow_mut();
            let mut s = Box::from(TestShape::new(expected_ray));
            s.set_transform(Matrix4::scaling(1.0, 0.5, 1.0) * Matrix4::rotation_z(PI / 5.0));
            id = reg.register(s);
        }
        let reg = registry.borrow();
        let s = reg.get(id).unwrap();
        assert_eq!(
            s.normal(&reg, Point::new(0.0, 2_f32.sqrt() / 2.0, -2_f32.sqrt() / 2.0)),
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
            let mut reg = Registry::new();
            let id = reg.register(Box::from(Plane::new()));
            let p = reg.get(id).unwrap();
            let r = Ray::new(Point::new(0.0, 10.0, 0.0), Vector::new(0.0, 0.0, 1.0));
            assert_eq!(
                p.local_intersect(&reg, &r)
                    .iter()
                    .map(|x| x.t)
                    .collect::<Vec<f32>>(),
                vec!()
            );
        }
        #[test]
        fn intersect_with_coplanar_ray() {
            let mut reg = Registry::new();
            let id = reg.register(Box::from(Plane::new()));
            let p = reg.get(id).unwrap();
            let r = Ray::new(Point::origin(), Vector::new(0.0, 0.0, 1.0));
            assert_eq!(
                p.local_intersect(&reg, &r)
                    .iter()
                    .map(|x| x.t)
                    .collect::<Vec<f32>>(),
                vec!()
            );
        }

        #[test]
        fn intersect_with_ray_from_above() {
            let mut reg = Registry::new();
            let id = reg.register(Box::from(Plane::new()));
            let p = reg.get(id).unwrap();
            let r = Ray::new(Point::new(0.0, 1.0, 0.0), Vector::new(0.0, -1.0, 0.0));
            assert_eq!(
                p.local_intersect(&reg, &r)
                    .iter()
                    .map(|x| x.t)
                    .collect::<Vec<f32>>(),
                vec!(1.0)
            );
        }

        #[test]
        fn plane_intersect_with_ray_from_below() {
            let mut reg = Registry::new();
            let id = reg.register(Box::from(Plane::new()));
            let p = reg.get(id).unwrap();
            let r = Ray::new(Point::new(0.0, -1.0, 0.0), Vector::new(0.0, 1.0, 0.0));
            assert_eq!(
                p.local_intersect(&reg, &r)
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
            let registry = Rc::new(RefCell::new(Registry::new()));
            let mut sph = Box::from(Sphere::new());
            sph.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
            let mut sph2 = Box::from(Sphere::new());
            sph2.set_transform(Matrix4::translation(5.0, 0.0, 0.0));
            let id;
            let id2;
            {
                let mut reg = registry.borrow_mut();
                id = reg.register(sph);
                id2 = reg.register(sph2);
            }
            let reg = registry.borrow();
            let s = reg.get(id).unwrap();
            let transformed = reg.get(id2).unwrap();

            assert_eq!(
                s.intersect(&reg, &r).iter().map(|x| x.t).collect::<Vec<f32>>(),
                vec!(3.0, 7.0)
            );

            assert_eq!(
                transformed.intersect(&reg, &r).iter().map(|x| x.t).collect::<Vec<f32>>(),
                vec!()
            );
        }

        #[test]
        fn normal() {
            let registry = Rc::new(RefCell::new(Registry::new()));
            let id;
            let tid;
            let srid;
            {
                let mut reg = registry.borrow_mut();
                id = reg.register(Box::from(Sphere::new()));
                let mut translated = Box::from(Sphere::new());
                let mut scaled_and_rotated = Box::from(Sphere::new());
                let translate = Matrix4::translation(0.0, 1.0, 0.0);
                let scale_and_rotate = Matrix4::scaling(1.0, 0.5, 1.0) * Matrix4::rotation_z(PI / 5.0);
                translated.set_transform(translate);
                scaled_and_rotated.set_transform(scale_and_rotate);

                tid = reg.register(translated);
                srid = reg.register(scaled_and_rotated);
            }
            let reg = registry.borrow();
            let s = reg.get(id).unwrap();
            assert_eq!(
                s.normal(&reg, Point::new(1.0, 0.0, 0.0)),
                Vector::new(1.0, 0.0, 0.0)
            );
            assert_eq!(
                s.normal(&reg, Point::new(0.0, 1.0, 0.0)),
                Vector::new(0.0, 1.0, 0.0)
            );
            assert_eq!(
                s.normal(&reg, Point::new(0.0, 0.0, 1.0)),
                Vector::new(0.0, 0.0, 1.0)
            );
            let v = s.normal(&reg, Point::new(
                3_f32.sqrt() / 3.0,
                3_f32.sqrt() / 3.0,
                3_f32.sqrt() / 3.0,
            ));
            assert_eq!(v, v.normalize());

            let t = reg.get(tid).unwrap();
            assert_eq!(
                t.normal(&reg, Point::new(0.0, 1.70711, -0.70711)),
                Vector::new(0.0, 0.70711, -0.70711)
            );

            let sr = reg.get(srid).unwrap();
            assert_eq!(
                sr.normal(&reg, Point::new(0.0, 2_f32.sqrt() / 2.0, -2_f32.sqrt() / 2.0)),
                Vector::new(0.0, 0.970140, -0.24254)
            );
        }

        #[test]
        fn intersection() {
            let mut r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
            let sph = Box::from(Sphere::new());
            let mut reg = Registry::new();
            let id = reg.register(sph);
            let s = reg.get(id).unwrap();

            let mut is = s.intersect(&reg, &r);
            assert_eq!(is.iter().map(|x| x.t).collect::<Vec<f32>>(), vec!(4.0, 6.0));

            r = Ray::new(Point::new(0.0, 1.0, -5.0), Vector::new(0.0, 0.0, 1.0));
            is = s.intersect(&reg, &r);
            assert_eq!(is.iter().map(|x| x.t).collect::<Vec<f32>>(), vec!(5.0, 5.0));

            r = Ray::new(Point::new(0.0, 2.0, -5.0), Vector::new(0.0, 0.0, 1.0));
            assert_eq!(
                s.intersect(&reg, &r).iter().map(|x| x.t).collect::<Vec<f32>>(),
                vec!()
            );

            r = Ray::new(Point::origin(), Vector::new(0.0, 0.0, 1.0));
            is = s.intersect(&reg, &r);
            assert_eq!(
                is.iter().map(|x| x.t).collect::<Vec<f32>>(),
                vec!(-1.0, 1.0)
            );

            r = Ray::new(Point::new(0.0, 0.0, 5.0), Vector::new(0.0, 0.0, 1.0));
            is = s.intersect(&reg, &r);
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
            let mut reg = Registry::new();
            let id = reg.register(Box::from(Cube::new()));
            let c = reg.get(id).unwrap();

            for (origin, direction, t1, t2) in &examples {
                let r = Ray::new(*origin, *direction);
                let hits = c
                    .local_intersect(&reg, &r)
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
            let reg = Registry::new();

            for (origin, direction) in &examples {
                let r = Ray::new(*origin, *direction);
                let hits = c
                    .local_intersect(&reg, &r)
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
            let mut reg = Registry::new();
            let id = reg.register(Box::from(Cylinder::new()));
            let c = reg.get(id).unwrap();

            for (origin, direction, t1, t2) in &examples {
                let r = Ray::new(*origin, direction.normalize());
                let hits = c
                    .local_intersect(&reg, &r)
                    .iter()
                    .map(|x| x.t)
                    .collect::<Vec<f32>>();
                assert_eq!(hits, [*t1, *t2]);
            }
        }

        #[test]
        fn ray_misses_cylinder() {
            let examples = [
                (Point::new(1.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
                (Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
                (Point::new(0.0, 0.0, -5.0), Vector::new(1.0, 1.0, 1.0)),
            ];
            let mut reg = Registry::new();
            let id = reg.register(Box::from(Cylinder::new()));
            let c = reg.get(id).unwrap();

            for (origin, direction) in &examples {
                let r = Ray::new(*origin, direction.normalize());
                let hits = c
                    .local_intersect(&reg, &r)
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
            let examples: &[(Point, Vector, usize)] = &[
                // origin, direction, intersection count
                (Point::new(0.0, 1.5, 0.0), Vector::new(0.1, 1.0, 0.0), 0),
                (Point::new(0.0, 3.0, -5.0), Vector::new(0.0, 0.0, 1.0), 0),
                (Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0), 0),
                (Point::new(0.0, 2.0, -5.0), Vector::new(0.0, 0.0, 1.0), 0),
                (Point::new(0.0, 1.0, -5.0), Vector::new(0.0, 0.0, 1.0), 0),
                (Point::new(0.0, 1.5, -2.0), Vector::new(0.0, 0.0, 1.0), 2),
            ];
            let mut reg = Registry::new();
            let id = reg.register(Box::from(Cylinder::open(1.0, 2.0)));
            let c = reg.get(id).unwrap();

            for (origin, direction, hit_count) in examples {
                let r = Ray::new(*origin, direction.normalize());
                let hits = c
                    .local_intersect(&reg, &r)
                    .iter()
                    .map(|x| x.t)
                    .collect::<Vec<f32>>();
                assert_eq!(hits.len(), *hit_count);
            }
        }

        #[test]
        fn intersecting_caps_of_closed_cylinder() {
            let examples: &[(Point, Vector, usize)] = &[
                // point, direction, intersection count
                (Point::new(0.0, 3.0, 0.0), Vector::new(0.0, -1.0, 0.0), 2),
                (Point::new(0.0, 3.0, -2.0), Vector::new(0.0, -1.0, 2.0), 2),
                (Point::new(0.0, 4.01, -2.0), Vector::new(0.0, -1.0, 1.0), 2), // edge case
                (Point::new(0.0, 0.0, -2.0), Vector::new(0.0, 1.0, 2.0), 2),
                (Point::new(0.0, -1.01, -2.0), Vector::new(0.0, 1.0, 1.0), 2), // edge case
            ];
            let mut reg = Registry::new();
            let id = reg.register(Box::from(Cylinder::closed(1.0, 2.0)));
            let c = reg.get(id).unwrap();

            for (origin, direction, hit_count) in examples {
                let r = Ray::new(*origin, direction.normalize());
                let hits = c
                    .local_intersect(&reg, &r)
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
            let examples = [
                // origin, direction, t1, t2
                (
                    Point::new(0.0, 0.0, -5.0),
                    Vector::new(0.0, 0.0, 1.0),
                    5.0,
                    5.0,
                ),
                (
                    Point::new(0.0, 0.0, -5.00001),
                    Vector::new(1.0, 1.0, 1.0),
                    8.660272,
                    8.660272,
                ),
                (
                    Point::new(1.0, 1.0, -5.0),
                    Vector::new(-0.5, -1.0, 1.0),
                    4.5500546,
                    49.449955,
                ),
            ];
            let mut reg = Registry::new();
            let id = reg.register(Box::from(Cone::new()));
            let c = reg.get(id).unwrap();

            for (origin, direction, t1, t2) in &examples {
                let r = Ray::new(*origin, direction.normalize());
                let hits = c
                    .local_intersect(&reg, &r)
                    .iter()
                    .map(|x| x.t)
                    .collect::<Vec<f32>>();
                assert_eq!(hits, [*t1, *t2]);
            }
        }

        #[test]
        fn ray_strikes_cone_parallel_to_one_of_its_halves() {
            let origin = Point::new(0.0, 0.0, -1.0);
            let direction = Vector::new(0.0, 1.0, 1.0).normalize();
            let r = Ray::new(origin, direction);
            let mut reg = Registry::new();
            let id = reg.register(Box::from(Cone::new()));
            let c = reg.get(id).unwrap();

            let hits = c
                .local_intersect(&reg, &r)
                .iter()
                .map(|x| x.t)
                .collect::<Vec<f32>>();
            assert_eq!(hits, [0.35355338]);
        }

        #[test]
        fn intersecting_cone_end_caps() {
            let examples = [
                // origin, direction, hit_count
                (Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 1.0, 0.0), 0),
                (Point::new(0.0, 0.0, -0.25), Vector::new(0.0, 1.0, 1.0), 2),
                (Point::new(0.0, 0.0, -0.25), Vector::new(0.0, 1.0, 0.0), 4),
            ];
            let mut reg = Registry::new();
            let id = reg.register(Box::from(Cone::closed(-0.5, 0.5)));
            let c = reg.get(id).unwrap();

            for (origin, direction, hit_count) in &examples {
                let r = Ray::new(*origin, direction.normalize());
                let hits = c
                    .local_intersect(&reg, &r)
                    .iter()
                    .map(|x| x.t)
                    .collect::<Vec<f32>>();
                assert_eq!(hits.len(), *hit_count);
            }
        }

        #[test]
        fn computing_normal_vector_on_cone() {
            let c = Cone::new();
            let examples = [
                // point, normal
                (Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 0.0)),
                (Point::new(1.0, 1.0, 1.0), Vector::new(1.0, -(2_f32.sqrt()), 1.0)),
                (Point::new(-1.0, -1.0, 0.0), Vector::new(-1.0, 1.0, 0.0)),
            ];

            for (point, direction) in &examples {
                assert_eq!(c.local_normal_at(&point), *direction);
            }
        }
    }

    mod group {
        use super::*;

        #[test]
        fn intersecting_ray_empty_group() {
            let r = Ray::new(Point::origin(), Vector::new(0.0, 0.0, 1.0));
            let mut reg = Registry::new();
            let id = reg.register(Box::from(Group::new()));
            let g = reg.get(id).unwrap();
            let hits = g
                .local_intersect(&reg, &r)
                .iter()
                .map(|x| x.t)
                .collect::<Vec<f32>>();
            assert_eq!(hits, vec![]);
        }

        #[test]
        fn intersecting_ray_non_empty_group() {
            let mut reg = Registry::new();

            let g = Group::new();
            let s1 = Sphere::new();
            let mut s2 = Sphere::new();
            s2.set_transform(Matrix4::translation(0.0, 0.0, -3.0));
            let mut s3 = Sphere::new();
            s3.set_transform(Matrix4::translation(5.0, 0.0, 0.0));

            let s1id = reg.register(Box::from(s1));
            let s2id = reg.register(Box::from(s2));
            let s3id = reg.register(Box::from(s3));
            let gid = reg.register(Box::from(g));

            reg.add(gid, s1id);
            reg.add(gid, s2id);
            reg.add(gid, s3id);

            let group = reg.get(gid).unwrap();

            let s2uid = reg.get(s2id).unwrap().id();
            let s1uid = reg.get(s1id).unwrap().id();

            let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
            let hits = group
                .local_intersect(&reg, &r)
                .iter()
                .map(|x| x.object.id())
                .collect::<Vec<i32>>();
            assert_eq!(hits, vec![s2uid, s2uid, s1uid, s1uid]);
        }

        #[test]
        fn intersecting_transformed_group() {
            let mut reg = Registry::new();
            let id = reg.register(Box::from(Group::new()));
            let mg = reg.get_mut(id).unwrap();
            mg.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
            let sid = reg.register(Box::from(Sphere::new()));
            let ms = reg.get_mut(sid).unwrap();
            ms.set_transform(Matrix4::translation(5.0, 0.0, 0.0));
            reg.add(id, sid);
            let r = Ray::new(Point::new(10.0, 0.0, -10.0), Vector::new(0.0, 0.0, 1.0));

            let g = reg.get(id).unwrap();

            let hits = g
                .intersect(&reg, &r)
                .iter()
                .map(|x| x.object.id())
                .collect::<Vec<i32>>();
            assert_eq!(hits.len(), 2);
        }

        #[test]
        fn converting_point_from_world_to_object_space() {
            let mut reg = Registry::new();

            let id = reg.register(Box::from(Group::new()));
            let mg = reg.get_mut(id).unwrap();
            mg.set_transform(Matrix4::rotation_y(PI/2.0));

            let id2 = reg.register(Box::from(Group::new()));
            let mg2 = reg.get_mut(id2).unwrap();
            mg2.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));

            let sid = reg.register(Box::from(Sphere::new()));
            let ms = reg.get_mut(sid).unwrap();
            ms.set_transform(Matrix4::translation(5.0, 0.0, 0.0));

            reg.add(id, id2);
            reg.add(id2, sid);

            let s = reg.get(sid).unwrap();
            let p = s.world_to_object(&reg, &Point::new(-2.0, 0.0, -10.0));
            assert_eq!(p, Point::new(0.0, 0.0, -1.0));
        }

        #[test]
        fn converting_normal_vector_from_object_to_world_space() {
            let mut reg = Registry::new();

            let id = reg.register(Box::from(Group::new()));
            let mg = reg.get_mut(id).unwrap();
            mg.set_transform(Matrix4::rotation_y(PI/2.0));

            let id2 = reg.register(Box::from(Group::new()));
            let mg2 = reg.get_mut(id2).unwrap();
            mg2.set_transform(Matrix4::scaling(1.0, 2.0, 3.0));

            let sid = reg.register(Box::from(Sphere::new()));
            let ms = reg.get_mut(sid).unwrap();
            ms.set_transform(Matrix4::translation(5.0, 0.0, 0.0));

            reg.add(id, id2);
            reg.add(id2, sid);

            let s = reg.get(sid).unwrap();
            let v = s.normal_to_world(&reg, &Vector::new(3_f32.sqrt()/3.0, 3_f32.sqrt()/3.0, 3_f32.sqrt()/3.0));
            assert_eq!(v, Vector::new(0.2857, 0.4286, -0.8571));
        }

        #[test]
        fn finding_normal_on_child_object() {
            let mut reg = Registry::new();

            let id = reg.register(Box::from(Group::new()));
            let mg = reg.get_mut(id).unwrap();
            mg.set_transform(Matrix4::rotation_y(PI/2.0));

            let id2 = reg.register(Box::from(Group::new()));
            let mg2 = reg.get_mut(id2).unwrap();
            mg2.set_transform(Matrix4::scaling(1.0, 2.0, 3.0));

            let sid = reg.register(Box::from(Sphere::new()));
            let ms = reg.get_mut(sid).unwrap();
            ms.set_transform(Matrix4::translation(5.0, 0.0, 0.0));

            reg.add(id, id2);
            reg.add(id2, sid);

            let s = reg.get(sid).unwrap();
            let v = s.normal(&reg, Point::new(1.7321, 1.1547, -5.5774));
            assert_eq!(v, Vector::new(0.2857, 0.4286, -0.8571));
        }
    }
}
