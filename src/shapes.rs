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

fn search_material(reg: &Registry, id: NodeId) -> &Material {
    let mut mat = Material::base();
    if let Some(p) = reg.parent(id) {
        mat = p.material(reg)
    }
    mat
}

#[derive(Debug, Copy, Clone)]
pub struct Bounds(Point, Point);

impl Bounds {
    fn union(&self, other: &Bounds) -> Bounds {
        Bounds(
            Point::new(
                self.0.x.min(other.0.x),
                self.0.y.min(other.0.y),
                self.0.z.min(other.0.z),
            ),
            Point::new(
                self.0.x.max(other.0.x),
                self.0.y.max(other.0.y),
                self.0.z.max(other.0.z),
            ),
        )
    }

    fn transform(&self, t: &Matrix4) -> Bounds {
        let vertices = vec![
            self.0,
            Point::new(self.0.x, self.0.y, self.1.z),
            Point::new(self.0.x, self.1.y, self.0.z),
            Point::new(self.0.x, self.1.y, self.1.z),
            Point::new(self.1.x, self.0.y, self.0.z),
            Point::new(self.1.x, self.0.y, self.1.z),
            Point::new(self.1.x, self.1.y, self.0.z),
            self.1,
        ];
        let mut out = Bounds(
            Point::new(INFINITY, INFINITY, INFINITY),
            Point::new(-INFINITY, -INFINITY, -INFINITY)
        );
        for v in vertices {
            let p = *t * v;
            out = out.union(&Bounds(p, p));
        }
        out
    }

    fn check_axis(origin: f32, invdir: f32, min: f32, max: f32) -> (f32, f32) {
        let tmin;
        let tmax;
        if invdir >= 0.0 {
            tmin = (min - origin) * invdir;
            tmax = (max - origin) * invdir;
        } else {
			tmax = (min - origin) * invdir;
			tmin = (max - origin) * invdir;
        }
        (tmin, tmax)
    }

    fn intersects_inverse(&self, r: &Ray) -> bool {
        let (xtmin, xtmax) = Bounds::check_axis(r.origin.x, r.direction.x, self.0.x, self.1.x);
        let (ytmin, ytmax) = Bounds::check_axis(r.origin.y, r.direction.y, self.0.y, self.1.y);
        let (ztmin, ztmax) = Bounds::check_axis(r.origin.z, r.direction.z, self.0.z, self.1.z);
        let mut tmin = -INFINITY;
        let mut tmax = INFINITY;
        if xtmin > tmin {
            tmin = xtmin;
        }
        if ytmin > tmin {
            tmin = ytmin;
        }
        if ztmin > tmin {
            tmin = ztmin;
        }
        if xtmax < tmax {
            tmax = xtmax;
        }
        if ytmax < tmax {
            tmax = ytmax;
        }
        if ztmax < tmax {
            tmax = ztmax;
        }

        tmin <= tmax
    }

    fn intersects(&self, r: &Ray) -> bool {
        let reverted = Ray::new(r.origin, Vector::new(1.0/r.direction.x, 1.0/r.direction.y, 1.0/r.direction.z));
        self.intersects_inverse(&reverted)
    }
}

pub trait Shape: Send + Sync {
    fn id(&self) -> i32;
    fn bounds(&self) -> &Bounds;
    fn set_bounds(&mut self, bounds: Bounds);
    fn set_tag(&mut self, id: NodeId);
    fn tag(&self) -> NodeId;
    fn casts_shadows(&self) -> bool;
    fn transform(&self) -> &Matrix4;
    fn inverse_transform(&self) -> &Matrix4;
    fn set_transform(&mut self, t: Matrix4);
    fn material<'a>(&'a self, reg: &'a Registry) -> &'a Material;
    fn set_material(&mut self, m: Material);
    fn normal(&self, reg: &Registry, world_point: Point, i: &Intersection) -> Vector {
        let local_point = self.world_to_object(&reg, &world_point);
        let local_normal = self.local_normal_at(&local_point, i);
        self.normal_to_world(&reg, &local_normal)
    }
    fn local_normal_at(&self, p: &Point, i: &Intersection) -> Vector;
    fn local_intersect<'a>(&'a self, reg: &'a Registry, ray: &Ray) -> Vec<Intersection>;
    fn intersect<'a>(&'a self, reg: &'a Registry, r: &Ray) -> Vec<Intersection> {
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
    pub bounds: Bounds,
    pub material: Option<Material>,
    pub casts_shadows: bool,
    pub inverse_transform: Matrix4,
}

impl Sphere {
    pub fn new() -> Self {
        Sphere {
            id: gen_id(),
            tag: None,
            material: None,
            transform: Matrix4::id(),
            casts_shadows: true,
            bounds: Bounds(Point::new(-1.0, -1.0, -1.0), Point::new(1.0, 1.0, 1.0)),
            inverse_transform: Matrix4::id().inverse(),
        }
    }

    pub fn glass() -> Self {
        Sphere {
            id: gen_id(),
            tag: None,
            material: Some(Material {
                transparency: 1.0,
                refractive_index: 1.5,
                ..Material::default()
            }),
            transform: Matrix4::id(),
            casts_shadows: true,
            bounds: Bounds(Point::new(-1.0, -1.0, -1.0), Point::new(1.0, 1.0, 1.0)),
            inverse_transform: Matrix4::id().inverse(),
        }
    }
}

impl Shape for Sphere {
    fn id(&self) -> i32 {
        self.id
    }
    fn bounds(&self) -> &Bounds {
        &self.bounds
    }
    fn set_bounds(&mut self, bounds: Bounds) {
        self.bounds = bounds;
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
    fn material<'a>(&'a self, reg: &'a Registry) -> &'a Material {
        if let Some(mat) = &self.material {
            &mat
        } else {
            search_material(reg, self.tag())
        }
    }
    fn set_material(&mut self, m: Material) {
        self.material = Some(m);
    }

    fn local_normal_at(&self, p: &Point, _i: &Intersection) -> Vector {
        *p - Point::origin()
    }

    fn local_intersect<'a>(&'a self, _reg: &'a Registry, ray: &Ray) -> Vec<Intersection> {
        let sphere_to_ray = ray.origin - Point::origin();
        let a = ray.direction.dot(&ray.direction);
        let b = ray.direction.dot(&sphere_to_ray) * 2.0;
        let c = sphere_to_ray.dot(&sphere_to_ray) - 1.0;
        let discriminant = b.powi(2) - 4.0 * a * c;

        let s = self.tag();

        if discriminant < 0.0 {
            vec![]
        } else {
            let t1 = (-b - discriminant.sqrt()) / (a * 2.0);
            let t2 = (-b + discriminant.sqrt()) / (a * 2.0);
            if t1 > t2 {
                vec![
                    Intersection { uv: None, t: t2, object: s },
                    Intersection { uv: None, t: t1, object: s },
                ]
            } else {
                vec![
                    Intersection { uv: None, t: t1, object: s },
                    Intersection { uv: None, t: t2, object: s },
                ]
            }
        }
    }
}

pub struct Plane {
    pub id: i32,
    pub tag: Option<NodeId>,
    pub transform: Matrix4,
    pub material: Option<Material>,
    pub casts_shadows: bool,
    pub bounds: Bounds,
    inverse_transform: Matrix4,
}

impl Plane {
    pub fn new() -> Plane {
        Plane {
            id: gen_id(),
            tag: None,
            transform: Matrix4::id(),
            material: None,
            casts_shadows: true,
            bounds: Bounds(
            Point::new(-INFINITY, -INFINITY, -INFINITY),
            Point::new(INFINITY, INFINITY, INFINITY)
            ),
            inverse_transform: Matrix4::id().inverse(),
        }
    }
}

impl Shape for Plane {
    fn id(&self) -> i32 {
        self.id
    }
    fn bounds(&self) -> &Bounds {
        &self.bounds
    }
    fn set_bounds(&mut self, bounds: Bounds) {
        self.bounds = bounds;
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
        self.material = Some(m);
    }
    fn material<'a>(&'a self, reg: &'a Registry) -> &'a Material {
        if let Some(mat) = &self.material {
            &mat
        } else {
            search_material(reg, self.tag())
        }
    }
    fn transform(&self) -> &Matrix4 {
        &self.transform
    }
    fn inverse_transform(&self) -> &Matrix4 {
        &self.inverse_transform
    }
    fn local_normal_at(&self, _p: &Point, _i: &Intersection) -> Vector {
        Vector::new(0.0, 1.0, 0.0)
    }
    fn local_intersect<'a>(&'a self, _reg: &'a Registry, ray: &Ray) -> Vec<Intersection> {
        if ray.direction.y.abs() < EPSILON {
            vec![]
        } else {
            let t = -(ray.origin.y) / ray.direction.y;
            vec![Intersection { uv: None, t, object: self.tag() }]
        }
    }
}

pub struct Cube {
    pub id: i32,
    pub tag: Option<NodeId>,
    pub transform: Matrix4,
    pub material: Option<Material>,
    pub casts_shadows: bool,
    pub bounds: Bounds,
    inverse_transform: Matrix4,
}

impl Cube {
    pub fn new() -> Cube {
        Cube {
            id: gen_id(),
            tag: None,
            transform: Matrix4::id(),
            material: None,
            casts_shadows: true,
            bounds: Bounds(
                Point::new(-1.0, -1.0, -1.0),
                Point::new(1.0, 1.0, 1.0)
            ),
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
    fn bounds(&self) -> &Bounds {
        &self.bounds
    }
    fn set_bounds(&mut self, bounds: Bounds) {
        self.bounds = bounds;
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
        self.material = Some(m);
    }
    fn material<'a>(&'a self, reg: &'a Registry) -> &'a Material {
        if let Some(mat) = &self.material {
            &mat
        } else {
            search_material(reg, self.tag())
        }
    }
    fn transform(&self) -> &Matrix4 {
        &self.transform
    }
    fn inverse_transform(&self) -> &Matrix4 {
        &self.inverse_transform
    }
    fn local_normal_at(&self, p: &Point, _i: &Intersection) -> Vector {
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
    fn local_intersect<'a>(&'a self, _reg: &'a Registry, ray: &Ray) -> Vec<Intersection> {
        let (xtmin, xtmax) = check_axis(ray.origin.x, ray.direction.x);
        let (ytmin, ytmax) = check_axis(ray.origin.y, ray.direction.y);
        let (ztmin, ztmax) = check_axis(ray.origin.z, ray.direction.z);

        let tmin = xtmin.max(ytmin).max(ztmin);
        let tmax = xtmax.min(ytmax).min(ztmax);

        if tmin > tmax {
            vec![]
        } else {
            let shape = self.tag();
            vec![
                Intersection { uv: None,
                    t: tmin,
                    object: shape,
                },
                Intersection { uv: None,
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
    pub material: Option<Material>,
    pub casts_shadows: bool,
    pub minimum: f32,
    pub maximum: f32,
    pub closed: bool,
    pub bounds: Bounds,
    inverse_transform: Matrix4,
}

impl Cylinder {
    pub fn new() -> Cylinder {
        Cylinder {
            id: gen_id(),
            tag: None,
            transform: Matrix4::id(),
            material: None,
            casts_shadows: true,
            minimum: -INFINITY,
            maximum: INFINITY,
            bounds: Bounds(Point::new(-1.0, -INFINITY, -1.0), Point::new(1.0, INFINITY, 1.0)),
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

    fn intersect_caps<'a>(&'a self, ray: &Ray) -> Vec<Intersection> {
        let mut is = vec![];
        if !self.closed || ray.direction.y.abs() < EPSILON {
            return is;
        }
        let shape = self.tag();
        let lower_t = (self.minimum - ray.origin.y) / ray.direction.y;
        if Self::check_cap(ray, lower_t) {
            is.push(Intersection { uv: None,
                t: lower_t,
                object: shape,
            })
        }

        let upper_t = (self.maximum - ray.origin.y) / ray.direction.y;
        if Self::check_cap(ray, upper_t) {
            is.push(Intersection { uv: None,
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
    fn bounds(&self) -> &Bounds {
        &self.bounds
    }
    fn set_bounds(&mut self, bounds: Bounds) {
        self.bounds = bounds;
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
        self.material = Some(m);
    }
    fn material<'a>(&'a self, reg: &'a Registry) -> &'a Material {
        if let Some(mat) = &self.material {
            &mat
        } else {
            search_material(reg, self.tag())
        }
    }
    fn transform(&self) -> &Matrix4 {
        &self.transform
    }
    fn inverse_transform(&self) -> &Matrix4 {
        &self.inverse_transform
    }
    fn local_normal_at(&self, p: &Point, _i: &Intersection) -> Vector {
        let dist = p.x.powi(2) + p.z.powi(2);
        if dist < 1.0 && p.y >= self.maximum - EPSILON {
            Vector::new(0.0, 1.0, 0.0)
        } else if dist < 1.0 && p.y <= self.minimum + EPSILON {
            Vector::new(0.0, -1.0, 0.0)
        } else {
            Vector::new(p.x, 0.0, p.z)
        }
    }
    fn local_intersect<'a>(&'a self, _reg: &'a Registry, ray: &Ray) -> Vec<Intersection> {
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
                let shape = self.tag();
                let y0 = ray.origin.y + t0 * ray.direction.y;
                if self.minimum < y0 && y0 < self.maximum {
                    intersections.push(Intersection { uv: None,
                        t: t0,
                        object: shape,
                    });
                }
                let y1 = ray.origin.y + t1 * ray.direction.y;
                if self.minimum < y1 && y1 < self.maximum {
                    intersections.push(Intersection { uv: None,
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
    pub tag: Option<NodeId>,
    pub transform: Matrix4,
    pub material: Option<Material>,
    pub casts_shadows: bool,
    pub minimum: f32,
    pub maximum: f32,
    pub closed: bool,
    pub bounds: Bounds,
    inverse_transform: Matrix4,
}

impl Cone {
    fn calculate_bounds(minimum: f32, maximum: f32) -> Bounds {
        let max = minimum.abs().max(maximum.abs());
        Bounds(Point::new(-max, minimum, -max), Point::new(max, maximum, max))
    }

    pub fn new() -> Cone {
        Cone {
            id: gen_id(),
            tag: None,
            transform: Matrix4::id(),
            material: None,
            casts_shadows: true,
            minimum: -INFINITY,
            maximum: INFINITY,
            closed: false,
            bounds: Cone::calculate_bounds(-INFINITY, INFINITY),
            inverse_transform: Matrix4::id().inverse(),
        }
    }

    pub fn open(minimum: f32, maximum: f32) -> Cone {
        let mut c = Cone::new();
        c.minimum = minimum;
        c.maximum = maximum;
        c.bounds = Cone::calculate_bounds(minimum, maximum);
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
    fn bounds(&self) -> &Bounds {
        &self.bounds
    }
    fn set_bounds(&mut self, bounds: Bounds) {
        self.bounds = bounds;
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
        self.material = Some(m);
    }
    fn material<'a>(&'a self, reg: &'a Registry) -> &'a Material {
        if let Some(mat) = &self.material {
            &mat
        } else {
            search_material(reg, self.tag())
        }
    }
    fn transform(&self) -> &Matrix4 {
        &self.transform
    }
    fn inverse_transform(&self) -> &Matrix4 {
        &self.inverse_transform
    }
    fn local_normal_at(&self, p: &Point, _i: &Intersection) -> Vector {
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
    fn local_intersect<'a>(&'a self, _reg: &'a Registry, ray: &Ray) -> Vec<Intersection> {
        let o = ray.origin;
        let d = ray.direction;
        let this = self.tag();
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
                    intersections.push(Intersection { uv: None,
                        t: t0,
                        object: this,
                    });
                }
                let y1 = ray.origin.y + t1 * ray.direction.y;
                if self.minimum < y1 && y1 < self.maximum {
                    intersections.push(Intersection { uv: None,
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
                    intersections.push(Intersection { uv: None,
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
                intersections.push(Intersection { uv: None,
                    t: t0,
                    object: this,
                });
            }

            let t1 = (self.maximum - ray.origin.y) / ray.direction.y;
            if Self::check_cap(ray, self.maximum.abs(), t1) {
                intersections.push(Intersection { uv: None,
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
    pub material: Option<Material>,
    pub casts_shadows: bool,
    inverse_transform: Matrix4,
    bounds: Bounds,
}

impl Group {
    pub fn new() -> Group {
        let transform = Matrix4::id();
        Group {
            id: gen_id(),
            tag: None,
            transform: transform,
            inverse_transform: transform.inverse(),
            material: None,
            casts_shadows: true,
            bounds: Bounds(
                Point::new(-INFINITY, -INFINITY, -INFINITY),
                Point::new(INFINITY, INFINITY, INFINITY),
                )
        }
    }
}

pub fn precompute_bounds(children: Vec<&Box<Shape>>) -> Bounds {
    let mut bounds = Bounds(
            Point::new(INFINITY, INFINITY, INFINITY),
            Point::new(-INFINITY, -INFINITY, -INFINITY),
            );
    for child in children.iter() {
        bounds = bounds.union(&child.bounds().transform(&child.transform()));
    }
    println!("New precomputed bounds: {:?}", bounds);
    bounds
}

impl Shape for Group {
    fn id(&self) -> i32 {
        self.id
    }
    fn bounds(&self) -> &Bounds {
        &self.bounds
    }
    fn set_bounds(&mut self, bounds: Bounds) {
        self.bounds = bounds;
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
        self.material = Some(m);
    }
    fn material<'a>(&'a self, reg: &'a Registry) -> &'a Material {
        if let Some(mat) = &self.material {
            &mat
        } else {
            search_material(reg, self.tag())
        }
    }
    fn transform(&self) -> &Matrix4 {
        &self.transform
    }
    fn inverse_transform(&self) -> &Matrix4 {
        &self.inverse_transform
    }
    fn local_normal_at(&self, _p: &Point, _i: &Intersection) -> Vector {
        panic!("Calling local_normal_at in group");
    }
    fn local_intersect<'a>(&'a self, reg: &'a Registry, ray: &Ray) -> Vec<Intersection> {
        if !self.bounds().intersects(&ray) { return vec![] };

        let children = reg.children(self.tag());
        let mut is = vec![];
        for child in children {
            is.append(&mut child.intersect(&reg, &ray));
        }
        is.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
        is
    }
}

pub struct Triangle {
    a: Point,
    b: Point,
    c: Point,
    e1: Vector,
    e2: Vector,
    pub normal: Vector,
    tag: Option<NodeId>,
    id: i32,
    casts_shadows: bool,
    transform: Matrix4,
    pub bounds: Bounds,
    inverse_transform: Matrix4,
    material: Option<Material>,
}

impl Triangle {
    pub fn new(a: Point, b: Point, c: Point) -> Triangle {
        let e1 = b - a;
        let e2 = c - a;
        let normal = e2.cross(&e1).normalize();

        let bounds = Bounds(
            Point::new(
                a.x.min(b.x).min(c.x),
                a.y.min(b.y).min(c.y),
                a.z.min(b.z).min(c.z),
            ),
            Point::new(
                a.x.max(b.x).max(c.x),
                a.y.max(b.y).max(c.y),
                a.z.max(b.z).max(c.z),
            ),
        );

        Triangle {
            a, b, c, e1, e2, normal,
            id: 0, //gen_id(),
            casts_shadows: true,
            transform: Matrix4::id(),
            tag: None,
            bounds: bounds,
            inverse_transform: Matrix4::id().inverse(),
            material: None,
            }
    }
}

impl Shape for Triangle {
    fn id(&self) -> i32 {
        self.id
    }
    fn bounds(&self) -> &Bounds {
        &self.bounds
    }
    fn set_bounds(&mut self, bounds: Bounds) {
        self.bounds = bounds;
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
        self.material = Some(m)
    }
    fn material<'a>(&'a self, reg: &'a Registry) -> &'a Material {
        if let Some(mat) = &self.material {
            &mat
        } else {
            search_material(reg, self.tag())
        }
    }
    fn transform(&self) -> &Matrix4 {
        &self.transform
    }
    fn inverse_transform(&self) -> &Matrix4 {
        &self.inverse_transform
    }
    fn local_normal_at(&self, _p: &Point, _i: &Intersection) -> Vector {
        self.normal
    }
    fn local_intersect<'a>(&'a self, _reg: &'a Registry, ray: &Ray) -> Vec<Intersection> {
        let dir_cross_e2 = ray.direction.cross(&self.e2);
        let det = self.e1.dot(&dir_cross_e2);
        if det.abs() < EPSILON {
            vec![]
        } else {
            let f = 1.0 / det;
            let p1_to_origin = ray.origin - self.a;
            let u = f * p1_to_origin.dot(&dir_cross_e2);
            if u < 0.0 || u > 1.0 {
                vec![]
            } else {
                let origin_cross_e1 = p1_to_origin.cross(&self.e1);
                let v = f * ray.direction.dot(&origin_cross_e1);
                if v < 0.0 || (u + v) > 1.0 {
                    vec![]
                } else {
                    let t = f * self.e2.dot(&origin_cross_e1);
                    let shape = self.tag();
                    vec![Intersection { uv: None, t: t, object: shape }]
                }
            }
        }
    }
}

pub struct SmoothTriangle {
    a: Point,
    b: Point,
    c: Point,
    n1: Vector,
    n2: Vector,
    n3: Vector,

    e1: Vector,
    e2: Vector,

    tag: Option<NodeId>,
    id: i32,
    casts_shadows: bool,
    transform: Matrix4,
    inverse_transform: Matrix4,
    pub bounds: Bounds,
    material: Option<Material>,
}

impl SmoothTriangle {
    pub fn new(a: Point, b: Point, c: Point, n1: Vector, n2: Vector, n3: Vector) -> SmoothTriangle {
        let e1 = b - a;
        let e2 = c - a;

        let bounds = Bounds(
            Point::new(
                a.x.min(b.x).min(c.x),
                a.y.min(b.y).min(c.y),
                a.z.min(b.z).min(c.z),
            ),
            Point::new(
                a.x.max(b.x).max(c.x),
                a.y.max(b.y).max(c.y),
                a.z.max(b.z).max(c.z),
            ),
        );

        SmoothTriangle {
            a, b, c, n1, n2, n3, e1, e2,
            id: 0, //gen_id(),
            casts_shadows: true,
            transform: Matrix4::id(),
            tag: None,
            bounds: bounds,
            inverse_transform: Matrix4::id().inverse(),
            material: None,
            }
    }
}

impl Shape for SmoothTriangle {
    fn id(&self) -> i32 {
        self.id
    }
    fn bounds(&self) -> &Bounds {
        &self.bounds
    }
    fn set_bounds(&mut self, bounds: Bounds) {
        self.bounds = bounds;
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
        self.material = Some(m)
    }
    fn material<'a>(&'a self, reg: &'a Registry) -> &'a Material {
        if let Some(mat) = &self.material {
            &mat
        } else {
            search_material(reg, self.tag())
        }
    }
    fn transform(&self) -> &Matrix4 {
        &self.transform
    }
    fn inverse_transform(&self) -> &Matrix4 {
        &self.inverse_transform
    }
    fn local_normal_at(&self, _p: &Point, i: &Intersection) -> Vector {
        let (u,v) = i.uv.unwrap();
        self.n2 * u + self.n3 * v + self.n1 * (1.0 - u - v)
    }
    fn local_intersect<'a>(&'a self, _reg: &'a Registry, ray: &Ray) -> Vec<Intersection> {
        let dir_cross_e2 = ray.direction.cross(&self.e2);
        let det = self.e1.dot(&dir_cross_e2);
        if det.abs() < EPSILON {
            vec![]
        } else {
            let f = 1.0 / det;
            let p1_to_origin = ray.origin - self.a;
            let u = f * p1_to_origin.dot(&dir_cross_e2);
            if u < 0.0 || u > 1.0 {
                vec![]
            } else {
                let origin_cross_e1 = p1_to_origin.cross(&self.e1);
                let v = f * ray.direction.dot(&origin_cross_e1);
                if v < 0.0 || (u + v) > 1.0 {
                    vec![]
                } else {
                    let t = f * self.e2.dot(&origin_cross_e1);
                    let shape = self.tag();
                    vec![Intersection { uv: Some((u,v)), t: t, object: shape }]
                }
            }
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::f32::consts::PI;
    use std::cell::RefCell;
    use std::rc::Rc;
    use crate::registry::id;

    struct TestShape {
        tag: Option<NodeId>,
        transform: Matrix4,
        inverse_transform: Matrix4,
        material: Option<Material>,
        bounds: Bounds,
        expected_ray: Ray,
    }

    impl TestShape {
        pub fn new(expected_ray: Ray) -> TestShape {
            TestShape {
                tag: None,
                transform: Matrix4::id(),
                material: None,
                inverse_transform: Matrix4::id().inverse(),
                bounds: Bounds(Point::origin(), Point::origin()),
                expected_ray,
            }
        }
    }

    impl Shape for TestShape {
        fn id(&self) -> i32 {
            0
        }
        fn bounds(&self) -> &Bounds {
            &self.bounds
        }
        fn set_bounds(&mut self, bounds: Bounds) {
            self.bounds = bounds;
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
            self.material = Some(m)
        }
        fn material<'a>(&'a self, _reg: &'a Registry) -> &'a Material {
            Material::base()
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

        fn local_intersect(&self, _reg: &Registry, ray: &Ray) -> Vec<Intersection> {
            if self.expected_ray.origin != ray.origin {
                panic!("Ray origin is not like expected")
            }
            if self.expected_ray.direction != ray.direction {
                panic!("Ray direction is not like expected")
            }
            vec![]
        }

        fn local_normal_at(&self, p: &Point, _i: &Intersection) -> Vector {
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
        let s = reg.get(id);
        let i = Intersection { uv: None, t: 0.0, object: id };

        assert_eq!(
            s.normal(&reg, Point::new(0.0, 1.70711, -0.70711), &i),
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
        let s = reg.get(id);
        let i = Intersection { uv: None, t: 0.0, object: id };

        assert_eq!(
            s.normal(&reg, Point::new(0.0, 2_f32.sqrt() / 2.0, -2_f32.sqrt() / 2.0), &i),
            Vector::new(0.0, 0.97014, -0.24254)
        );
    }

    mod plane {
        use super::*;

        #[test]
        fn normal_is_constant() {
            let registry = Rc::new(RefCell::new(Registry::new()));
            let id = {
                let mut reg = registry.borrow_mut();
                reg.register(Box::from(Plane::new()))
            };
            let reg = registry.borrow();
            let p = reg.get(id);
            let i = Intersection { uv: None, t: 0.0, object: id };

            assert_eq!(
                p.local_normal_at(&Point::origin(), &i),
                Vector::new(0.0, 1.0, 0.0)
            );
            assert_eq!(
                p.local_normal_at(&Point::new(10.0, 0.0, -10.0), &i),
                Vector::new(0.0, 1.0, 0.0)
            );
            assert_eq!(
                p.local_normal_at(&Point::new(-5.0, 0.0, 150.0), &i),
                Vector::new(0.0, 1.0, 0.0)
            );
        }

        #[test]
        fn intersect_with_ray_parallel_to_it() {
            let mut reg = Registry::new();
            let id = reg.register(Box::from(Plane::new()));
            let p = reg.get(id);
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
            let p = reg.get(id);
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
            let p = reg.get(id);
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
            let p = reg.get(id);
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
            let s = reg.get(id);
            let transformed = reg.get(id2);

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
            let s = reg.get(id);
            let i = Intersection { uv: None, t: 0.0, object: id };

            assert_eq!(
                s.normal(&reg, Point::new(1.0, 0.0, 0.0), &i),
                Vector::new(1.0, 0.0, 0.0)
            );
            assert_eq!(
                s.normal(&reg, Point::new(0.0, 1.0, 0.0), &i),
                Vector::new(0.0, 1.0, 0.0)
            );
            assert_eq!(
                s.normal(&reg, Point::new(0.0, 0.0, 1.0), &i),
                Vector::new(0.0, 0.0, 1.0)
            );
            let v = s.normal(&reg, Point::new(
                3_f32.sqrt() / 3.0,
                3_f32.sqrt() / 3.0,
                3_f32.sqrt() / 3.0,
            ), &i);
            assert_eq!(v, v.normalize());

            let t = reg.get(tid);
            assert_eq!(
                t.normal(&reg, Point::new(0.0, 1.70711, -0.70711), &i),
                Vector::new(0.0, 0.70711, -0.70711)
            );

            let sr = reg.get(srid);
            assert_eq!(
                sr.normal(&reg, Point::new(0.0, 2_f32.sqrt() / 2.0, -2_f32.sqrt() / 2.0), &i),
                Vector::new(0.0, 0.970140, -0.24254)
            );
        }

        #[test]
        fn intersection() {
            let mut r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
            let sph = Box::from(Sphere::new());
            let mut reg = Registry::new();
            let id = reg.register(sph);
            let s = reg.get(id);

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
            let c = reg.get(id);

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
            let registry = Rc::new(RefCell::new(Registry::new()));
            let id = {
                let mut reg = registry.borrow_mut();
                reg.register(Box::from(Cube::new()))
            };
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

            let reg = registry.borrow();
            let c = reg.get(id);
            let i = Intersection { uv: None, t: 0.0, object: id };

            for (point, normal) in &examples {
                assert_eq!(c.local_normal_at(&point, &i), *normal);
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
            let c = reg.get(id);

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
            let c = reg.get(id);

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
            let mut reg = Registry::new();
            let id = reg.register(Box::from(Cylinder::new()));
            let c = reg.get(id);

            let examples = [
                (Point::new(1.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)),
                (Point::new(0.0, 5.0, -1.0), Vector::new(0.0, 0.0, -1.0)),
                (Point::new(0.0, -2.0, 1.0), Vector::new(0.0, 0.0, 1.0)),
                (Point::new(-1.0, 1.0, 0.0), Vector::new(-1.0, 0.0, 0.0)),
            ];

            let i = Intersection { uv: None, t: 0.0, object: id };

            for (point, normal) in &examples {
                assert_eq!(c.local_normal_at(&point, &i), *normal);
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
            let c = reg.get(id);

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
            let c = reg.get(id);

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
            let mut reg = Registry::new();
            let id = reg.register(Box::from(Cylinder::closed(1.0, 2.0)));
            let c = reg.get(id);
            let examples = [
                // point, normal
                (Point::new(0.0, 1.0, 0.0), Vector::new(0.0, -1.0, 0.0)),
                (Point::new(0.5, 1.0, 0.0), Vector::new(0.0, -1.0, 0.0)),
                (Point::new(0.0, 1.0, 0.5), Vector::new(0.0, -1.0, 0.0)),
                (Point::new(0.0, 2.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
                (Point::new(0.5, 2.0, 0.0), Vector::new(0.0, 1.0, 0.0)),
                (Point::new(0.0, 2.0, 0.5), Vector::new(0.0, 1.0, 0.0)),
            ];

            let i = Intersection { uv: None, t: 0.0, object: id };

            for (point, direction) in &examples {
                assert_eq!(c.local_normal_at(&point, &i), *direction);
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
            let c = reg.get(id);

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
            let c = reg.get(id);

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
            let c = reg.get(id);

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
            let mut reg = Registry::new();
            let id = reg.register(Box::from(Cone::new()));
            let c = reg.get(id);
            let examples = [
                // point, normal
                (Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 0.0)),
                (Point::new(1.0, 1.0, 1.0), Vector::new(1.0, -(2_f32.sqrt()), 1.0)),
                (Point::new(-1.0, -1.0, 0.0), Vector::new(-1.0, 1.0, 0.0)),
            ];

            let i = Intersection { uv: None, t: 0.0, object: id };

            for (point, direction) in &examples {
                assert_eq!(c.local_normal_at(&point, &i), *direction);
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
            let g = reg.get(id);
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

            let group = reg.get(gid);

            let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
            let hits = group
                .local_intersect(&reg, &r)
                .iter()
                .map(|x| x.object)
                .collect::<Vec<NodeId>>();
            assert_eq!(hits, vec![s2id, s2id, s1id, s1id]);
        }

        #[test]
        fn intersecting_transformed_group() {
            let mut reg = Registry::new();
            let id = reg.register(Box::from(Group::new()));
            let mg = reg.get_mut(id);
            mg.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
            let sid = reg.register(Box::from(Sphere::new()));
            let ms = reg.get_mut(sid);
            ms.set_transform(Matrix4::translation(5.0, 0.0, 0.0));
            reg.add(id, sid);
            let r = Ray::new(Point::new(10.0, 0.0, -10.0), Vector::new(0.0, 0.0, 1.0));

            let g = reg.get(id);

            let hits = g
                .intersect(&reg, &r)
                .iter()
                .map(|x| x.t)
                .collect::<Vec<f32>>();
            assert_eq!(hits.len(), 2);
        }

        #[test]
        fn converting_point_from_world_to_object_space() {
            let mut reg = Registry::new();

            let id = reg.register(Box::from(Group::new()));
            let mg = reg.get_mut(id);
            mg.set_transform(Matrix4::rotation_y(PI/2.0));

            let id2 = reg.register(Box::from(Group::new()));
            let mg2 = reg.get_mut(id2);
            mg2.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));

            let sid = reg.register(Box::from(Sphere::new()));
            let ms = reg.get_mut(sid);
            ms.set_transform(Matrix4::translation(5.0, 0.0, 0.0));

            reg.add(id, id2);
            reg.add(id2, sid);

            let s = reg.get(sid);
            let p = s.world_to_object(&reg, &Point::new(-2.0, 0.0, -10.0));
            assert_eq!(p, Point::new(0.0, 0.0, -1.0));
        }

        #[test]
        fn converting_normal_vector_from_object_to_world_space() {
            let mut reg = Registry::new();

            let id = reg.register(Box::from(Group::new()));
            let mg = reg.get_mut(id);
            mg.set_transform(Matrix4::rotation_y(PI/2.0));

            let id2 = reg.register(Box::from(Group::new()));
            let mg2 = reg.get_mut(id2);
            mg2.set_transform(Matrix4::scaling(1.0, 2.0, 3.0));

            let sid = reg.register(Box::from(Sphere::new()));
            let ms = reg.get_mut(sid);
            ms.set_transform(Matrix4::translation(5.0, 0.0, 0.0));

            reg.add(id, id2);
            reg.add(id2, sid);

            let s = reg.get(sid);
            let v = s.normal_to_world(&reg, &Vector::new(3_f32.sqrt()/3.0, 3_f32.sqrt()/3.0, 3_f32.sqrt()/3.0));
            assert_eq!(v, Vector::new(0.2857, 0.4286, -0.8571));
        }

        #[test]
        fn finding_normal_on_child_object() {
            let mut reg = Registry::new();

            let id = reg.register(Box::from(Group::new()));
            let mg = reg.get_mut(id);
            mg.set_transform(Matrix4::rotation_y(PI/2.0));

            let id2 = reg.register(Box::from(Group::new()));
            let mg2 = reg.get_mut(id2);
            mg2.set_transform(Matrix4::scaling(1.0, 2.0, 3.0));

            let sid = reg.register(Box::from(Sphere::new()));
            let ms = reg.get_mut(sid);
            ms.set_transform(Matrix4::translation(5.0, 0.0, 0.0));

            reg.add(id, id2);
            reg.add(id2, sid);

            let s = reg.get(sid);
            let i = Intersection { uv: None, t: 0.0, object: sid };

            let v = s.normal(&reg, Point::new(1.7321, 1.1547, -5.5774), &i);
            assert_eq!(v, Vector::new(0.2857, 0.4286, -0.8571));
        }
    }

    mod triangle {
        use super::*;

        #[test]
        fn normal() {
            let t = Triangle::new(
                Point::new(0.0, 1.0, 0.0),
                Point::new(-1.0, 0.0, 0.0),
                Point::new(1.0, 0.0, 0.0),
            );
            let i = Intersection { uv: None, t: 0.0, object: id() };

            let n1 = t.local_normal_at(&Point::new(0.0, 0.5, 0.0), &i);
            let n2 = t.local_normal_at(&Point::new(-0.5, 0.75, 0.0), &i);
            let n3 = t.local_normal_at(&Point::new(0.5, 0.25, 0.0), &i);
            assert_eq!(n1, t.normal);
            assert_eq!(n2, t.normal);
            assert_eq!(n3, t.normal);
        }

        #[test]
        fn intersecting_parallel_ray() {
            let registry = Rc::new(RefCell::new(Registry::new()));
            let id;
            {
                let mut reg = registry.borrow_mut();
                id = reg.register(Box::from(Triangle::new(
                    Point::new(0.0, 1.0, 0.0),
                    Point::new(-1.0, 0.0, 0.0),
                    Point::new(1.0, 0.0, 0.0),
                )));
            }
            let reg = registry.borrow();
            let t = reg.get(id);
            let r = Ray::new(Point::new(0.0, -1.0, -2.0), Vector::new(0.0, 1.0, 0.0));

            let hits = t
                .intersect(&reg, &r)
                .iter()
                .map(|x| x.object)
                .collect::<Vec<NodeId>>();
            assert_eq!(hits.len(), 0);
        }

        #[test]
        fn ray_misses_p1_p3_edge() {
            let registry = Rc::new(RefCell::new(Registry::new()));
            let id;
            {
                let mut reg = registry.borrow_mut();
                id = reg.register(Box::from(Triangle::new(
                    Point::new(0.0, 1.0, 0.0),
                    Point::new(-1.0, 0.0, 0.0),
                    Point::new(1.0, 0.0, 0.0),
                )));
            }
            let reg = registry.borrow();
            let t = reg.get(id);
            let r = Ray::new(Point::new(-1.0, 1.0, -2.0), Vector::new(0.0, 0.0, 1.0));

            let hits = t
                .intersect(&reg, &r)
                .iter()
                .map(|x| x.object)
                .collect::<Vec<NodeId>>();
            assert_eq!(hits.len(), 0);
        }

        #[test]
        fn ray_misses_p2_p3_edge() {
            let registry = Rc::new(RefCell::new(Registry::new()));
            let id;
            {
                let mut reg = registry.borrow_mut();
                id = reg.register(Box::from(Triangle::new(
                    Point::new(0.0, 1.0, 0.0),
                    Point::new(-1.0, 0.0, 0.0),
                    Point::new(1.0, 0.0, 0.0),
                )));
            }
            let reg = registry.borrow();
            let t = reg.get(id);
            let r = Ray::new(Point::new(0.0, -1.0, -2.0), Vector::new(0.0, 0.0, 1.0));

            let hits = t
                .intersect(&reg, &r)
                .iter()
                .map(|x| x.object)
                .collect::<Vec<NodeId>>();
            assert_eq!(hits.len(), 0);
        }

        #[test]
        fn ray_strikes_triangle() {
            let registry = Rc::new(RefCell::new(Registry::new()));
            let id;
            {
                let mut reg = registry.borrow_mut();
                id = reg.register(Box::from(Triangle::new(
                    Point::new(0.0, 1.0, 0.0),
                    Point::new(-1.0, 0.0, 0.0),
                    Point::new(1.0, 0.0, 0.0),
                )));
            }
            let reg = registry.borrow();
            let t = reg.get(id);
            let r = Ray::new(Point::new(0.0, 0.5, -2.0), Vector::new(0.0, 0.0, 1.0));

            let hits = t
                .intersect(&reg, &r)
                .iter()
                .map(|x| x.t)
                .collect::<Vec<f32>>();
            assert_eq!(hits, vec![2.0]);
        }
    }

    mod smooth_triangle {
        use super::*;

        fn tri() -> SmoothTriangle {
            SmoothTriangle::new(
                Point::new(0.0, 1.0, 0.0),
                Point::new(-1.0, 0.0, 0.0),
                Point::new(1.0, 0.0, 0.0),
                Vector::new(0.0, 1.0, 0.0),
                Vector::new(-1.0, 0.0, 0.0),
                Vector::new(1.0, 0.0, 0.0),
            )
        }

        #[test]
        fn local_intersect() {
            let registry = Rc::new(RefCell::new(Registry::new()));
            let id;
            {
                let mut reg = registry.borrow_mut();
                id = reg.register(Box::from(tri()));
            }
            let reg = registry.borrow();
            let t = reg.get(id);
            let r = Ray::new(Point::new(-0.2, 0.3, -2.0), Vector::new(0.0, 0.0, 1.0));
            let hits = t
                .intersect(&reg, &r)
                .iter()
                .map(|x| x.uv.unwrap())
                .collect::<Vec<(f32, f32)>>();
            assert_eq!(hits, vec![(0.45, 0.25)]);
        }

        #[test]
        fn normal() {
            let registry = Rc::new(RefCell::new(Registry::new()));
            let id;
            {
                let mut reg = registry.borrow_mut();
                id = reg.register(Box::from(tri()));
            }
            let reg = registry.borrow();
            let t = reg.get(id);
            let i = Intersection {
                uv: Some((0.45, 0.25)),
                t: 1.0,
                object: id
            };
            let v = t.normal(&reg, Point::origin(), &i);
            assert_eq!(v, Vector::new(-0.5547, 0.83205, 0.0));
        }
    }
}
