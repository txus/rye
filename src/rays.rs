use crate::linear::{Matrix4, Point, Vector, EPSILON};
use crate::shapes::Shape;

pub struct Ray {
    pub origin: Point,
    pub direction: Vector,
}

impl Ray {
    pub fn new(origin: Point, direction: Vector) -> Self {
        Ray { origin, direction }
    }

    pub fn position(&self, t: f32) -> Point {
        self.origin + self.direction * t
    }

    pub fn transform(&self, transform: &Matrix4) -> Ray {
        Ray::new(*transform * self.origin, *transform * self.direction)
    }
}

#[derive(Clone, Copy)]
pub struct Intersection<'a> {
    pub t: f32,
    pub object: &'a Shape,
}

impl<'a> std::cmp::PartialEq for Intersection<'a> {
    fn eq(&self, other: &Intersection) -> bool {
        self.t == other.t && self.object.id() == other.object.id()
    }
}

pub struct Precomputation<'a> {
    pub t: f32,
    pub object: &'a Shape,
    pub point: Point,
    pub over_point: Point,
    pub under_point: Point,
    pub eye: Vector,
    pub normal: Vector,
    pub reflect: Vector,
    pub n1: f32,
    pub n2: f32,
    pub inside: bool,
}

impl<'a> Intersection<'a> {
    pub fn hit(xs: &'a [Intersection<'a>]) -> Option<Intersection> {
        let mut out: Vec<Intersection> = vec![];
        for i in xs {
            out.push(i.clone());
        }
        out.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
        let is: Vec<&Intersection> = out.iter().filter(|x| x.t >= 0.0).collect();
        match is.first() {
            Some(i) => Some(*i.clone()),
            None => None,
        }
    }

    pub fn precompute(&self, r: &Ray, is: &'a [Intersection<'a>]) -> Precomputation {
        let point = r.position(self.t);
        let normal = self.object.normal(point);
        let eye = -r.direction;
        let inside = normal.dot(&eye) < 0.0;
        let n = if inside { -normal } else { normal };

        let mut n1: f32 = 0.0;
        let mut n2: f32 = 0.0;

        let mut containers: Vec<&Shape> = vec![];
        for i in is {
            if i == self {
                n1 = if containers.is_empty() {
                    1.0
                } else {
                    containers.last().unwrap().material().refractive_index
                }
            }
            if let Some(_) = containers.iter().position(|x| x.id() == i.object.id()) {
                containers.retain(|&x| x.id() != i.object.id());
            } else {
                containers.push(i.object);
            }

            if i == self {
                n2 = if containers.is_empty() {
                    1.0
                } else {
                    containers.last().unwrap().material().refractive_index
                }
            }
        }

        Precomputation {
            t: self.t,
            object: self.object,
            point,
            over_point: point + n * EPSILON,
            under_point: point - n * EPSILON,
            eye,
            normal: n,
            reflect: r.direction.reflect(&n),
            n1,
            n2,
            inside,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shapes::{Sphere, Plane};

    #[test]
    fn initialize() {
        let origin = Point::new(1.0, 2.0, 3.0);
        let direction = Vector::new(4.0, 5.0, 6.0);
        let r = Ray::new(origin, direction);
        assert_eq!(r.origin, origin);
        assert_eq!(r.direction, direction);
    }

    #[test]
    fn position() {
        let origin = Point::new(2.0, 3.0, 4.0);
        let direction = Vector::new(1.0, 0.0, 0.0);
        let r = Ray::new(origin, direction);

        assert_eq!(r.position(0.0), origin);
        assert_eq!(r.position(1.0), Point::new(3.0, 3.0, 4.0));
        assert_eq!(r.position(-1.0), Point::new(1.0, 3.0, 4.0));
        assert_eq!(r.position(2.5), Point::new(4.5, 3.0, 4.0));
    }

    #[test]
    fn intersection_hits() {
        let s = Sphere::new();
        let mut i1 = Intersection { t: 1.0, object: &s };
        let mut i2 = Intersection { t: 2.0, object: &s };
        assert_eq!(Intersection::hit(&vec!(i1, i2)).unwrap().t, 1.0);

        i1 = Intersection {
            t: -1.0,
            object: &s,
        };
        i2 = Intersection { t: 1.0, object: &s };
        assert_eq!(Intersection::hit(&vec!(i1, i2)).unwrap().t, 1.0);

        i1 = Intersection {
            t: -2.0,
            object: &s,
        };
        i2 = Intersection {
            t: -1.0,
            object: &s,
        };
        if let Some(_) = Intersection::hit(&vec![i1, i2]) {
            assert!(false, "Something intersected when it shouldn't have")
        }

        i1 = Intersection { t: 5.0, object: &s };
        i2 = Intersection { t: 7.0, object: &s };
        let i3 = Intersection {
            t: -3.0,
            object: &s,
        };
        let i4 = Intersection { t: 2.0, object: &s };
        assert_eq!(Intersection::hit(&vec!(i1, i2, i3, i4)).unwrap().t, 2.0);
    }

    #[test]
    fn hit_when_intersection_occurs_on_outside() {
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let s = Sphere::new();
        let i = Intersection { t: 4.0, object: &s };
        let is = [i];
        let c = i.precompute(&r, &is);
        assert_eq!(c.inside, false);
    }

    #[test]
    fn hit_when_intersection_occurs_on_inside() {
        let r = Ray::new(Point::origin(), Vector::new(0.0, 0.0, 1.0));
        let s = Sphere::new();
        let i = Intersection { t: 1.0, object: &s };
        let mut is = [i];
        let c = i.precompute(&r, &mut is);
        assert_eq!(c.point, Point::new(0.0, 0.0, 1.0));
        assert_eq!(c.eye, Vector::new(0.0, 0.0, -1.0));
        assert_eq!(c.normal, Vector::new(0.0, 0.0, -1.0));
        assert_eq!(c.inside, true);
    }

    #[test]
    fn precomputing_state_intersection() {
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let s = Sphere::new();
        let i = Intersection { t: 4.0, object: &s };
        let mut is = [i];
        let c = i.precompute(&r, &mut is);
        assert_eq!(c.t, i.t);
        assert_eq!(c.point, Point::new(0.0, 0.0, -1.0));
        assert_eq!(c.eye, Vector::new(0.0, 0.0, -1.0));
        assert_eq!(c.normal, Vector::new(0.0, 0.0, -1.0));
    }

    #[test]
    fn precomputing_reflect_vector() {
        let s = Plane::new();
        let r = Ray::new(Point::new(0.0, 1.0, -1.0), Vector::new(0.0, -2_f32.sqrt()/2.0, 2_f32.sqrt()/2.0));
        let i = Intersection { t: 2_f32.sqrt(), object: &s };
        let mut is = [i];
        let c = i.precompute(&r, &mut is);
        assert_eq!(c.reflect, Vector::new(0.0, 2_f32.sqrt()/2.0, 2_f32.sqrt()/2.0));
    }

    #[test]
    fn translation() {
        let r = Ray::new(Point::new(1.0, 2.0, 3.0), Vector::new(0.0, 1.0, 0.0));
        let m = Matrix4::translation(3.0, 4.0, 5.0);
        let r2 = r.transform(&m);
        assert_eq!(r2.origin, Point::new(4.0, 6.0, 8.0));
        assert_eq!(r2.direction, Vector::new(0.0, 1.0, 0.0));
    }

    #[test]
    fn scaling() {
        let r = Ray::new(Point::new(1.0, 2.0, 3.0), Vector::new(0.0, 1.0, 0.0));
        let m = Matrix4::scaling(2.0, 3.0, 4.0);
        let r2 = r.transform(&m);
        assert_eq!(r2.origin, Point::new(2.0, 6.0, 12.0));
        assert_eq!(r2.direction, Vector::new(0.0, 3.0, 0.0));
    }

    #[test]
    fn finding_n1_and_n2_at_various_intersection() {
        let mut a = Sphere::glass();
        a.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
        a.material.refractive_index = 1.5;

        let mut b = Sphere::glass();
        b.set_transform(Matrix4::translation(0.0, 0.0, -0.25));
        b.material.refractive_index = 2.0;

        let mut c = Sphere::glass();
        c.set_transform(Matrix4::translation(0.0, 0.0, 0.25));
        c.material.refractive_index = 2.5;

        let r = Ray::new(Point::new(0.0, 0.0, -4.0), Vector::new(0.0, 0.0, 1.0));
        let is = vec![
            Intersection { t: 2.0, object: &a },
            Intersection { t: 2.75, object: &b },
            Intersection { t: 3.25, object: &c },
            Intersection { t: 4.75, object: &b },
            Intersection { t: 5.25, object: &c },
            Intersection { t: 6.0, object: &a },
        ];

        let results = vec![
            (1.0, 1.5),
            (1.5, 2.0),
            (2.0, 2.5),
            (2.5, 2.5),
            (2.5, 1.5),
            (1.5, 1.0),
        ];

        for idx in 0..6 {
            let i = is[idx];
            let comps = i.precompute(&r, &is);
            let (n1, n2) = results[idx];
            assert_eq!(comps.n1, n1);
            assert_eq!(comps.n2, n2);
        }
    }

    #[test]
    fn under_point_is_offset_below_the_surface() {
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let mut shape = Sphere::glass();
        shape.set_transform(Matrix4::translation(0.0, 0.0, 1.0));
        let i = Intersection {
            t: 5.0, object: &shape
        };
        let is = [i];
        let comps = i.precompute(&r, &is);
        assert!(comps.under_point.z > EPSILON/2.0);
        assert!(comps.point.z < comps.under_point.z);
    }
}
