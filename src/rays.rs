use crate::primitives::{Matrix4, Point, Vector, EPSILON};
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

    pub fn transform(&self, transform: Matrix4) -> Ray {
        Ray::new(transform * self.origin, transform * self.direction)
    }
}

pub struct Intersection<'a> {
    pub t: f32,
    pub object: &'a Shape,
}

pub struct Precomputation<'a> {
    pub t: f32,
    pub object: &'a Shape,
    pub point: Point,
    pub over_point: Point,
    pub eye: Vector,
    pub normal: Vector,
    pub inside: bool,
}

impl<'a> Intersection<'a> {
    pub fn hit(xs: &'a mut Vec<Intersection<'a>>) -> Option<&'a Intersection<'a>> {
        xs.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
        let is: Vec<&Intersection> = xs.iter().filter(|x| x.t >= 0.0).collect();
        match is.first() {
            Some(i) => Some(*i),
            None => None,
        }
    }

    pub fn precompute(&self, r: &Ray) -> Precomputation {
        let point = r.position(self.t);
        let normal = self.object.normal(point);
        let eye = -r.direction;
        let inside = normal.dot(&eye) < 0.0;
        let n = if inside { -normal } else { normal };
        Precomputation {
            t: self.t,
            object: self.object,
            point,
            over_point: point + n * EPSILON,
            eye,
            normal: n,
            inside,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::Matrix;
    use crate::shapes::Sphere;
    use std::f32::consts::PI;

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
    fn sphere_intersection() {
        let mut r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let s = Sphere::unit();
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

        r = Ray::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0));
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

    #[test]
    fn intersection_hits() {
        let s = Sphere::unit();
        let mut i1 = Intersection { t: 1.0, object: &s };
        let mut i2 = Intersection { t: 2.0, object: &s };
        assert_eq!(Intersection::hit(&mut vec!(i1, i2)).unwrap().t, 1.0);

        i1 = Intersection {
            t: -1.0,
            object: &s,
        };
        i2 = Intersection { t: 1.0, object: &s };
        assert_eq!(Intersection::hit(&mut vec!(i1, i2)).unwrap().t, 1.0);

        i1 = Intersection {
            t: -2.0,
            object: &s,
        };
        i2 = Intersection {
            t: -1.0,
            object: &s,
        };
        // dont' know how to do it
        //let x: Option<Intersection> = None;
        //assert_eq!(Intersection::hit(&mut vec!(i1, i2)), x);

        i1 = Intersection { t: 5.0, object: &s };
        i2 = Intersection { t: 7.0, object: &s };
        let i3 = Intersection {
            t: -3.0,
            object: &s,
        };
        let i4 = Intersection { t: 2.0, object: &s };
        assert_eq!(Intersection::hit(&mut vec!(i1, i2, i3, i4)).unwrap().t, 2.0);
    }

    #[test]
    fn hit_when_intersection_occurs_on_outside() {
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let s = Sphere::unit();
        let i = Intersection { t: 4.0, object: &s };
        let c = i.precompute(&r);
        assert_eq!(c.inside, false);
    }

    #[test]
    fn hit_when_intersection_occurs_on_inside() {
        let r = Ray::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0));
        let s = Sphere::unit();
        let i = Intersection { t: 1.0, object: &s };
        let c = i.precompute(&r);
        assert_eq!(c.point, Point::new(0.0, 0.0, 1.0));
        assert_eq!(c.eye, Vector::new(0.0, 0.0, -1.0));
        assert_eq!(c.normal, Vector::new(0.0, 0.0, -1.0));
        assert_eq!(c.inside, true);
    }

    #[test]
    fn precomputing_state_intersection() {
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let s = Sphere::unit();
        let i = Intersection { t: 4.0, object: &s };
        let c = i.precompute(&r);
        assert_eq!(c.t, i.t);
        assert_eq!(c.point, Point::new(0.0, 0.0, -1.0));
        assert_eq!(c.eye, Vector::new(0.0, 0.0, -1.0));
        assert_eq!(c.normal, Vector::new(0.0, 0.0, -1.0));
    }

    #[test]
    fn translation() {
        let r = Ray::new(Point::new(1.0, 2.0, 3.0), Vector::new(0.0, 1.0, 0.0));
        let m = Matrix4::translation(3.0, 4.0, 5.0);
        let r2 = r.transform(m);
        assert_eq!(r2.origin, Point::new(4.0, 6.0, 8.0));
        assert_eq!(r2.direction, Vector::new(0.0, 1.0, 0.0));
    }

    #[test]
    fn scaling() {
        let r = Ray::new(Point::new(1.0, 2.0, 3.0), Vector::new(0.0, 1.0, 0.0));
        let m = Matrix4::scaling(2.0, 3.0, 4.0);
        let r2 = r.transform(m);
        assert_eq!(r2.origin, Point::new(2.0, 6.0, 12.0));
        assert_eq!(r2.direction, Vector::new(0.0, 3.0, 0.0));
    }

    #[test]
    fn sphere_transform() {
        let mut s = Sphere::unit();
        assert_eq!(s.transform(), Matrix4::id());
        let t = Matrix4::translation(2.0, 3.0, 4.0);
        s.set_transform(t);
        assert_eq!(s.transform(), t);
    }

    #[test]
    fn scaled_sphere_intersection() {
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let mut s = Sphere::unit();
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
    fn sphere_normal() {
        let mut s = Sphere::unit();
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
}
