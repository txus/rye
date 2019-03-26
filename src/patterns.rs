use crate::color::Color;
use crate::linear::{Matrix, Matrix4, Point};
use crate::shapes::Shape;

pub trait Pattern: Send + Sync {
    fn transform(&self) -> &Matrix4;
    fn set_transform(&mut self, t: Matrix4);
    fn color_at(&self, p: &Point) -> Color;

    fn color_at_object(&self, s: &Shape, p: &Point) -> Color {
        let object_point = s.transform().inverse() * *p;
        let pattern_point = self.transform().inverse() * object_point;
        self.color_at(&pattern_point)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StripePattern {
    pub a: Color,
    pub b: Color,
    pub transform: Matrix4,
}

impl StripePattern {
    pub fn new(a: Color, b: Color) -> StripePattern {
        StripePattern {
            a,
            b,
            transform: Matrix4::id(),
        }
    }
}

impl Pattern for StripePattern {
    fn transform(&self) -> &Matrix4 {
        &self.transform
    }

    fn set_transform(&mut self, t: Matrix4) {
        self.transform = t;
    }

    fn color_at(&self, p: &Point) -> Color {
        if (p.x.floor() as u32 % 2) == 0 {
            self.a
        } else {
            self.b
        }
    }
}

pub struct TestPattern {
    pub transform: Matrix4,
}

impl TestPattern {
    pub fn new() -> TestPattern {
        TestPattern {
            transform: Matrix4::id(),
        }
    }
}

impl Pattern for TestPattern {
    fn transform(&self) -> &Matrix4 {
        &self.transform
    }

    fn set_transform(&mut self, t: Matrix4) {
        self.transform = t;
    }

    fn color_at(&self, p: &Point) -> Color {
        Color::new(p.x, p.y, p.z)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GradientPattern {
    pub a: Color,
    pub b: Color,
    pub transform: Matrix4,
}

impl GradientPattern {
    pub fn new(a: Color, b: Color) -> GradientPattern {
        GradientPattern {
            a,
            b,
            transform: Matrix4::id(),
        }
    }
}

impl Pattern for GradientPattern {
    fn transform(&self) -> &Matrix4 {
        &self.transform
    }

    fn set_transform(&mut self, t: Matrix4) {
        self.transform = t;
    }

    fn color_at(&self, p: &Point) -> Color {
        let distance = self.b - self.a;
        let fraction = p.x - p.x.floor();
        self.a + distance * fraction
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RingPattern {
    pub a: Color,
    pub b: Color,
    pub transform: Matrix4,
}

impl RingPattern {
    pub fn new(a: Color, b: Color) -> RingPattern {
        RingPattern {
            a,
            b,
            transform: Matrix4::id(),
        }
    }
}

impl Pattern for RingPattern {
    fn transform(&self) -> &Matrix4 {
        &self.transform
    }

    fn set_transform(&mut self, t: Matrix4) {
        self.transform = t;
    }

    fn color_at(&self, p: &Point) -> Color {
        let px2 = p.x.powi(2);
        let pz2 = p.z.powi(2);
        if (px2 + pz2).sqrt().floor() % 2.0 == 0.0 {
            self.a
        } else {
            self.b
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CheckerPattern {
    pub a: Color,
    pub b: Color,
    pub transform: Matrix4,
}

impl CheckerPattern {
    pub fn new(a: Color, b: Color) -> CheckerPattern {
        CheckerPattern {
            a,
            b,
            transform: Matrix4::id(),
        }
    }
}

impl Pattern for CheckerPattern {
    fn transform(&self) -> &Matrix4 {
        &self.transform
    }

    fn set_transform(&mut self, t: Matrix4) {
        self.transform = t;
    }

    fn color_at(&self, p: &Point) -> Color {
        if (p.x.floor() + p.y.floor() + p.z.floor()) % 2.0 == 0.0 {
            self.a
        } else {
            self.b
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    mod stripe {
        use super::*;
        use crate::linear::Matrix4;
        use crate::shapes::Sphere;

        #[test]
        fn initialize() {
            let p = StripePattern::new(Color::white(), Color::black());
            assert_eq!(p.a, Color::white());
            assert_eq!(p.b, Color::black());
        }

        #[test]
        fn constant_in_y() {
            let p = StripePattern::new(Color::white(), Color::black());
            assert_eq!(p.color_at(&Point::origin()), Color::white());
            assert_eq!(p.color_at(&Point::new(0.0, 1.0, 0.0)), Color::white());
            assert_eq!(p.color_at(&Point::new(0.0, 2.0, 0.0)), Color::white());
        }

        #[test]
        fn constant_in_z() {
            let p = StripePattern::new(Color::white(), Color::black());
            assert_eq!(p.color_at(&Point::origin()), Color::white());
            assert_eq!(p.color_at(&Point::new(0.0, 0.0, 1.0)), Color::white());
            assert_eq!(p.color_at(&Point::new(0.0, 0.0, 2.0)), Color::white());
        }

        #[test]
        fn alternates_in_x() {
            let p = StripePattern::new(Color::white(), Color::black());
            assert_eq!(p.color_at(&Point::origin()), Color::white());
            assert_eq!(p.color_at(&Point::new(0.9, 0.0, 0.0)), Color::white());
            assert_eq!(p.color_at(&Point::new(1.0, 0.0, 0.0)), Color::black());
            assert_eq!(p.color_at(&Point::new(-0.1, 0.0, 0.0)), Color::black());
            assert_eq!(p.color_at(&Point::new(-1.0, 0.0, 0.0)), Color::black());
            assert_eq!(p.color_at(&Point::new(-1.1, 0.0, 0.0)), Color::white());
        }

        #[test]
        fn stripe_with_object_transformation() {
            let mut s = Sphere::new();
            s.set_transform(Matrix4::scaling(2.2, 2.2, 2.2));
            let p = StripePattern::new(Color::white(), Color::black());
            let c = p.color_at_object(&s, &Point::new(1.5, 0.0, 0.0));
            assert_eq!(c, Color::white());
        }

        #[test]
        fn stripe_with_pattern_transformation() {
            let s = Sphere::new();
            let mut p = StripePattern::new(Color::white(), Color::black());
            p.set_transform(Matrix4::scaling(2.0, 2.0, 2.0));
            let c = p.color_at_object(&s, &Point::new(1.5, 0.0, 0.0));
            assert_eq!(c, Color::white());
        }

        #[test]
        fn stripe_with_both_object_and_pattern_transformation() {
            let mut s = Sphere::new();
            s.set_transform(Matrix4::scaling(2.2, 2.2, 2.2));
            let mut p = StripePattern::new(Color::white(), Color::black());
            p.set_transform(Matrix4::translation(0.5, 0.0, 0.0));
            let c = p.color_at_object(&s, &Point::new(2.5, 0.0, 0.0));
            assert_eq!(c, Color::white());
        }
    }

    mod gradient {
        use super::*;

        #[test]
        fn linear_interpolation() {
            let p = GradientPattern::new(Color::white(), Color::black());
            assert_eq!(p.color_at(&Point::origin()), Color::white());
            assert_eq!(
                p.color_at(&Point::new(0.25, 0.0, 0.0)),
                Color::new(0.75, 0.75, 0.75)
            );
            assert_eq!(
                p.color_at(&Point::new(0.5, 0.0, 0.0)),
                Color::new(0.5, 0.5, 0.5)
            );
            assert_eq!(
                p.color_at(&Point::new(0.75, 0.0, 0.0)),
                Color::new(0.25, 0.25, 0.25)
            );
        }
    }

    mod ring {
        use super::*;

        #[test]
        fn extends_in_both_x_and_z() {
            let p = RingPattern::new(Color::white(), Color::black());
            assert_eq!(p.color_at(&Point::origin()), Color::white());
            assert_eq!(p.color_at(&Point::new(1.0, 0.0, 0.0)), Color::black());
            assert_eq!(p.color_at(&Point::new(0.0, 0.0, 1.0)), Color::black());
            assert_eq!(p.color_at(&Point::new(0.708, 0.0, 0.708)), Color::black());
        }
    }

    mod checker {
        use super::*;

        #[test]
        fn repeat_in_x() {
            let p = CheckerPattern::new(Color::white(), Color::black());
            assert_eq!(p.color_at(&Point::origin()), Color::white());
            assert_eq!(p.color_at(&Point::new(0.99, 0.0, 0.0)), Color::white());
            assert_eq!(p.color_at(&Point::new(1.01, 0.0, 0.0)), Color::black());
        }

        #[test]
        fn repeat_in_y() {
            let p = CheckerPattern::new(Color::white(), Color::black());
            assert_eq!(p.color_at(&Point::origin()), Color::white());
            assert_eq!(p.color_at(&Point::new(0.0, 0.99, 0.0)), Color::white());
            assert_eq!(p.color_at(&Point::new(0.0, 1.01, 0.0)), Color::black());
        }

        #[test]
        fn repeat_in_z() {
            let p = CheckerPattern::new(Color::white(), Color::black());
            assert_eq!(p.color_at(&Point::origin()), Color::white());
            assert_eq!(p.color_at(&Point::new(0.0, 0.0, 0.99)), Color::white());
            assert_eq!(p.color_at(&Point::new(0.0, 0.0, 1.01)), Color::black());
        }
    }
}
