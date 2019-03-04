const EPSILON: f32 = 0.00001;

#[derive(Debug, Clone, Copy)]
pub struct Point {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Vector {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

#[derive(Debug, Clone, Copy)]
pub struct Color {
    pub red: f32,
    pub green: f32,
    pub blue: f32,
}

impl Point {
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self { x: x, y: y, z: z }
    }
}

impl Color {
    pub fn new(red: f32, green: f32, blue: f32) -> Self {
        Self {
            red: red,
            green: green,
            blue: blue,
        }
    }
}

impl Vector {
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self { x: x, y: y, z: z }
    }
    pub fn magnitude(&self) -> f32 {
        (self.x.powf(2.0) + self.y.powf(2.0) + self.z.powf(2.0)).sqrt()
    }
    pub fn normalize(&self) -> Self {
        let m = self.magnitude();
        Self::new(self.x / m, self.y / m, self.z / m)
    }
    pub fn dot(&self, other: &Self) -> f32 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
    pub fn cross(&self, other: &Self) -> Self {
        Self::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }
}

fn close(a: f32, b: f32) -> bool {
    (a - b).abs() < EPSILON
}

impl std::cmp::PartialEq for Point {
    fn eq(&self, other: &Self) -> bool {
        close(self.x, other.x) && close(self.y, other.y) && close(self.z, other.z)
    }
}

impl std::cmp::PartialEq for Color {
    fn eq(&self, other: &Self) -> bool {
        close(self.red, other.red) && close(self.green, other.green) && close(self.blue, other.blue)
    }
}

impl std::ops::Add<Vector> for Point {
    type Output = Self;
    fn add(self, other: Vector) -> Self {
        Point::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

impl std::ops::Add<Vector> for Vector {
    type Output = Self;
    fn add(self, other: Vector) -> Self {
        Self::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }
}

impl std::ops::Sub<Vector> for Point {
    type Output = Self;
    fn sub(self, other: Vector) -> Self {
        Self::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl std::ops::Sub<Vector> for Vector {
    type Output = Self;
    fn sub(self, other: Vector) -> Self {
        Self::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }
}

impl std::ops::Neg for Vector {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

impl std::ops::Mul<f32> for Vector {
    type Output = Self;
    fn mul(self, scalar: f32) -> Self {
        Self::new(self.x * scalar, self.y * scalar, self.z * scalar)
    }
}

impl std::ops::Div<f32> for Vector {
    type Output = Self;
    fn div(self, scalar: f32) -> Self {
        Self::new(self.x / scalar, self.y / scalar, self.z / scalar)
    }
}

impl std::ops::Add<Color> for Color {
    type Output = Self;
    fn add(self, other: Color) -> Self {
        Self::new(
            self.red + other.red,
            self.green + other.green,
            self.blue + other.blue,
        )
    }
}

impl std::ops::Sub<Color> for Color {
    type Output = Self;
    fn sub(self, other: Color) -> Self {
        Self::new(
            self.red - other.red,
            self.green - other.green,
            self.blue - other.blue,
        )
    }
}

impl std::ops::Mul<f32> for Color {
    type Output = Self;
    fn mul(self, scalar: f32) -> Self {
        Self::new(self.red * scalar, self.green * scalar, self.blue * scalar)
    }
}

impl std::ops::Mul<Color> for Color {
    type Output = Self;
    fn mul(self, other: Color) -> Self {
        Self::new(
            self.red * other.red,
            self.green * other.green,
            self.blue * other.blue,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn points() {
        let p = Point::new(4.3, -4.2, 3.1);
        assert_eq!(p.x, 4.3);
        assert_eq!(p.y, -4.2);
        assert_eq!(p.z, 3.1);
    }

    #[test]
    fn vectors() {
        let v = Vector::new(4.3, -4.2, 3.1);
        assert_eq!(v.x, 4.3);
        assert_eq!(v.y, -4.2);
        assert_eq!(v.z, 3.1);
    }

    #[test]
    fn points_equality() {
        let p1 = Point::new(4.3, -4.2, 3.1);
        let p2 = Point::new(4.3 + (EPSILON / 2.0), -4.2, 3.1);
        let p3 = Point::new(4.3 + EPSILON, -4.2, 3.1);

        assert!(p1 == p2);
        assert!(p1 != p3);
        assert!(p2 == p3);
    }

    #[test]
    fn addition() {
        let p = Point::new(3.0, -2.0, 5.0);
        let v = Vector::new(-2.0, 3.0, 1.0);

        assert_eq!(p + v, Point::new(1.0, 1.0, 6.0));
    }

    #[test]
    fn subtraction() {
        let p = Point::new(3.0, 2.0, 1.0);
        let v = Vector::new(5.0, 6.0, 7.0);

        assert_eq!(p - v, Point::new(-2.0, -4.0, -6.0));
    }

    #[test]
    fn vector_subtraction() {
        let v1 = Vector::new(3.0, 2.0, 1.0);
        let v2 = Vector::new(5.0, 6.0, 7.0);

        assert_eq!(v1 - v2, Vector::new(-2.0, -4.0, -6.0));
    }

    #[test]
    fn vector_negation() {
        let v = Vector::new(1.0, -2.0, 3.0);

        assert_eq!(-v, Vector::new(-1.0, 2.0, -3.0));
    }

    #[test]
    fn vector_scaling() {
        let v = Vector::new(1.0, -2.0, 3.0);

        assert_eq!(v * 3.5, Vector::new(3.5, -7.0, 10.5));
        assert_eq!(v / 2.0, Vector::new(0.5, -1.0, 1.5));
    }

    #[test]
    fn vector_magnitude() {
        assert_eq!(Vector::new(1.0, 0.0, 0.0).magnitude(), 1.0);
        assert_eq!(Vector::new(0.0, 1.0, 0.0).magnitude(), 1.0);
        assert_eq!(Vector::new(0.0, 0.0, 1.0).magnitude(), 1.0);
        assert_eq!(Vector::new(1.0, 2.0, 3.0).magnitude(), 14_f32.sqrt());
        assert_eq!(Vector::new(-1.0, -2.0, -3.0).magnitude(), 14_f32.sqrt());
    }

    #[test]
    fn vector_normalization() {
        assert_eq!(
            Vector::new(4.0, 0.0, 0.0).normalize(),
            Vector::new(1.0, 0.0, 0.0)
        );
        assert_eq!(
            Vector::new(1.0, 2.0, 3.0).normalize(),
            Vector::new(
                1.0 / 14_f32.sqrt(),
                2.0 / 14_f32.sqrt(),
                3.0 / 14_f32.sqrt()
            )
        );
        assert!(
            (Vector::new(1.0, 2.0, 3.0).normalize().magnitude() - 1.0).abs() < EPSILON,
            "A normalized vector's magnitude is 1"
        );
    }

    #[test]
    fn vector_dot_product() {
        let v1 = Vector::new(1.0, 2.0, 3.0);
        let v2 = Vector::new(2.0, 3.0, 4.0);
        assert_eq!(v1.dot(&v2), 20.0);
    }

    #[test]
    fn vector_cross_product() {
        let v1 = Vector::new(1.0, 2.0, 3.0);
        let v2 = Vector::new(2.0, 3.0, 4.0);
        assert_eq!(v1.cross(&v2), Vector::new(-1.0, 2.0, -1.0));
    }

    #[test]
    fn colors() {
        let c = Color::new(-0.5, 0.4, 1.7);
        assert_eq!(c.red, -0.5);
        assert_eq!(c.green, 0.4);
        assert_eq!(c.blue, 1.7);
    }

    #[test]
    fn color_addition() {
        let v1 = Color::new(0.9, 0.6, 0.75);
        let v2 = Color::new(0.7, 0.1, 0.25);

        assert_eq!(v1 + v2, Color::new(1.6, 0.7, 1.0));
    }

    #[test]
    fn color_subtraction() {
        let v1 = Color::new(0.9, 0.6, 0.75);
        let v2 = Color::new(0.7, 0.1, 0.25);

        assert_eq!(v1 - v2, Color::new(0.2, 0.5, 0.5));
    }

    #[test]
    fn color_scaling() {
        let v = Color::new(0.2, 0.3, 0.4);

        assert_eq!(v * 2.0, Color::new(0.4, 0.6, 0.8));
    }

    #[test]
    fn color_hadamard_product() {
        let v1 = Color::new(1.0, 0.2, 0.4);
        let v2 = Color::new(0.9, 1.0, 0.1);

        assert_eq!(v1 * v2, Color::new(0.9, 0.2, 0.04));
    }
}
