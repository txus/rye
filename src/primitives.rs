pub const EPSILON: f32 = 0.0001;

#[derive(Debug, Clone, Copy)]
pub struct Point {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

#[derive(Debug, Clone, Copy)]
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

pub trait Matrix {
    fn at(&self, row: usize, column: usize) -> f32;
    fn id() -> Self;
}

#[derive(Debug, Clone, Copy)]
pub struct Matrix4 {
    m: [[f32; 4]; 4],
}

impl Matrix4 {
    pub fn new(m: [[f32; 4]; 4]) -> Self {
        Self { m }
    }

    pub fn translation(x: f32, y: f32, z: f32) -> Self {
        Self::new([
            [1.0, 0.0, 0.0, x],
            [0.0, 1.0, 0.0, y],
            [0.0, 0.0, 1.0, z],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }

    pub fn scaling(x: f32, y: f32, z: f32) -> Self {
        Self::new([
            [x, 0.0, 0.0, 0.0],
            [0.0, y, 0.0, 0.0],
            [0.0, 0.0, z, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }

    pub fn rotation_x(radians: f32) -> Self {
        Self::new([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, radians.cos(), -(radians.sin()), 0.0],
            [0.0, radians.sin(), radians.cos(), 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }

    pub fn rotation_y(radians: f32) -> Self {
        Self::new([
            [radians.cos(), 0.0, radians.sin(), 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [-(radians.sin()), 0.0, radians.cos(), 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }

    pub fn rotation_z(radians: f32) -> Self {
        Self::new([
            [radians.cos(), -(radians.sin()), 0.0, 0.0],
            [radians.sin(), radians.cos(), 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }

    pub fn shearing(xy: f32, xz: f32, yx: f32, yz: f32, zx: f32, zy: f32) -> Self {
        Self::new([
            [1.0, xy, xz, 0.0],
            [yx, 1.0, yz, 0.0],
            [zx, zy, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }

    pub fn submatrix(&self, row: usize, column: usize) -> Matrix3 {
        let v: Vec<Vec<f32>> = self
            .m
            .iter()
            .enumerate()
            .filter(|&(i, _)| i != row)
            .map(|(_, cs)| {
                cs.iter()
                    .enumerate()
                    .filter(|&(i, _)| i != column)
                    .map(|(_, v)| *v)
                    .collect::<Vec<f32>>()
            })
            .collect();

        Matrix3::new([
            [v[0][0], v[0][1], v[0][2]],
            [v[1][0], v[1][1], v[1][2]],
            [v[2][0], v[2][1], v[2][2]],
        ])
    }

    pub fn transpose(&self) -> Self {
        Self::new([
            [self.at(0, 0), self.at(1, 0), self.at(2, 0), self.at(3, 0)],
            [self.at(0, 1), self.at(1, 1), self.at(2, 1), self.at(3, 1)],
            [self.at(0, 2), self.at(1, 2), self.at(2, 2), self.at(3, 2)],
            [self.at(0, 3), self.at(1, 3), self.at(2, 3), self.at(3, 3)],
        ])
    }

    pub fn minor(&self, row: usize, column: usize) -> f32 {
        self.submatrix(row, column).determinant()
    }

    pub fn cofactor(&self, row: usize, column: usize) -> f32 {
        let m = self.minor(row, column);
        if (row + column) % 2 == 0 {
            m
        } else {
            -m
        }
    }

    pub fn determinant(&self) -> f32 {
        let row = self.m[0];
        row[0] * self.cofactor(0, 0)
            + row[1] * self.cofactor(0, 1)
            + row[2] * self.cofactor(0, 2)
            + row[3] * self.cofactor(0, 3)
    }

    pub fn inverse(&self) -> Self {
        let d = self.determinant();
        if d == 0.0 {
            panic!("Can't invert matrix");
        }

        Self::new([
            [
                self.cofactor(0, 0) / d,
                self.cofactor(1, 0) / d,
                self.cofactor(2, 0) / d,
                self.cofactor(3, 0) / d,
            ],
            [
                self.cofactor(0, 1) / d,
                self.cofactor(1, 1) / d,
                self.cofactor(2, 1) / d,
                self.cofactor(3, 1) / d,
            ],
            [
                self.cofactor(0, 2) / d,
                self.cofactor(1, 2) / d,
                self.cofactor(2, 2) / d,
                self.cofactor(3, 2) / d,
            ],
            [
                self.cofactor(0, 3) / d,
                self.cofactor(1, 3) / d,
                self.cofactor(2, 3) / d,
                self.cofactor(3, 3) / d,
            ],
        ])
    }
}

impl Matrix for Matrix4 {
    fn at(&self, row: usize, column: usize) -> f32 {
        self.m[row][column]
    }
    fn id() -> Self {
        Self::new([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ])
    }
}

impl std::cmp::PartialEq for Matrix4 {
    fn eq(&self, other: &Self) -> bool {
        for row in 0..4 {
            for column in 0..4 {
                if !close(self.at(row, column), other.at(row, column)) {
                    return false;
                }
            }
        }
        true
    }
}

impl std::cmp::PartialEq for Matrix3 {
    fn eq(&self, other: &Self) -> bool {
        for row in 0..3 {
            for column in 0..3 {
                if !close(self.at(row, column), other.at(row, column)) {
                    return false;
                }
            }
        }
        true
    }
}

impl std::cmp::PartialEq for Matrix2 {
    fn eq(&self, other: &Self) -> bool {
        for row in 0..2 {
            for column in 0..2 {
                if !close(self.at(row, column), other.at(row, column)) {
                    return false;
                }
            }
        }
        true
    }
}

fn c(a: &Matrix4, b: &Matrix4, row: usize, column: usize) -> f32 {
    a.at(row, 0) * b.at(0, column)
        + a.at(row, 1) * b.at(1, column)
        + a.at(row, 2) * b.at(2, column)
        + a.at(row, 3) * b.at(3, column)
}

impl std::ops::Mul<Matrix4> for Matrix4 {
    type Output = Self;
    fn mul(self, other: Matrix4) -> Self {
        Self::new([
            [
                c(&self, &other, 0, 0),
                c(&self, &other, 0, 1),
                c(&self, &other, 0, 2),
                c(&self, &other, 0, 3),
            ],
            [
                c(&self, &other, 1, 0),
                c(&self, &other, 1, 1),
                c(&self, &other, 1, 2),
                c(&self, &other, 1, 3),
            ],
            [
                c(&self, &other, 2, 0),
                c(&self, &other, 2, 1),
                c(&self, &other, 2, 2),
                c(&self, &other, 2, 3),
            ],
            [
                c(&self, &other, 3, 0),
                c(&self, &other, 3, 1),
                c(&self, &other, 3, 2),
                c(&self, &other, 3, 3),
            ],
        ])
    }
}

impl std::ops::Mul<(f32, f32, f32, f32)> for Matrix4 {
    type Output = (f32, f32, f32, f32);
    fn mul(self, other: (f32, f32, f32, f32)) -> (f32, f32, f32, f32) {
        let (a, b, c, d) = other;
        (
            self.at(0, 0) * a + self.at(0, 1) * b + self.at(0, 2) * c + self.at(0, 3) * d,
            self.at(1, 0) * a + self.at(1, 1) * b + self.at(1, 2) * c + self.at(1, 3) * d,
            self.at(2, 0) * a + self.at(2, 1) * b + self.at(2, 2) * c + self.at(2, 3) * d,
            self.at(3, 0) * a + self.at(3, 1) * b + self.at(3, 2) * c + self.at(3, 3) * d,
        )
    }
}

impl std::ops::Mul<Point> for Matrix4 {
    type Output = Point;
    fn mul(self, Point { x, y, z }: Point) -> Point {
        let (x2, y2, z2, _) = self * (x, y, z, 1.0);
        Point::new(x2, y2, z2)
    }
}

impl std::ops::Mul<Vector> for Matrix4 {
    type Output = Vector;
    fn mul(self, Vector { x, y, z }: Vector) -> Vector {
        let (x2, y2, z2, _) = self * (x, y, z, 0.0);
        Vector::new(x2, y2, z2)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Matrix3 {
    m: [[f32; 3]; 3],
}

impl Matrix3 {
    pub fn new(m: [[f32; 3]; 3]) -> Self {
        Self { m }
    }

    pub fn submatrix(&self, row: usize, column: usize) -> Matrix2 {
        let v: Vec<Vec<f32>> = self
            .m
            .iter()
            .enumerate()
            .filter(|&(i, _)| i != row)
            .map(|(_, cs)| {
                cs.iter()
                    .enumerate()
                    .filter(|&(i, _)| i != column)
                    .map(|(_, v)| *v)
                    .collect::<Vec<f32>>()
            })
            .collect();

        Matrix2::new([[v[0][0], v[0][1]], [v[1][0], v[1][1]]])
    }

    pub fn minor(&self, row: usize, column: usize) -> f32 {
        self.submatrix(row, column).determinant()
    }

    pub fn cofactor(&self, row: usize, column: usize) -> f32 {
        let m = self.minor(row, column);
        if (row + column) % 2 == 0 {
            m
        } else {
            -m
        }
    }

    pub fn determinant(&self) -> f32 {
        let row = self.m[0];
        row[0] * self.cofactor(0, 0) + row[1] * self.cofactor(0, 1) + row[2] * self.cofactor(0, 2)
    }
}

impl Matrix for Matrix3 {
    fn at(&self, row: usize, column: usize) -> f32 {
        self.m[row][column]
    }
    fn id() -> Self {
        Self::new([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Matrix2 {
    m: [[f32; 2]; 2],
}

impl Matrix2 {
    pub fn new(m: [[f32; 2]; 2]) -> Self {
        Self { m }
    }

    pub fn determinant(&self) -> f32 {
        self.at(0, 0) * self.at(1, 1) - self.at(0, 1) * self.at(1, 0)
    }
}

impl Matrix for Matrix2 {
    fn at(&self, row: usize, column: usize) -> f32 {
        self.m[row][column]
    }
    fn id() -> Self {
        Self::new([[1.0, 0.0], [0.0, 1.0]])
    }
}

impl Point {
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self { x, y, z }
    }
}

impl Color {
    pub fn new(red: f32, green: f32, blue: f32) -> Self {
        Self { red, green, blue }
    }
    pub fn black() -> Color {
        Self::new(0.0, 0.0, 0.0)
    }

    pub fn white() -> Color {
        Self::new(1.0, 1.0, 1.0)
    }

    pub fn red() -> Color {
        Self::new(1.0, 0.0, 0.0)
    }

    pub fn green() -> Color {
        Self::new(0.0, 1.0, 0.0)
    }

    pub fn blue() -> Color {
        Self::new(0.0, 0.0, 1.0)
    }
}

impl Vector {
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self { x, y, z }
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

    pub fn reflect(&self, normal: &Vector) -> Vector {
        *self - *normal * 2.0 * self.dot(normal)
    }
}

fn close(a: f32, b: f32) -> bool {
    (a - b).abs() < (EPSILON * 2.0)
}

impl std::cmp::PartialEq for Point {
    fn eq(&self, other: &Self) -> bool {
        close(self.x, other.x) && close(self.y, other.y) && close(self.z, other.z)
    }
}

impl std::cmp::PartialEq for Vector {
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

impl std::ops::Sub<Point> for Point {
    type Output = Vector;
    fn sub(self, other: Point) -> Vector {
        Vector::new(self.x - other.x, self.y - other.y, self.z - other.z)
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
    use std::f32::consts::PI;

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
        let p2 = Point {
            x: 4.3 + (2.0 * EPSILON),
            ..p1
        };
        let p3 = Point {
            x: 4.3 + (2.01 * EPSILON),
            ..p1
        };

        assert!(p1 == p2);
        assert!(p1 != p3);
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

    #[test]
    fn matrix_init() {
        let m = Matrix4::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.5, 6.5, 7.5, 8.5],
            [9.0, 10.0, 11.0, 12.0],
            [13.5, 14.5, 15.5, 16.5],
        ]);

        assert_eq!(m.at(0, 0), 1.0);
        assert_eq!(m.at(0, 3), 4.0);
        assert_eq!(m.at(1, 0), 5.5);
        assert_eq!(m.at(1, 2), 7.5);
        assert_eq!(m.at(2, 2), 11.0);
        assert_eq!(m.at(3, 0), 13.5);
        assert_eq!(m.at(3, 2), 15.5);
    }

    #[test]
    fn matrix_equality() {
        let m = Matrix4::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.5, 6.5, 7.5, 8.5],
            [9.0, 10.0, 11.0, 12.0],
            [13.5, 14.5, 15.5, 16.5],
        ]);
        let n = Matrix4::new([
            [1.0 + EPSILON / 2.0, 2.0, 3.0, 4.0],
            [5.5, 6.5, 7.5, 8.5],
            [9.0, 10.0, 11.0 - EPSILON / 2.0, 12.0],
            [13.5, 14.5, 15.5, 16.5],
        ]);
        let other = Matrix4::new([
            [1.0, 2.0, 3.0, 4.0],
            [9.0, 10.0, 11.0, 12.0],
            [5.5, 6.5, 7.5, 8.5],
            [13.5, 14.5, 15.5, 16.5],
        ]);

        assert_eq!(m, m);
        assert_eq!(m, n);
        assert!(m != other);
    }
    #[test]
    fn matrix_multiplication() {
        let a = Matrix4::new([
            [1.0, 2.0, 3.0, 4.0],
            [5.0, 6.0, 7.0, 8.0],
            [9.0, 8.0, 7.0, 6.0],
            [5.0, 4.0, 3.0, 2.0],
        ]);
        let b = Matrix4::new([
            [-2.0, 1.0, 2.0, 3.0],
            [3.0, 2.0, 1.0, -1.0],
            [4.0, 3.0, 6.0, 5.0],
            [1.0, 2.0, 7.0, 8.0],
        ]);
        assert_eq!(
            (a * b).m,
            [
                [20.0, 22.0, 50.0, 48.0],
                [44.0, 54.0, 114.0, 108.0],
                [40.0, 58.0, 110.0, 102.0],
                [16.0, 26.0, 46.0, 42.0]
            ]
        );
    }

    #[test]
    fn matrix_tuple_multiplication() {
        let a = Matrix4::new([
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 4.0, 4.0, 2.0],
            [8.0, 6.0, 4.0, 1.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);
        let t = (1.0, 2.0, 3.0, 1.0);
        assert_eq!(a * t, (18.0, 24.0, 33.0, 1.0));
    }

    #[test]
    fn matrix_identity() {
        let a = Matrix4::new([
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 4.0, 4.0, 2.0],
            [8.0, 6.0, 4.0, 1.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);
        let t = (1.0, 2.0, 3.0, 4.0);

        assert_eq!(a * Matrix4::id(), a);
        assert_eq!(Matrix4::id() * t, t);
    }

    #[test]
    fn matrix_transposition() {
        let a = Matrix4::new([
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 4.0, 4.0, 2.0],
            [8.0, 6.0, 4.0, 1.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);

        assert_eq!(
            a.transpose(),
            Matrix4::new([
                [1.0, 2.0, 8.0, 0.0],
                [2.0, 4.0, 6.0, 0.0],
                [3.0, 4.0, 4.0, 0.0],
                [4.0, 2.0, 1.0, 1.0],
            ])
        );

        assert_eq!(Matrix4::id().transpose(), Matrix4::id());
    }

    #[test]
    fn matrix_determinant() {
        let a = Matrix2::new([[1.0, 5.0], [-3.0, 2.0]]);
        assert_eq!(a.determinant(), 17.0);

        let b = Matrix3::new([[1.0, 2.0, 6.0], [-5.0, 8.0, -4.0], [2.0, 6.0, 4.0]]);
        assert_eq!(b.cofactor(0, 0), 56.0);
        assert_eq!(b.cofactor(0, 1), 12.0);
        assert_eq!(b.cofactor(0, 2), -46.0);
        assert_eq!(b.determinant(), -196.0);

        let c = Matrix4::new([
            [-2.0, -8.0, 3.0, 5.0],
            [-3.0, 1.0, 7.0, 3.0],
            [1.0, 2.0, -9.0, 6.0],
            [-6.0, 7.0, 7.0, -9.0],
        ]);
        assert_eq!(c.cofactor(0, 0), 690.0);
        assert_eq!(c.cofactor(0, 1), 447.0);
        assert_eq!(c.cofactor(0, 2), 210.0);
        assert_eq!(c.cofactor(0, 3), 51.0);
        assert_eq!(c.determinant(), -4071.0);
    }

    #[test]
    fn matrix_submatrix() {
        let a = Matrix4::new([
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 4.0, 4.0, 2.0],
            [8.0, 6.0, 4.0, 1.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);

        assert_eq!(
            a.submatrix(1, 2),
            Matrix3::new([[1.0, 2.0, 4.0], [8.0, 6.0, 1.0], [0.0, 0.0, 1.0]])
        );

        let b = Matrix3::new([[1.0, 2.0, 4.0], [8.0, 6.0, 1.0], [0.0, 0.0, 1.0]]);

        assert_eq!(b.submatrix(1, 1), Matrix2::new([[1.0, 4.0], [0.0, 1.0]]));
    }

    #[test]
    fn matrix_minor() {
        let b = Matrix3::new([[3.0, 5.0, 0.0], [2.0, -1.0, -7.0], [6.0, -1.0, 5.0]]);
        assert_eq!(b.minor(1, 0), 25.0);
    }

    #[test]
    fn matrix_cofactor() {
        let b = Matrix3::new([[3.0, 5.0, 0.0], [2.0, -1.0, -7.0], [6.0, -1.0, 5.0]]);
        assert_eq!(b.cofactor(0, 0), -12.0);
        assert_eq!(b.cofactor(1, 0), -25.0);
    }

    #[test]
    fn matrix_inversion() {
        let a = Matrix4::new([
            [8.0, -5.0, 9.0, 2.0],
            [7.0, 5.0, 6.0, 1.0],
            [-6.0, 0.0, 9.0, 6.0],
            [-3.0, 0.0, -9.0, -4.0],
        ]);
        assert_eq!(
            a.inverse(),
            Matrix4::new([
                [-0.15384616, -0.15384616, -0.2820513, -0.53846157],
                [-0.07692308, 0.12307692, 0.025641026, 0.03076923],
                [0.35897437, 0.35897437, 0.43589744, 0.9230769],
                [-0.6923077, -0.6923077, -0.7692308, -1.9230769]
            ])
        );
    }

    #[test]
    fn matrix_inversion_property() {
        let a = Matrix4::new([
            [3.0, -9.0, 7.0, 3.0],
            [3.0, -8.0, 2.0, -9.0],
            [-4.0, 4.0, 4.0, 1.0],
            [-6.0, 5.0, -1.0, 1.0],
        ]);
        let b = Matrix4::new([
            [8.0, 2.0, 2.0, 2.0],
            [3.0, -1.0, 7.0, 0.0],
            [7.0, 0.0, 5.0, 4.0],
            [6.0, -2.0, 0.0, 5.0],
        ]);

        let c = a * b;
        assert_eq!(c * b.inverse(), a);
    }

    #[test]
    fn matrix_translation() {
        let transform = Matrix4::translation(5.0, -3.0, 2.0);
        let point = Point::new(-3.0, 4.0, 5.0);
        assert_eq!(transform * point, Point::new(2.0, 1.0, 7.0));

        let inv = transform.inverse();
        assert_eq!(inv * point, Point::new(-8.0, 7.0, 3.0));
    }

    #[test]
    fn matrix_translation_does_not_affect_vectors() {
        let transform = Matrix4::translation(5.0, -3.0, 2.0);
        let v = Vector::new(-3.0, 4.0, 5.0);
        assert_eq!(transform * v, v);
    }

    #[test]
    fn matrix_scaling() {
        let transform = Matrix4::scaling(2.0, 3.0, 4.0);
        let point = Point::new(-4.0, 6.0, 8.0);
        assert_eq!(transform * point, Point::new(-8.0, 18.0, 32.0));

        let vector = Vector::new(-4.0, 6.0, 8.0);
        assert_eq!(transform * vector, Vector::new(-8.0, 18.0, 32.0));

        let inv = transform.inverse();
        assert_eq!(inv * vector, Vector::new(-2.0, 2.0, 2.0));
    }

    #[test]
    fn matrix_reflection() {
        let transform = Matrix4::scaling(-1.0, 1.0, 1.0);
        let point = Point::new(2.0, 3.0, 4.0);
        assert_eq!(transform * point, Point::new(-2.0, 3.0, 4.0));
    }

    #[test]
    fn matrix_rotation_x() {
        let point = Point::new(0.0, 1.0, 0.0);
        let half_quarter = Matrix4::rotation_x(PI / 4.0);
        let full_quarter = Matrix4::rotation_x(PI / 2.0);
        assert_eq!(
            half_quarter * point,
            Point::new(0.0, 2_f32.sqrt() / 2.0, 2_f32.sqrt() / 2.0)
        );
        assert_eq!(full_quarter * point, Point::new(0.0, 0.0, 1.0));

        let inv = half_quarter.inverse();
        assert_eq!(
            inv * point,
            Point::new(0.0, 2_f32.sqrt() / 2.0, -2_f32.sqrt() / 2.0)
        );
    }

    #[test]
    fn matrix_rotation_y() {
        let point = Point::new(0.0, 0.0, 1.0);
        let half_quarter = Matrix4::rotation_y(PI / 4.0);
        let full_quarter = Matrix4::rotation_y(PI / 2.0);
        assert_eq!(
            half_quarter * point,
            Point::new(2_f32.sqrt() / 2.0, 0.0, 2_f32.sqrt() / 2.0)
        );
        assert_eq!(full_quarter * point, Point::new(1.0, 0.0, 0.0));
    }

    #[test]
    fn matrix_rotation_z() {
        let point = Point::new(0.0, 1.0, 0.0);
        let half_quarter = Matrix4::rotation_z(PI / 4.0);
        let full_quarter = Matrix4::rotation_z(PI / 2.0);
        assert_eq!(
            half_quarter * point,
            Point::new(-2_f32.sqrt() / 2.0, 2_f32.sqrt() / 2.0, 0.0)
        );
        assert_eq!(full_quarter * point, Point::new(-1.0, 0.0, 0.0));
    }

    #[test]
    fn matrix_shearing() {
        let point = Point::new(2.0, 3.0, 4.0);

        let mut transform = Matrix4::shearing(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        assert_eq!(transform * point, Point::new(5.0, 3.0, 4.0));

        transform = Matrix4::shearing(0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
        assert_eq!(transform * point, Point::new(6.0, 3.0, 4.0));

        transform = Matrix4::shearing(0.0, 0.0, 1.0, 0.0, 0.0, 0.0);
        assert_eq!(transform * point, Point::new(2.0, 5.0, 4.0));

        transform = Matrix4::shearing(0.0, 0.0, 0.0, 1.0, 0.0, 0.0);
        assert_eq!(transform * point, Point::new(2.0, 7.0, 4.0));

        transform = Matrix4::shearing(0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
        assert_eq!(transform * point, Point::new(2.0, 3.0, 6.0));

        transform = Matrix4::shearing(0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        assert_eq!(transform * point, Point::new(2.0, 3.0, 7.0));
    }

    #[test]
    fn vector_reflection_45_deg() {
        let v = Vector::new(1.0, -1.0, 0.0);
        let n = Vector::new(0.0, 1.0, 0.0);
        assert_eq!(v.reflect(&n), Vector::new(1.0, 1.0, 0.0));
    }

    #[test]
    fn vector_reflection_slanted_surface() {
        let v = Vector::new(0.0, -1.0, 0.0);
        let n = Vector::new(2_f32.sqrt() / 2.0, 2_f32.sqrt() / 2.0, 0.0);
        assert_eq!(v.reflect(&n), Vector::new(1.0, 0.0, 0.0));
    }
}
