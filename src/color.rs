use crate::linear::EPSILON;

#[derive(Debug, Clone, Copy)]
pub struct Color {
    pub red: f32,
    pub green: f32,
    pub blue: f32,
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

fn close(a: f32, b: f32) -> bool {
    (a - b).abs() < (EPSILON * 2.0)
}

impl std::cmp::PartialEq for Color {
    fn eq(&self, other: &Self) -> bool {
        close(self.red, other.red) && close(self.green, other.green) && close(self.blue, other.blue)
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

impl std::ops::Div<f32> for Color {
    type Output = Self;
    fn div(self, scalar: f32) -> Self {
        Self::new(self.red / scalar, self.green / scalar, self.blue / scalar)
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
    fn initialize() {
        let c = Color::new(-0.5, 0.4, 1.7);
        assert_eq!(c.red, -0.5);
        assert_eq!(c.green, 0.4);
        assert_eq!(c.blue, 1.7);
    }

    #[test]
    fn addition() {
        let v1 = Color::new(0.9, 0.6, 0.75);
        let v2 = Color::new(0.7, 0.1, 0.25);

        assert_eq!(v1 + v2, Color::new(1.6, 0.7, 1.0));
    }

    #[test]
    fn subtraction() {
        let v1 = Color::new(0.9, 0.6, 0.75);
        let v2 = Color::new(0.7, 0.1, 0.25);

        assert_eq!(v1 - v2, Color::new(0.2, 0.5, 0.5));
    }

    #[test]
    fn scaling() {
        let v = Color::new(0.2, 0.3, 0.4);

        assert_eq!(v * 2.0, Color::new(0.4, 0.6, 0.8));
    }

    #[test]
    fn hadamard_product() {
        let v1 = Color::new(1.0, 0.2, 0.4);
        let v2 = Color::new(0.9, 1.0, 0.1);

        assert_eq!(v1 * v2, Color::new(0.9, 0.2, 0.04));
    }
}
