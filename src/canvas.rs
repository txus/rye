use crate::primitives::Color;

pub trait Canvas {
    fn width(&self) -> usize;
    fn height(&self) -> usize;
    fn write(&mut self, y: usize, x: usize, color: Color);
    fn at(&self, y: usize, x: usize) -> Color;
}

pub struct SmallCanvas {
    pixels: Box<[[Color; 5]; 5]>,
}

pub struct LargeCanvas {
    pixels: Box<[[Color; 800]; 600]>,
}

impl SmallCanvas {
    pub fn new() -> Self {
        let pixels = Box::new([[Color::new(0.0, 0.0, 0.0); 5]; 5]);
        Self { pixels: pixels }
    }
}

impl LargeCanvas {
    pub fn new() -> Self {
        let pixels = Box::new([[Color::new(0.0, 0.0, 0.0); 800]; 600]);
        Self { pixels: pixels }
    }
}

impl Canvas for SmallCanvas {
    fn width(&self) -> usize {
        5
    }
    fn height(&self) -> usize {
        5
    }

    fn write(&mut self, y: usize, x: usize, color: Color) {
        assert!(y < 5 && x < 5);
        self.pixels[y][x] = color;
    }

    fn at(&self, y: usize, x: usize) -> Color {
        assert!(
            y < 5 && x < 5,
            "Can't read X {} Y {} -- out of bounds",
            x,
            y
        );
        self.pixels[y][x]
    }
}

impl Canvas for LargeCanvas {
    fn width(&self) -> usize {
        800
    }
    fn height(&self) -> usize {
        600
    }

    fn write(&mut self, y: usize, x: usize, color: Color) {
        assert!(y < 600 && x < 800);
        self.pixels[y][x] = color;
    }

    fn at(&self, y: usize, x: usize) -> Color {
        assert!(
            y < 600 && x < 800,
            "Can't read X {} Y {} -- out of bounds",
            x,
            y
        );
        self.pixels[y][x]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn initialize_to_black() {
        let c = SmallCanvas::new();
        for x in c.pixels.iter() {
            for &y in x.iter() {
                assert_eq!(y, Color::new(0.0, 0.0, 0.0));
            }
        }
    }

    #[test]
    fn write_pixel() {
        let mut c = SmallCanvas::new();
        c.write(2, 4, Color::new(1.0, 0.0, 0.0));
        assert_eq!(c.at(2, 4), Color::new(1.0, 0.0, 0.0));
    }
}
