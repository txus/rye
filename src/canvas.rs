use crate::primitives::Color;

pub trait Canvas {
    fn width(&self) -> usize;
    fn height(&self) -> usize;
    fn write(&mut self, y: usize, x: usize, color: Color);
    fn at(&self, y: usize, x: usize) -> Color;
}

pub struct DynamicCanvas {
    width: usize,
    height: usize,
    pixels: Vec<Vec<Color>>,
}

impl DynamicCanvas {
    pub fn new(width: usize, height: usize) -> Self {
        let pixels = vec![vec![Color::black(); width]; height];
        Self {
            width,
            height,
            pixels,
        }
    }
}

impl Canvas for DynamicCanvas {
    fn width(&self) -> usize {
        self.width
    }
    fn height(&self) -> usize {
        self.height
    }
    fn write(&mut self, y: usize, x: usize, color: Color) {
        assert!(y < self.height && x < self.width);
        self.pixels[y][x] = color;
    }

    fn at(&self, y: usize, x: usize) -> Color {
        assert!(
            y < self.height && x < self.width,
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
        let c = DynamicCanvas::new(5, 5);
        for x in c.pixels.iter() {
            for &y in x.iter() {
                assert_eq!(y, Color::new(0.0, 0.0, 0.0));
            }
        }
    }

    #[test]
    fn write_pixel() {
        let mut c = DynamicCanvas::new(5, 5);
        c.write(2, 4, Color::new(1.0, 0.0, 0.0));
        assert_eq!(c.at(2, 4), Color::new(1.0, 0.0, 0.0));
    }
}
