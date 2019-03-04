use super::{Canvas, Output};
use crate::primitives::Color;
use std::cmp;

pub struct PPM;

fn clamp(min: u8, max: u8, number: u8) -> u8 {
    cmp::max(min, cmp::min(max, number))
}

fn scale(color: Color) -> (u8, u8, u8) {
    return (
        clamp(0, 255, (color.red * 255.0) as u8),
        clamp(0, 255, (color.green * 255.0) as u8),
        clamp(0, 255, (color.blue * 255.0) as u8),
    );
}

impl Output for PPM {
    fn render<C: Canvas>(&self, c: C) -> String {
        let w = c.width();
        let h = c.height();
        let mut s = String::new();
        s.push_str("P3\n");
        s.push_str(&format!("{} {}\n", w, h));
        s.push_str("255\n");
        for y in 0..h {
            for x in 0..w {
                let (red, green, blue) = scale(c.at(y, x));
                s.push_str(&format!("{} {} {}\n", red, green, blue));
            }
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::canvas::SmallCanvas;
    use crate::primitives::Color;

    #[test]
    fn render() {
        let f = PPM {};
        let mut c = SmallCanvas::new();
        c.write(2, 4, Color::new(1.0, 0.5, 0.0));
        let output = f.render(c);
        let mut lines = output.lines();
        assert_eq!(lines.next(), Some("P3"));
        assert_eq!(lines.next(), Some("5 5"));
        assert_eq!(lines.next(), Some("255"));
        for y in 0..5 {
            for x in 0..5 {
                if y == 2 && x == 4 {
                    assert_eq!(
                        lines.next(),
                        Some("255 127 0"),
                        "Pixel at {} {} should be colored",
                        y,
                        x
                    );
                } else {
                    assert_eq!(
                        lines.next(),
                        Some("0 0 0"),
                        "Pixel at {} {} should be black",
                        y,
                        x
                    );
                }
            }
        }
    }
}
