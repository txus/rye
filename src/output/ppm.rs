use super::{Canvas, Output};
use crate::primitives::Color;
use std::cmp;

pub struct PPM;

fn scale(color: Color) -> (u8, u8, u8) {
    let r = color.red * 255.0;
    let g = color.green * 255.0;
    let b = color.blue * 255.0;
    (
        if r > 255.0 {
            255
        } else if r < 0.0 {
            0
        } else {
            r as u8
        },
        if g > 255.0 {
            255
        } else if g < 0.0 {
            0
        } else {
            g as u8
        },
        if b > 255.0 {
            255
        } else if b < 0.0 {
            0
        } else {
            b as u8
        },
    )
}

impl Output for PPM {
    fn render<C: Canvas>(&self, c: &C) -> String {
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
    use crate::canvas::DynamicCanvas;
    use crate::primitives::Color;

    #[test]
    fn scale_test() {
        assert_eq!(scale(Color::new(2.0, 0.5, 1.0)), (255, 127, 255));
    }

    #[test]
    fn render() {
        let f = PPM {};
        let mut c = DynamicCanvas::new(5, 5);
        c.write(2, 4, Color::new(1.0, 0.5, 0.0));
        let output = f.render(&c);
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
