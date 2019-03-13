pub mod ppm;

use std::fs;

use crate::canvas::Canvas;

pub trait Output {
    fn render<C: Canvas>(&self, c: &C) -> String;
}

pub fn render<C: Canvas, O: Output>(c: &C, o: O, filename: &str) {
    fs::write(filename, o.render(c)).expect("Unable to write file")
}
