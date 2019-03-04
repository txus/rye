mod canvas;
mod output;
mod primitives;
use canvas::{Canvas, LargeCanvas};
use primitives::{Color, Point, Vector};

// Test

#[derive(Debug)]
struct Projectile {
    pos: Point,
    velocity: Vector,
}

struct Environment {
    gravity: Vector,
    wind: Vector,
}

fn tick(env: &Environment, proj: &Projectile) -> Projectile {
    let position = proj.pos + proj.velocity;
    let velocity = proj.velocity + env.gravity + env.wind;
    Projectile {
        pos: position,
        velocity: velocity,
    }
}

fn main() {
    let start = Point::new(0.0, 1.0, 0.0);
    let velocity = Vector::new(1.0, 1.8, 0.0).normalize() * 11.25;
    let mut p = Projectile {
        pos: start,
        velocity: velocity,
    };
    let e = Environment {
        gravity: Vector::new(0.0, -0.1, 0.0),
        wind: Vector::new(-0.01, 0.0, 0.0),
    };

    let mut canvas = LargeCanvas::new();

    while p.pos.y > 0.0 {
        p = tick(&e, &p);
        let y = canvas.height() as f32 - p.pos.y;
        let x = p.pos.x;
        if y > 0 as f32 && x < canvas.width() as f32 {
            canvas.write(y as usize, x as usize, Color::new(1.0, 0.0, 0.0));
        }
    }
    let ppm = output::ppm::PPM {};
    output::render(canvas, ppm, "/tmp/output.ppm");
}

#[cfg(test)]
mod tests {}
