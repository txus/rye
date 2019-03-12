mod canvas;
mod output;
mod primitives;
mod rays;
mod light;
mod materials;
mod world;
mod objects;
use canvas::{Canvas, LargeCanvas, MediumCanvas};
use primitives::{Color, Point, Vector, Matrix4};
use materials::Material;
use light::PointLight;
use rays::{Ray, Sphere, Intersection};

// Projectile test

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
    let pos = proj.pos + proj.velocity;
    let velocity = proj.velocity + env.gravity + env.wind;
    Projectile { pos, velocity }
}

fn projectile_test() {
    let start = Point::new(0.0, 1.0, 0.0);
    let velocity = Vector::new(1.0, 1.8, 0.0).normalize() * 11.25;
    let mut p = Projectile {
        pos: start,
        velocity,
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
    output::render(canvas, ppm, "/tmp/projectile.ppm");
}

// Rendering a sphere test

fn sphere_test() {
    let ray_origin = Point::new(0.0, 0.0, -5.0);
    let wall_z = 10.0;
    let wall_size = 7.0;
    let mut canvas = MediumCanvas::new();
 
    let pixel_size = wall_size / canvas.width() as f32;

    let half = wall_size / 2.0;

    let mut shape = Sphere::unit();

    let mut material = Material::default();
    material.color = Color::new(1.0, 0.2, 1.0);
    shape.material = material;

    let light_position = Point::new(-10.0, 10.0, -10.0);
    let light_color = Color::white();
    let light = PointLight::new(light_position, light_color);

    shape.transform = Matrix4::scaling(0.5, 1.0, 1.0) * Matrix4::shearing(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    for y in 0..100 {
        let world_y = half - pixel_size * y as f32;

        for x in 0..100 {
          let world_x = -half + pixel_size * x as f32;

          let position = Point::new(world_x, world_y, wall_z);
          let r = Ray::new(ray_origin, (position - ray_origin).normalize());

          if let Some(hit) = Intersection::hit(&mut shape.intersect(&r)) {
              let point = r.position(hit.t);
              let normal = hit.object.normal(point);
              let eye = -r.direction;
              let color = hit.object.material.lighting(&light, &point, &eye, &normal);
              canvas.write(y, x, color);
          }
        }
    }

    let ppm = output::ppm::PPM {};
    output::render(canvas, ppm, "/tmp/sphere.ppm");
}

fn main() {
    sphere_test();
}
