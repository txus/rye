mod canvas;
mod color;
mod light;
mod linear;
mod materials;
mod output;
mod patterns;
mod rays;
mod shapes;
mod world;
use color::Color;
use light::PointLight;
use linear::{Matrix4, Point, Vector};
use materials::Material;
use patterns::{CheckerPattern, GradientPattern, RingPattern, StripePattern};
use shapes::{Plane, Shape, Sphere};
use world::{view_transform, Camera, World};

use std::f32::consts::PI;

fn main() {
    let mut floor: Box<Shape> = Box::from(Plane::new());
    floor.set_material(Material {
        color: Color::new(1.0, 0.9, 0.9),
        specular: 0.0,
        pattern: Some(Box::from(CheckerPattern::new(
            Color::white(),
            Color::black(),
        ))),
        ..Material::default()
    });

    let mut middle: Box<Shape> = Box::from(Sphere::new());
    middle.set_transform(Matrix4::translation(-0.5, 1.0, 0.5));
    middle.set_material(Material {
        color: Color::new(0.1, 1.0, 0.5),
        diffuse: 0.7,
        specular: 0.3,
        pattern: Some(Box::from(GradientPattern::new(
            Color::green(),
            Color::blue(),
        ))),
        ..Material::default()
    });

    let mut right: Box<Shape> = Box::from(Sphere::new());
    right.set_transform(Matrix4::translation(1.5, 0.5, -0.5) * Matrix4::scaling(2.5, 2.5, 2.5));
    right.set_material(Material {
        color: Color::new(0.5, 1.0, 0.1),
        diffuse: 0.7,
        specular: 0.3,
        pattern: Some(Box::from(RingPattern::new(Color::red(), Color::white()))),
        ..Material::default()
    });

    let mut left: Box<Shape> = Box::from(Sphere::new());
    left.set_transform(
        Matrix4::translation(-1.5, 0.33, -0.75) * Matrix4::scaling(0.33, 0.33, 0.33),
    );
    left.set_material(Material {
        color: Color::new(1.0, 0.8, 0.1),
        diffuse: 0.7,
        specular: 0.3,
        pattern: Some(Box::from(StripePattern::new(Color::green(), Color::blue()))),
        ..Material::default()
    });

    let mut world = World::default();
    world.light_source = PointLight::new(Point::new(-10.0, 10.0, -10.0), Color::new(1.0, 1.0, 1.0));
    world.objects = vec![middle, right, left];

    let mut camera = Camera::new(200, 100, PI / 3.0);
    camera.transform = view_transform(
        Point::new(0.0, 1.5, -5.0),
        Point::new(0.0, 1.0, 0.0),
        Vector::new(0.0, 1.0, 0.0),
    );

    let canvas = camera.render(&world);

    output::render(&canvas, output::ppm::PPM {}, "/tmp/scene.ppm");
}
