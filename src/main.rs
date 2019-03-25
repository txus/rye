mod camera;
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
use camera::{view_transform, Camera};
use color::Color;
use light::PointLight;
use linear::{Matrix4, Point, Vector};
use materials::Material;
use patterns::{CheckerPattern, GradientPattern, RingPattern, StripePattern};
use shapes::{Cube, Cylinder, Plane, Shape, Sphere};
use world::World;

use std::f32::consts::PI;

fn scene1() {
    let mut floor: Box<Shape> = Box::from(Plane::new());
    floor.set_material(Material {
        color: Color::new(1.0, 0.9, 0.9),
        specular: 0.0,
        reflective: 0.5,
        pattern: Some(Box::from(CheckerPattern::new(
            Color::white(),
            Color::black(),
        ))),
        ..Material::default()
    });

    let mut wall: Box<Shape> = Box::from(Plane::new());
    wall.set_transform(Matrix4::rotation_x(PI / 2.0) * Matrix4::translation(0.0, 0.0, 8.0));
    wall.set_material(Material {
        color: Color::new(1.0, 0.9, 0.9),
        specular: 0.0,
        pattern: Some(Box::from(CheckerPattern::new(
            Color::white(),
            Color::black(),
        ))),
        ..Material::default()
    });

    let mut right_wall: Box<Shape> = Box::from(Plane::new());
    right_wall.set_transform(
        Matrix4::rotation_x(PI / 2.0)
            * Matrix4::rotation_z(PI / 2.0)
            * Matrix4::translation(-15.0, 0.0, 0.0),
    );
    right_wall.set_material(Material {
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
        reflective: 1.0,
        pattern: Some(Box::from(GradientPattern::new(
            Color::green(),
            Color::blue(),
        ))),
        ..Material::default()
    });

    let mut right: Box<Shape> = Box::from(Sphere::new());
    right.set_transform(Matrix4::translation(1.5, 0.5, -0.5) * Matrix4::scaling(0.5, 0.5, 0.5));
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
    world.light_source = PointLight::new(Point::new(-10.0, 10.0, -10.0), Color::white());
    world.objects = vec![middle, right, left, floor, wall];

    let mut camera = Camera::new(800, 600, PI / 1.8);
    camera.transform = view_transform(
        Point::new(-2.0, 1.5, -3.5),
        Point::new(1.0, 1.0, 0.0),
        Vector::new(0.0, 0.8, 0.3),
    );

    let canvas = camera.render(&world);

    output::render(&canvas, output::ppm::PPM {}, "/tmp/scene.ppm");
}

fn main() {
    let mut cylinder = Cylinder::closed(-1.5, 0.5);
    cylinder.set_transform(Matrix4::translation(0.0, 1.5, 1.5));
    cylinder.casts_shadows = false;
    cylinder.material.shininess = 0.1;
    cylinder.material.reflective = 1.0;
    cylinder.material.color = Color::blue();

    let mut floor = Plane::new();
    floor.material.pattern = Some(Box::from(CheckerPattern::new(
        Color::white(),
        Color::black(),
    )));

    let mut ball = Sphere::new();
    ball.set_material(Material {
        color: Color::new(1.0, 0.0, 0.0),
        ambient: 0.5,
        ..Material::default()
    });
    ball.set_transform(Matrix4::translation(-1.5, 0.5, -0.5));

    let mut world = World::default();
    world.objects = vec![Box::from(floor), Box::from(cylinder)];

    let mut camera = Camera::new(200, 150, PI / 2.2);
    camera.transform = view_transform(
        Point::new(0.0, 2.5, -3.5),
        Point::new(0.0, 0.0, 3.0),
        Vector::new(0.0, 1.0, 0.0),
    );

    let canvas = camera.render(&world);

    output::render(&canvas, output::ppm::PPM {}, "/tmp/refraction.ppm");
}
