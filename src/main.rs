mod canvas;
mod light;
mod materials;
mod output;
mod primitives;
mod rays;
mod shapes;
mod world;
use light::PointLight;
use materials::Material;
use primitives::{Color, Matrix4, Point, Vector};
use shapes::{Shape, Plane, Sphere};
use world::{view_transform, Camera, World};

use std::f32::consts::PI;

fn main() {
    let mut floor: Box<Shape> = Box::from(Plane::new());
    floor.set_material(Material {
        color: Color::new(1.0, 0.9, 0.9),
        specular: 0.0,
        ..Material::default()
    });

    let mut wall: Box<Shape> = Box::from(Plane::new());
    wall.set_transform(Matrix4::translation(0.0, 0.0, 10.0) * Matrix4::rotation_x(PI/2.0) * Matrix4::rotation_y(PI/2.0));
    wall.set_material(Material {
        color: Color::new(0.5, 0.9, 0.9),
        specular: 0.5,
        ..Material::default()
    });

    let mut middle: Box<Shape> = Box::from(Sphere::unit());
    middle.set_transform(Matrix4::translation(-0.5, 0.0, 0.5));
    middle.set_material(Material {
        color: Color::new(0.1, 1.0, 0.5),
        diffuse: 0.7,
        specular: 0.3,
        ..Material::default()
    });

    let mut right: Box<Shape> = Box::from(Sphere::unit());
    right.set_transform(Matrix4::translation(1.5, 0.5, -0.5) * Matrix4::scaling(0.5, 0.5, 0.5));
    right.set_material(Material {
        color: Color::new(0.5, 1.0, 0.1),
        diffuse: 0.7,
        specular: 0.3,
        ..Material::default()
    });

    let mut left: Box<Shape> = Box::from(Sphere::unit());
    left.set_transform(
        Matrix4::translation(-1.5, 0.33, -0.75) * Matrix4::scaling(0.33, 0.33, 0.33),
    );
    left.set_material(Material {
        color: Color::new(1.0, 0.8, 0.1),
        diffuse: 0.7,
        specular: 0.3,
        ..Material::default()
    });

    let mut world = World::default();
    world.light_source = PointLight::new(Point::new(-10.0, 10.0, -10.0), Color::new(1.0, 1.0, 1.0));
    world.objects = vec![floor, wall, middle, right, left];

    let mut camera = Camera::new(100, 50, PI / 3.0);
    camera.transform = view_transform(
        Point::new(0.0, 1.5, -5.0),
        Point::new(0.0, 1.0, 0.0),
        Vector::new(0.0, 1.0, 0.0),
    );

    let canvas = camera.render(&world);

    output::render(&canvas, output::ppm::PPM {}, "/tmp/scene.ppm");
}
