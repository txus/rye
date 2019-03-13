mod canvas;
mod light;
mod materials;
mod objects;
mod output;
mod primitives;
mod rays;
mod world;
use canvas::Canvas;
use light::PointLight;
use materials::Material;
use primitives::{Color, Matrix4, Point, Vector};
use rays::{Intersection, Ray, Sphere};
use world::{view_transform, Camera, World};

use std::f32::consts::PI;

fn main() {
    let mut floor = Sphere::unit();
    floor.transform = Matrix4::scaling(10.0, 0.01, 10.0);
    floor.material = Material::default();
    floor.material.color = Color::new(1.0, 0.9, 0.9);
    floor.material.specular = 0.0;

    let mut left_wall = Sphere::unit();
    left_wall.transform = Matrix4::translation(0.0, 0.0, 5.0)
        * Matrix4::rotation_y(-PI / 4.0)
        * Matrix4::rotation_x(PI / 2.0)
        * Matrix4::scaling(10.0, 0.01, 10.0);
    left_wall.material = floor.material;

    let mut right_wall = Sphere::unit();
    right_wall.transform = Matrix4::translation(0.0, 0.0, 5.0)
        * Matrix4::rotation_y(PI / 4.0)
        * Matrix4::rotation_x(PI / 2.0)
        * Matrix4::scaling(10.0, 0.01, 10.0);
    right_wall.material = floor.material;

    let mut middle = Sphere::unit();
    middle.transform = Matrix4::translation(-0.5, 1.0, 0.5);
    middle.material = Material::default();
    middle.material.color = Color::new(0.1, 1.0, 0.5);
    middle.material.diffuse = 0.7;
    middle.material.specular = 0.3;

    let mut right = Sphere::unit();
    right.transform = Matrix4::translation(1.5, 0.5, -0.5) * Matrix4::scaling(0.5, 0.5, 0.5);
    right.material = Material::default();
    right.material.color = Color::new(0.5, 1.0, 0.1);
    right.material.diffuse = 0.7;
    right.material.specular = 0.3;

    let mut left = Sphere::unit();
    left.transform = Matrix4::translation(-1.5, 0.33, -0.75) * Matrix4::scaling(0.33, 0.33, 0.33);
    left.material = Material::default();
    left.material.color = Color::new(1.0, 0.8, 0.1);
    left.material.diffuse = 0.7;
    left.material.specular = 0.3;

    let mut world = World::default();
    world.light_source = PointLight::new(Point::new(-10.0, 10.0, -10.0), Color::new(1.0, 1.0, 1.0));
    world.objects = vec![floor, left_wall, right_wall, middle, right, left];

    let mut camera = Camera::new(200, 100, PI / 3.0);
    camera.transform = view_transform(
        Point::new(0.0, 1.5, -5.0),
        Point::new(0.0, 1.0, 0.0),
        Vector::new(0.0, 1.0, 0.0),
    );

    let canvas = camera.render(&world);

    output::render(&canvas, output::ppm::PPM {}, "/tmp/scene.ppm");
}
