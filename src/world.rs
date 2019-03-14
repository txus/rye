use crate::canvas::{Canvas, DynamicCanvas};
use crate::color::Color;
use crate::light::PointLight;
use crate::linear::{Matrix, Matrix4, Point, Vector};
use crate::materials::Material;
use crate::rays::{Intersection, Precomputation, Ray};
use crate::shapes::{Shape, Sphere};

pub struct World {
    pub objects: Vec<Box<Shape>>,
    pub light_source: PointLight,
}

impl World {
    pub fn default() -> Self {
        let light = PointLight::new(Point::new(-10.0, 10.0, -10.0), Color::white());
        let mut s1: Box<Shape> = Box::from(Sphere::new());
        s1.set_material(Material {
            color: Color::new(0.8, 1.0, 0.6),
            diffuse: 0.7,
            specular: 0.2,
            ..Material::default()
        });
        let mut s2: Box<Shape> = Box::from(Sphere::new());
        s2.set_transform(Matrix4::scaling(0.5, 0.5, 0.5));

        World {
            objects: vec![s1, s2],
            light_source: light,
        }
    }

    pub fn intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let mut out: Vec<Intersection> = vec![];
        for object in &self.objects {
            out.append(&mut object.intersect(&ray))
        }
        out.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
        out
    }

    pub fn shade(&self, c: &Precomputation) -> Color {
        c.object.material().lighting(
            c.object,
            &self.light_source,
            &c.over_point,
            &c.eye,
            &c.normal,
            self.is_shadowed(c.over_point),
        )
    }

    pub fn color_at(&self, r: &Ray) -> Color {
        if let Some(i) = self.intersect(&r).first() {
            self.shade(&i.precompute(&r))
        } else {
            Color::black()
        }
    }

    pub fn is_shadowed(&self, p: Point) -> bool {
        let v = self.light_source.position - p;
        let distance = v.magnitude();
        let direction = v.normalize();

        let r = Ray::new(p, direction);
        let mut intersections = self.intersect(&r);

        match Intersection::hit(&mut intersections) {
            Some(h) if h.t < distance => true,
            _ => false,
        }
    }
}

pub struct Camera {
    hsize: u32,
    vsize: u32,
    fov: f32,
    pub transform: Matrix4,
    half_width: f32,
    half_height: f32,
    pixel_size: f32,
}

impl Camera {
    pub fn new(hsize: u32, vsize: u32, fov: f32) -> Camera {
        let half_view = (fov / 2.0).tan();
        let aspect_ratio = hsize as f32 / vsize as f32;
        let half_width = if aspect_ratio >= 1.0 {
            half_view
        } else {
            half_view * aspect_ratio
        };
        let half_height = if aspect_ratio >= 1.0 {
            half_view / aspect_ratio
        } else {
            half_view
        };
        let pixel_size = (half_width * 2.0) / hsize as f32;

        Camera {
            hsize,
            vsize,
            fov,
            pixel_size,
            half_width,
            half_height,
            transform: Matrix4::id(),
        }
    }

    pub fn ray_for_pixel(&self, px: u32, py: u32) -> Ray {
        let pixel_size = self.pixel_size;
        let xoffset = (px as f32 + 0.5) * pixel_size;
        let yoffset = (py as f32 + 0.5) * pixel_size;

        let world_x = self.half_width - xoffset;
        let world_y = self.half_height - yoffset;

        let inv = self.transform.inverse();

        let pixel = inv * Point::new(world_x, world_y, -1.0);
        let origin = inv * Point::origin();
        let direction = (pixel - origin).normalize();

        Ray::new(origin, direction)
    }

    pub fn render(&self, world: &World) -> DynamicCanvas {
        let mut canvas = DynamicCanvas::new(self.hsize as usize, self.vsize as usize);
        for y in 0..self.vsize {
            for x in 0..self.hsize {
                let ray = self.ray_for_pixel(x, y);
                let color = world.color_at(&ray);
                canvas.write(y as usize, x as usize, color);
            }
        }
        canvas
    }
}

pub fn view_transform(from: Point, to: Point, up: Vector) -> Matrix4 {
    let forward = (to - from).normalize();
    let left = forward.cross(&up.normalize());
    let true_up = left.cross(&forward);
    let orientation = Matrix4::new([
        [left.x, left.y, left.z, 0.0],
        [true_up.x, true_up.y, true_up.z, 0.0],
        [-forward.x, -forward.y, -forward.z, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]);
    orientation * Matrix4::translation(-from.x, -from.y, -from.z)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::color::Color;
    use crate::light::PointLight;
    use crate::linear::{Matrix4, Point, EPSILON};
    use crate::materials::Material;

    use std::f32::consts::PI;

    #[test]
    fn default_world() {
        let light = PointLight::new(Point::new(-10.0, 10.0, -10.0), Color::white());
        let mut s1 = Sphere::new();
        s1.material = Material {
            color: Color::new(0.8, 1.0, 0.6),
            diffuse: 0.7,
            specular: 0.2,
            ..s1.material
        };
        let mut s2 = Sphere::new();
        s2.set_transform(Matrix4::scaling(0.5, 0.5, 0.5));

        let w = World::default();
        assert_eq!(w.light_source, light);
    }

    #[test]
    fn intersect_world() {
        let w = World::default();
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let hits = w.intersect(&r).iter().map(|x| x.t).collect::<Vec<f32>>();
        assert_eq!(hits, vec!(4.0, 4.5, 5.5, 6.0));
    }

    #[test]
    fn shading_intersection() {
        let w = World::default();
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let shape: &Box<Shape> = w.objects.first().unwrap();
        let i = Intersection {
            t: 4.0,
            object: &**shape,
        };
        let comps = i.precompute(&r);
        assert_eq!(w.shade(&comps), Color::new(0.38066, 0.47583, 0.2855));
    }

    #[test]
    fn shading_intersection_from_inside() {
        let mut w = World::default();
        w.light_source = PointLight::new(Point::new(0.0, 0.25, 0.0), Color::new(1.0, 1.0, 1.0));
        let r = Ray::new(Point::origin(), Vector::new(0.0, 0.0, 1.0));
        let shape: &Box<Shape> = w.objects.last().unwrap();
        let i = Intersection {
            t: 0.5,
            object: &**shape,
        };
        let comps = i.precompute(&r);
        assert_eq!(w.shade(&comps), Color::new(0.90498, 0.90498, 0.90498));
    }

    #[test]
    fn color_at_miss() {
        let w = World::default();
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 1.0, 0.0));
        let c = w.color_at(&r);
        assert_eq!(c, Color::new(0.0, 0.0, 0.0));
    }

    #[test]
    fn color_at_hit() {
        let w = World::default();
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let c = w.color_at(&r);
        assert_eq!(c, Color::new(0.38066, 0.47582, 0.2855));
    }

    #[test]
    fn view_transform_default() {
        let from = Point::origin();
        let to = Point::new(0.0, 0.0, -10.0);
        let up = Vector::new(0.0, 1.0, 0.0);
        assert_eq!(view_transform(from, to, up), Matrix4::id());
    }

    #[test]
    fn view_transform_looking_pos_z_direction() {
        let from = Point::origin();
        let to = Point::new(0.0, 0.0, 1.0);
        let up = Vector::new(0.0, 1.0, 0.0);
        assert_eq!(
            view_transform(from, to, up),
            Matrix4::scaling(-1.0, 1.0, -1.0)
        );
    }

    #[test]
    fn view_transform_moves_world() {
        let from = Point::new(0.0, 0.0, 8.0);
        let to = Point::origin();
        let up = Vector::new(0.0, 1.0, 0.0);
        assert_eq!(
            view_transform(from, to, up),
            Matrix4::translation(0.0, 0.0, -8.0)
        );
    }

    #[test]
    fn view_transform_arbitrary() {
        let from = Point::new(1.0, 3.0, 2.0);
        let to = Point::new(4.0, -2.0, 8.0);
        let up = Vector::new(1.0, 1.0, 0.0);
        assert_eq!(
            view_transform(from, to, up),
            Matrix4::new([
                [-0.50709, 0.50709, 0.67612, -2.36643],
                [0.76772, 0.60609, 0.12122, -2.82843],
                [-0.35857, 0.59761, -0.71714, 0.00000],
                [0.00000, 0.00000, 0.00000, 1.00000]
            ])
        );
    }

    #[test]
    fn camera_pixel_size() {
        let c = Camera::new(200, 125, PI / 2.0);
        assert_eq!(c.pixel_size, 0.01);

        let c2 = Camera::new(125, 200, PI / 2.0);
        assert_eq!(c2.pixel_size, 0.01);
    }

    #[test]
    fn ray_through_center() {
        let c = Camera::new(201, 101, PI / 2.0);
        let r = c.ray_for_pixel(100, 50);
        assert_eq!(r.origin, Point::origin());
        assert_eq!(r.direction, Vector::new(0.0, 0.0, -1.0));
    }

    #[test]
    fn ray_through_corner() {
        let c = Camera::new(201, 101, PI / 2.0);
        let r = c.ray_for_pixel(0, 0);
        assert_eq!(r.origin, Point::origin());
        assert_eq!(r.direction, Vector::new(0.66519, 0.33259, -0.66851));
    }

    #[test]
    fn ray_with_transformed_camera() {
        let mut c = Camera::new(201, 101, PI / 2.0);
        c.transform = Matrix4::rotation_y(PI / 4.0) * Matrix4::translation(0.0, -2.0, 5.0);
        let r = c.ray_for_pixel(100, 50);
        assert_eq!(r.origin, Point::new(0.0, 2.0, -5.0));
        assert_eq!(
            r.direction,
            Vector::new(2_f32.sqrt() / 2.0, 0.0, -2_f32.sqrt() / 2.0)
        );
    }

    #[test]
    fn camera_render_world() {
        let w = World::default();
        let from = Point::new(0.0, 0.0, -5.0);
        let to = Point::origin();
        let up = Vector::new(0.0, 1.0, 0.0);
        let mut c = Camera::new(11, 11, PI / 2.0);
        c.transform = view_transform(from, to, up);
        let image = c.render(&w);
        assert_eq!(image.at(5, 5), Color::new(0.38066, 0.47583, 0.2855));
    }

    #[test]
    fn no_shadow_when_nothing_is_collinear_with_point_and_light() {
        let w = World::default();
        let p = Point::new(0.0, 10.0, 0.0);
        assert_eq!(w.is_shadowed(p), false);
    }

    #[test]
    fn shadow_when_object_is_between_light_and_point() {
        let w = World::default();
        let p = Point::new(10.0, -10.0, 10.0);
        assert_eq!(w.is_shadowed(p), true);
    }

    #[test]
    fn no_shadow_when_object_is_behind_the_light() {
        let w = World::default();
        let p = Point::new(-20.0, 20.0, -20.0);
        assert_eq!(w.is_shadowed(p), false);
    }

    #[test]
    fn no_shadow_when_object_is_behind_the_point() {
        let w = World::default();
        let p = Point::new(-2.0, 2.0, -2.0);
        assert_eq!(w.is_shadowed(p), false);
    }

    #[test]
    fn hit_should_offset_the_point() {
        let r = Ray::new(Point::new(0.0, 0.0, -5.0), Vector::new(0.0, 0.0, 1.0));
        let mut s = Sphere::new();
        s.set_transform(Matrix4::translation(0.0, 0.0, 1.0));
        let i = Intersection { t: 5.0, object: &s };
        let comps = i.precompute(&r);
        assert!(comps.over_point.z < -EPSILON / 2.0);
        assert!(comps.point.z > comps.over_point.z);
    }
}
