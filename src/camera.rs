use crate::canvas::{Canvas, DynamicCanvas};
use crate::linear::{Matrix, Matrix4, Point, Vector};
use crate::rays::Ray;
use crate::world::{World, MAX_REFLECTIONS};
use rand::Rng;

pub struct Camera {
    hsize: u32,
    vsize: u32,
    pub transform: Matrix4,
    supersampling: usize,
    half_width: f32,
    half_height: f32,
    pixel_size: f32,
}

impl Camera {
    pub fn new(hsize: u32, vsize: u32, fov: f32, supersampling: usize) -> Camera {
        let half_view = (fov / 2.0).tan();
        let aspect_ratio = hsize as f32 / vsize as f32;
        let landscape = aspect_ratio >= 1.0;
        let half_width = if landscape {
            half_view
        } else {
            half_view * aspect_ratio
        };
        let half_height = if landscape {
            half_view / aspect_ratio
        } else {
            half_view
        };
        let pixel_size = (half_width * 2.0) / hsize as f32;

        Camera {
            hsize,
            vsize,
            supersampling,
            pixel_size,
            half_width,
            half_height,
            transform: Matrix4::id(),
        }
    }

    pub fn ray_for_pixel(&self, px: f32, py: f32) -> Ray {
        let xoffset = (px + 0.5) * self.pixel_size;
        let yoffset = (py + 0.5) * self.pixel_size;

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

        let samples_per_pixel = self.supersampling * self.supersampling;

        let mut rng = rand::thread_rng();

        for y in 0..self.vsize {
            for x in 0..self.hsize {
                for _ in 0..=samples_per_pixel {
                    let offset: f32 = rng.gen();
                    let ray = self.ray_for_pixel((x as f32) + offset, (y as f32) + offset);
                    let color = world.color_at(&ray, MAX_REFLECTIONS);
                    canvas.add(y as usize, x as usize, color);
                }
            }
        }

        canvas.divide(samples_per_pixel as f32);

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
    use crate::linear::{Matrix4, Point};

    use std::f32::consts::PI;

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
        let c = Camera::new(200, 125, PI / 2.0, 1);
        assert_eq!(c.pixel_size, 0.01);

        let c2 = Camera::new(125, 200, PI / 2.0, 1);
        assert_eq!(c2.pixel_size, 0.01);
    }

    #[test]
    fn ray_through_center() {
        let c = Camera::new(201, 101, PI / 2.0, 1);
        let r = c.ray_for_pixel(100.0, 50.0);
        assert_eq!(r.origin, Point::origin());
        assert_eq!(r.direction, Vector::new(0.0, 0.0, -1.0));
    }

    #[test]
    fn ray_through_corner() {
        let c = Camera::new(201, 101, PI / 2.0, 1);
        let r = c.ray_for_pixel(0.0, 0.0);
        assert_eq!(r.origin, Point::origin());
        assert_eq!(r.direction, Vector::new(0.66519, 0.33259, -0.66851));
    }

    #[test]
    fn ray_with_transformed_camera() {
        let mut c = Camera::new(201, 101, PI / 2.0, 1);
        c.transform = Matrix4::rotation_y(PI / 4.0) * Matrix4::translation(0.0, -2.0, 5.0);
        let r = c.ray_for_pixel(100.0, 50.0);
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
        let mut c = Camera::new(11, 11, PI / 2.0, 1);
        c.transform = view_transform(from, to, up);
        let image = c.render(&w);
        assert!(image.at(5, 5) != Color::black());
    }
}
