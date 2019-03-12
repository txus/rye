use crate::rays::{Ray, Intersection, Sphere, Precomputation};
use crate::materials::Material;
use crate::light::{PointLight};
use crate::primitives::{Point, Color, Matrix4, Matrix, Vector};

pub struct World {
    objects: Vec<Sphere>,
    light_source: PointLight
}

impl World {
    pub fn default() -> Self {
        let light = PointLight::new(Point::new(-10.0, 10.0, -10.0), Color::white());
        let mut s1 = Sphere::unit();
        s1.material = Material {
            color: Color::new(0.8, 1.0, 0.6),
            diffuse: 0.7,
            specular: 0.2,
            ..s1.material
        };
        let mut s2 = Sphere::unit();
        s2.set_transform(Matrix4::scaling(0.5, 0.5, 0.5));

        World {
            objects: vec!(s1, s2),
            light_source: light
        }
    }

    pub fn intersect(&self, ray: &Ray) -> Vec<Intersection> {
        let mut out : Vec<Intersection> = vec!();
        for object in &self.objects {
            out.append(&mut object.intersect(&ray))
        }
        out.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap_or(std::cmp::Ordering::Equal));
        out
    }

    pub fn shade(&self, c: &Precomputation) -> Color {
        c.object.material.lighting(&self.light_source, &c.point, &c.eye, &c.normal)
    }

    pub fn color_at(&self, r: &Ray) -> Color {
        if let Some(i) = self.intersect(&r).first() {
            self.shade(&i.precompute(&r))
        } else {
            Color::black()
        }
    }
}

pub struct Camera {
    hsize: u32,
    vsize: u32,
    fov: f32,
    transform: Matrix4
}

impl Camera {
    pub fn new(hsize: u32, vsize: u32, fov: f32) -> Camera {
        Camera { hsize, vsize, fov, transform: Matrix4::id() }
    }

    pub fn pixel_size(&self) -> f32 {
        let half_view = (self.fov / 2.0).tan();
        let aspect_ratio = self.hsize as f32 / self.vsize as f32;
        let half_width: f32;
        let half_height: f32;
        if aspect_ratio >= 1.0 {
            half_width = half_view;
            half_height = half_view / aspect_ratio;
        } else {
            half_width = half_view * aspect_ratio;
            half_height = half_view;
        }
        (half_width * 2.0) / self.hsize as f32
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
        [0.0, 0.0, 0.0, 1.0]
    ]);
    orientation * Matrix4::translation(-from.x, -from.y, -from.z)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primitives::{Point, Color};
    use crate::materials::Material;
    use crate::light::PointLight;

    use std::f32::consts::PI;

    #[test]
    fn default_world() {
        let light = PointLight::new(Point::new(-10.0, 10.0, -10.0), Color::white());
        let mut s1 = Sphere::unit();
        s1.material = Material {
            color: Color::new(0.8, 1.0, 0.6),
            diffuse: 0.7,
            specular: 0.2,
            ..s1.material
        };
        let mut s2 = Sphere::unit();
        s2.set_transform(Matrix4::scaling(0.5, 0.5, 0.5));

        let w = World::default();
        assert_eq!(w.light_source, light);
        assert_eq!(w.objects.iter().filter(|x| **x == s1).count(), 1);
        assert_eq!(w.objects.iter().filter(|x| **x == s2).count(), 1);
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
        let shape = w.objects.first().unwrap();
        let i = Intersection { t: 4.0, object: shape };
        let comps = i.precompute(&r);
        assert_eq!(w.shade(&comps), Color::new(0.38066, 0.47583, 0.2855));
    }

    #[test]
    fn shading_intersection_from_inside() {
        let mut w = World::default();
        w.light_source = PointLight::new(Point::new(0.0, 0.25, 0.0), Color::new(1.0, 1.0, 1.0));
        let r = Ray::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 1.0));
        let shape = w.objects.last().unwrap();
        let i = Intersection { t: 0.5, object: shape };
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
        let from = Point::new(0.0, 0.0, 0.0);
        let to = Point::new(0.0, 0.0, -10.0);
        let up = Vector::new(0.0, 1.0, 0.0);
        assert_eq!(view_transform(from, to, up), Matrix4::id());
    }

    #[test]
    fn view_transform_looking_pos_z_direction() {
        let from = Point::new(0.0, 0.0, 0.0);
        let to = Point::new(0.0, 0.0, 1.0);
        let up = Vector::new(0.0, 1.0, 0.0);
        assert_eq!(view_transform(from, to, up), Matrix4::scaling(-1.0, 1.0, -1.0));
    }

    #[test]
    fn view_transform_moves_world() {
        let from = Point::new(0.0, 0.0, 8.0);
        let to = Point::new(0.0, 0.0, 0.0);
        let up = Vector::new(0.0, 1.0, 0.0);
        assert_eq!(view_transform(from, to, up), Matrix4::translation(0.0, 0.0, -8.0));
    }

    #[test]
    fn view_transform_arbitrary() {
        let from = Point::new(1.0, 3.0, 2.0);
        let to = Point::new(4.0, -2.0, 8.0);
        let up = Vector::new(1.0, 1.0, 0.0);
        assert_eq!(
            view_transform(from, to, up),
            Matrix4::new([
                [-0.50709 , 0.50709 ,  0.67612 , -2.36643 ],
                [ 0.76772 , 0.60609 ,  0.12122 , -2.82843 ],
                [-0.35857 , 0.59761 , -0.71714 ,  0.00000 ],
                [ 0.00000 , 0.00000 ,  0.00000 ,  1.00000 ]
            ])
        );
    }

    #[test]
    fn camera_pixel_size() {
        let c = Camera::new(200, 125, PI/2.0);
        assert_eq!(c.pixel_size(), 0.01);

        let c2 = Camera::new(125, 200, PI/2.0);
        assert_eq!(c2.pixel_size(), 0.01);
    }
}