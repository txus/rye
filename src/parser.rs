use yaml_rust::{YamlLoader, Yaml};
use std::fs;
use indextree::NodeId;
use crate::light::PointLight;
use crate::world::World;
use crate::shapes::{Shape, Cube, Sphere, Cylinder, Cone, Plane, Group, precompute_bounds, CSG, CSGOperation};
use crate::linear::{Point, Vector, Matrix4, Matrix};
use crate::materials::Material;
use crate::patterns::{Pattern, StripePattern, GradientPattern, RingPattern, CheckerPattern};
use crate::color::Color;
use crate::camera::{Camera, view_transform};
use crate::registry::Registry;
use crate::obj_parser;

use std::rc::Rc;
use std::cell::RefCell;

#[derive(Debug)]
pub struct ParseError(String);

#[derive(Debug)]
pub enum Error {
    Parse(String),
    File(String),
    EmptyDocument,
}

fn parse_f32(doc: &Yaml, context: &Yaml) -> Result<f32, Error> {
    match doc.as_f64().map(|x| x as f32).or(doc.as_i64().map(|x| x as f32)) {
        Some(f) => Ok(f),
        _ => Err(Error::Parse(format!("Error parsing number: {:?} -- context: {:?}", doc, context)))
    }
}

fn parse_tuple(doc: &Yaml) -> Result<(f32, f32, f32), Error> {
    let a = parse_f32(&doc[0], doc)?;
    let b = parse_f32(&doc[1], doc)?;
    let c = parse_f32(&doc[2], doc)?;
    Ok((a, b, c))
}

fn parse_point(doc: &Yaml) -> Result<Point, Error> {
    parse_tuple(doc).map(|(a, b, c)| Point::new(a, b, c))
}

fn parse_vector(doc: &Yaml) -> Result<Vector, Error> {
    parse_tuple(doc).map(|(a, b, c)| Vector::new(a, b, c))
}

fn parse_color(doc: &Yaml) -> Result<Color, Error> {
    parse_tuple(doc).map(|(a, b, c)| Color::new(a, b, c))
}

fn parse_light(doc: &Yaml) -> Result<PointLight, Error> {
    let position = parse_point(&doc["at"])?;
    let color = parse_color(&doc["intensity"])?;
    match doc["type"].as_str() {
        Some("PointLight") => Ok(PointLight {
            position: position,
            intensity: color,
        }),
        _ => Err(unknown("Light", doc, "type", "PointLight"))
    }
}

fn unknown(kind: &'static str, document: &Yaml, key: &'static str, options: &'static str) -> Error {
    Error::Parse(format!("{} not recognized: {:?} (available: {})-- context: {:?}", kind, document[key], options, document))

}

fn parse_pattern(doc: &Yaml) -> Result<Box<Pattern>, Error> {
    let a = parse_color(&doc["a"])?;
    let b = parse_color(&doc["b"])?;
    let (translation, rotation, scale) = parse_transform(doc);
    let transform = rotation * scale * translation;
    match doc["type"].as_str() {
        Some("Stripe") => Ok(Box::from(StripePattern { a, b, transform })),
        Some("Gradient") => Ok(Box::from(GradientPattern { a, b, transform })),
        Some("Ring") => Ok(Box::from(RingPattern { a, b, transform })),
        Some("Checker") => Ok(Box::from(CheckerPattern { a, b, transform })),
        _ => Err(unknown("Pattern", doc, "type", "Stripe, Gradient, Ring, Checker"))
    }
}

fn parse_material(doc: &Yaml) -> Result<Material, Error> {
    let color = parse_color(&doc["color"])?;
    let ambient = parse_f32(&doc["ambient"], doc);
    let diffuse = parse_f32(&doc["diffuse"], doc);
    let specular = parse_f32(&doc["specular"], doc);
    let shininess = parse_f32(&doc["shininess"], doc);
    let reflective = parse_f32(&doc["reflective"], doc);
    let transparency = parse_f32(&doc["transparency"], doc);
    let refractive_index = parse_f32(&doc["refractive_index"], doc);
    let pattern = parse_pattern(&doc["pattern"]);
    let mut m = Material::default();
    m.color = color;
    if let Ok(x) = ambient { m.ambient = x };
    if let Ok(x) = diffuse { m.diffuse = x };
    if let Ok(x) = specular { m.specular = x };
    if let Ok(x) = shininess { m.shininess = x };
    if let Ok(x) = reflective { m.reflective = x };
    if let Ok(x) = transparency { m.transparency = x };
    if let Ok(x) = refractive_index { m.refractive_index = x };
    if let Ok(x) = pattern { m.pattern = Some(x) };
    Ok(m)
}

fn parse_transform(doc: &Yaml) -> (Matrix4, Matrix4, Matrix4) {
    let translation = parse_tuple(&doc["at"]).map(|(x, y, z)| Matrix4::translation(x, y, z)).unwrap_or(Matrix4::id());
    let rotation = parse_tuple(&doc["rotate"]).map(|(x, y, z)| Matrix4::rotation_x(x.to_radians()) * Matrix4::rotation_y(y.to_radians()) * Matrix4::rotation_z(z.to_radians())).unwrap_or(Matrix4::id());
    let scale = parse_tuple(&doc["scale"]).map(|(x, y, z)| Matrix4::scaling(x, y, z)).unwrap_or(Matrix4::id());

    (translation, rotation, scale)
}

fn parse_object(registry: Rc<RefCell<Registry>>, doc: &Yaml) -> Result<NodeId, Error> {
    let (translation, rotation, scale) = parse_transform(&doc);
    let transform = translation * scale * rotation;
    let material = parse_material(&doc["material"]);
    let casts_shadows = doc["casts_shadows"].as_bool().unwrap_or(true);
    match doc["shape"].as_str() {
        Some("Cube") => {
            let mut c = Cube::new();
            c.set_transform(transform);
            if let Ok(m) = material { c.set_material(m) }
            c.casts_shadows = casts_shadows;
            Ok({
                let mut reg = registry.borrow_mut();
                reg.register(Box::from(c))
            })
        },
        Some("Sphere") => {
            let mut s = Sphere::new();
            s.set_transform(transform);
            if let Ok(m) = material { s.set_material(m) }
            s.casts_shadows = casts_shadows;
            Ok({
                let mut reg = registry.borrow_mut();
                reg.register(Box::from(s))
            })
        },
        Some("Plane") => {
            let mut p = Plane::new();
            p.set_transform(transform);
            if let Ok(m) = material { p.set_material(m) }
            p.casts_shadows = casts_shadows;
            Ok({
                let mut reg = registry.borrow_mut();
                reg.register(Box::from(p))
            })
        },
        Some("Cylinder") => {
            let mut c = Cylinder::new();
            c.set_transform(transform);
            if let Ok(m) = material { c.set_material(m) }
            if let Ok(min) = parse_f32(&doc["minimum"], doc) {
                c.minimum = min;
            }
            if let Ok(max) = parse_f32(&doc["maximum"], doc) {
                c.maximum = max;
            }
            if let Ok(max) = parse_f32(&doc["maximum"], doc) {
                c.maximum = max;
            }
            let closed = match doc["closed"].as_str() {
                Some("true") => true,
                _ => false
            };
            c.closed = closed;
            c.casts_shadows = casts_shadows;
            Ok({
                let mut reg = registry.borrow_mut();
                reg.register(Box::from(c))
            })
        },
        Some("Cone") => {
            let mut c = Cone::new();
            c.set_transform(transform);
            if let Ok(m) = material { c.set_material(m) }
            if let Ok(min) = parse_f32(&doc["minimum"], doc) {
                c.minimum = min;
            }
            if let Ok(max) = parse_f32(&doc["maximum"], doc) {
                c.maximum = max;
            }
            if let Ok(max) = parse_f32(&doc["maximum"], doc) {
                c.maximum = max;
            }
            let closed = match doc["closed"].as_str() {
                Some("true") => true,
                _ => false
            };
            c.closed = closed;
            c.casts_shadows = casts_shadows;
            Ok({
                let mut reg = registry.borrow_mut();
                reg.register(Box::from(c))
            })
        },
        Some("CSG") => {
            let left = parse_object(registry.clone(), &doc["left"]).unwrap();
            let right = parse_object(registry.clone(), &doc["right"]).unwrap();
            let combined_bounds = {
                let reg = registry.borrow();
                let l = reg.get(left);
                let r = reg.get(right);
                precompute_bounds(vec![l, r])
            };
            let operation: CSGOperation = match &doc["operation"].as_str() {
                Some("Union") => CSGOperation::Union,
                Some("Intersection") => CSGOperation::Intersection,
                Some("Difference") => CSGOperation::Difference,
                _ => panic!("unknown CSG operation")
            };
            Ok({
                let mut reg = registry.borrow_mut();
                let mut csg = Box::from(CSG::new(operation, left, right));
                csg.set_transform(transform);
                csg.set_bounds(combined_bounds);
                let id = reg.register(csg);
                reg.add(id, left);
                reg.add(id, right);
                id
            })
        },
        Some("Group") => {
            let gid = {
                let gid = {
                    let mut reg = registry.borrow_mut();
                    let mut g = Box::from(Group::new());
                    g.set_transform(transform);
                    if let Ok(m) = material { g.set_material(m) }
                    reg.register(g)
                };

                let obj_array = match &doc["children"] {
                    Yaml::Array(a) => Ok(a),
                    _ => Err(Error::Parse("'children' is not an array".to_owned()))
                }?;

                for o in obj_array.iter() {
                    let id = parse_object(registry.clone(), &o).unwrap();
                    let mut reg = registry.borrow_mut();
                    reg.add(gid, id);
                }
                gid
            };
            let bounds = {
                let reg = registry.borrow();
                precompute_bounds(reg.children(gid))
            };
            {
                let mut reg = registry.borrow_mut();
                let g = reg.get_mut(gid);
                g.set_bounds(bounds);
            }
            Ok(gid)
        },
        Some("Obj") => {
            let filename = doc["filename"].as_str().expect("no filename in OBJ clause");
            let results = {
                let mut reg = registry.borrow_mut();
                let results = obj_parser::read_filename(&mut reg, &filename).expect("can't parse OBJ file");
                let root_group_id = results.root;
                let g = reg.get_mut(root_group_id);
                g.set_transform(transform);
                if let Ok(m) = material { g.set_material(m) }
                results
            };

            for group_id in results.group_ids {
                let bounds = {
                    let reg = registry.borrow();
                    precompute_bounds(reg.children(group_id))
                };
                {
                    let mut reg = registry.borrow_mut();
                    let group = reg.get_mut(group_id);
                    group.set_bounds(bounds);
                }
            }

            let root_group_bounds = {
                let reg = registry.borrow();
                precompute_bounds(reg.children(results.root))
            };
            {
                let mut reg = registry.borrow_mut();
                let g = reg.get_mut(results.root);
                g.set_bounds(root_group_bounds);
            }
            Ok(results.root)
        }
        _ => Err(unknown("Shape", doc, "type", "Cube, Sphere, Plane, Cylinder, Cone, Group, Obj, CSG"))
    }
}

pub fn read_string(s: &str) -> Result<(World, Point, Point, Vector, f32), Error> {
    let docs = YamlLoader::load_from_str(s).or(Err(Error::EmptyDocument))?;
    let doc = &docs[0];
    let light = parse_light(&doc["light"])?;
    let mut w = World::empty();

    let obj_array = match &doc["objects"] {
        Yaml::Array(a) => Ok(a),
        _ => Err(Error::Parse("'objects' is not an array".to_owned()))
    }?;

    for o in obj_array.iter() {
        parse_object(w.registry.clone(), &o).unwrap();
    }

    w.light_source = light;

    let at = parse_point(&doc["camera"]["at"])?;
    let look_at = parse_point(&doc["camera"]["look_at"])?;
    let up = parse_vector(&doc["camera"]["look_at"])?;
    let fov = parse_f32(&doc["camera"]["fov"], &doc["camera"])?;

    Ok((w, at, look_at, up, fov.to_radians()))
}

pub fn read_filename(filename: &str, width: u32, height: u32, supersampling: usize) -> Result<(World, Camera), Error> {
    let s = fs::read_to_string(filename).or(Err(Error::File(format!("Can't read file {:?}", filename))))?;
    let (world, at, look_at, up, fov) = read_string(&s)?;
    let mut c = Camera::new(width, height, fov, supersampling);
    c.transform = view_transform(at, look_at, up);
    Ok((world, c))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn world_parsing() {
        let s =
"
camera:
    at: [0, 2.5, -3.5]
    look_at: [0, 0, 3]
    up: [0, 1, 0]
    fov: 50
light:
    type: PointLight
    at: [1, 6, 5]
    intensity: [1, 1, 1]

objects:
    - shape: Cube
      at: [1.0, 6.0, 5.0]
      rotate: [0, 45, 0]
      material:
        color: [1, 0, 1]
        pattern:
            type: Checker
            a: [0.5, 0, 1]
            b: [0, 0, 0]
            at: [3, 1, 2]
            rotation: [0.5, 0.1, 0.3]
            scale: [0.5, 0.5, 0.5]
    - shape: Plane
      at: [0, 0, 0]
    - shape: Cone
      at: [0, 0, 0]
      closed: true
      minimum: 10
    - shape: Cylinder
      at: [0, 0, 0]
      minimum: 10
      maximum: 20
    - shape: Sphere
      at: [1, 0, 0]
      scale: [1, 2, 6]
";
        let (world, _camera_at, _camera_look_at, _camera_up, _camera_fov) = read_string(&s).unwrap();
        assert_eq!(world.light_source.position, Point::new(1.0, 6.0, 5.0));
        assert_eq!(world.light_source.intensity, Color::white());
        assert_eq!(world.registry.borrow().all().len(), 5);
    }
}