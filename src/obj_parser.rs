use crate::registry::Registry;
use crate::shapes::{Group, Triangle, SmoothTriangle, Shape};
use crate::linear::{Point, Vector};

use indextree::NodeId;

use std::fs;

#[derive(Debug)]
pub enum Error {
    File(String),
}

pub struct ParseResult {
    pub ignored_count: usize,
    pub triangles_count: usize,
    pub vertices: Vec<Point>,
    pub root: NodeId,
    pub group_ids: Vec<NodeId>,
}

fn fan_triangulation(faces: &[(Point, Option<Point>, Option<Vector>)]) -> Vec<Box<Shape>> {
    let mut triangles: Vec<Box<Shape>> = vec![];

    for index in 1..(faces.len() - 1) {
        let (a, _, mna) = faces[0];
        let (b, _, mnb) = faces[index];
        let (c, _, mnc) = faces[index + 1];
        let tri: Box<Shape> = if let (Some(na), Some(nb), Some(nc)) = (mna, mnb, mnc) {
            Box::from(SmoothTriangle::new(a, b, c, na, nb, nc))
        } else {
            Box::from(Triangle::new(a, b, c))
        };
        triangles.push(tri);
    }
    triangles
}

pub fn read_string(reg: &mut Registry, s: &str) -> Result<ParseResult, Error> {
    let mut ignored_count: usize = 0;
    let mut triangles_count: usize = 0;
    let mut group_ids: Vec<NodeId> = vec![];
    let mut vertices: Vec<Point> = vec![];
    let mut normals: Vec<Vector> = vec![];

    let group = Box::from(Group::new());
    let root = reg.register(group);
    let mut current_group = root;

    for line in s.lines() {
        let mut chars = line.chars();
        match chars.next() {
            Some('v') => {
                if let Some('n') = chars.next() {
                    let numbers = chars.collect::<String>().trim().split(' ').flat_map(str::parse::<f32>).collect::<Vec<_>>();
                    let (x, y, z) = (numbers[0], numbers[1], numbers[2]);
                    let v = Vector::new(x, y, z);
                    normals.push(v);
                } else {
                    let numbers = chars.collect::<String>().trim().split(' ').flat_map(str::parse::<f32>).collect::<Vec<_>>();
                    let (x, y, z) = (numbers[0], numbers[1], numbers[2]);
                    let v = Point::new(x, y, z);
                    vertices.push(v);
                }
            },
            Some('f') => {
                let faces: Vec<(Point, Option<Point>, Option<Vector>)> = chars.collect::<String>().trim().split(' ').map(|s| {
                    let parts = s.split('/').map(|part| str::parse::<usize>(part).unwrap_or(0)).collect::<Vec<_>>();
                    let (a, _b, c) = (parts.get(0).unwrap(), parts.get(1), parts.get(2));
                    let texture = None; // from b
                    let normal = c.map(|n| normals[n-1]);
                    (vertices[a - 1], texture, normal)
                }).collect::<Vec<_>>();
                let triangles = fan_triangulation(&faces);
                for t in triangles {
                    let triangle = Box::from(t);
                    let id = reg.register(triangle);
                    reg.add(current_group, id);
                    triangles_count += 1;
                }
            },
            Some('g') => {
                chars.next(); //skip space
                let _group_name = chars.collect::<String>().trim();
                let group = Box::from(Group::new());
                let gid = reg.register(group);
                reg.add(root, gid);
                current_group = gid;
                group_ids.push(gid);
            },
            _ => {
                ignored_count += 1;
            }
        }
    }

    Ok(ParseResult {
        ignored_count,
        group_ids,
        triangles_count,
        vertices,
        root,
    })
}

pub fn read_filename(reg: &mut Registry, filename: &str) -> Result<ParseResult, Error> {
    let s = fs::read_to_string(filename).or(Err(Error::File(format!("Can't read file {:?}", filename))))?;
    let result = read_string(reg, &s)?;
    Ok(result)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ignoring_unrecognized_lines() {
        let mut registry = Registry::new();
        let s = String::from("FROTFO\nurstpnp");
        let result = read_string(&mut registry, &s).unwrap();
        assert_eq!(result.ignored_count, 2);
    }

    #[test]
    fn vertex_records() {
        let mut registry = Registry::new();
        let s = String::from("
v -1 1 0
v -1.0000 0.5000 0.0000
v 1 0 0
v 1 1 0");
        let result = read_string(&mut registry, &s).unwrap();
        assert_eq!(result.vertices[0], Point::new(-1.0, 1.0, 0.0));
        assert_eq!(result.vertices[1], Point::new(-1.0, 0.5, 0.0));
        assert_eq!(result.vertices[2], Point::new(1.0, 0.0, 0.0));
        assert_eq!(result.vertices[3], Point::new(1.0, 1.0, 0.0));
    }

    #[test]
    fn triangle_faces() {
        let mut registry = Registry::new();
        let s = String::from("
v -1 1 0
v -1 0 0
v 1 0 0
v 1 1 0

f 1 2 3
f 1 3 4");
        let result = read_string(&mut registry, &s).unwrap();
        let id = result.root;
        let children = registry.children(id);
        assert_eq!(children.len(), 2);
    }

    #[test]
    fn triangulating_polygons() {
        let mut registry = Registry::new();
        let s = String::from("
v -1 1 0
v -1 0 0
v 1 0 0
v 1 1 0
v 0 2 0
f 1 2 3 4 5");
        let result = read_string(&mut registry, &s).unwrap();
        let id = result.root;
        let children = registry.children(id);
        assert_eq!(children.len(), 3);
    }

    #[test]
    fn triangles_in_groups() {
        let mut registry = Registry::new();
        let s = String::from("
v -1 1 0
v -1 0 0
v 1 0 0
v 1 1 0
g FirstGroup
f 1 2 3
g SecondGroup
f 1 3 4");
        let result = read_string(&mut registry, &s).unwrap();
        assert_eq!(result.group_ids.len(), 2);
        assert_eq!(result.triangles_count, 2);
        let id = result.root;
        let children = registry.children(id);
        assert_eq!(children.len(), 2); // groups
    }

    #[test]
    fn vertex_normal_data() {
        let mut registry = Registry::new();
        let s = String::from("
v 0 1 0
v -1 0 0
v 1 0 0
vn -1 0 0
vn 1 0 0
vn 0 1 0
f 1//3 2//1 3//2
f 1/0/3 2/102/1 3/14/2");
        let result = read_string(&mut registry, &s).unwrap();
        let id = result.root;
        let children = registry.children(id);
        assert_eq!(children.len(), 2);
    }
}