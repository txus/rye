use indextree::{Arena, NodeId};
use crate::shapes::Shape;
use std::num::NonZeroUsize;

pub struct Registry {
    arena: Arena<Box<Shape>>
}

pub fn id() -> NodeId {
    id_from(1)
}
pub fn id_from(id: usize) -> NodeId {
    NodeId::from_non_zero_usize(NonZeroUsize::new(id).unwrap())
}

impl Registry {
    pub fn new() -> Registry {
        Registry {
            arena: Arena::new()
        }
    }

    pub fn register(&mut self, mut shape: Box<Shape>) -> NodeId {
        let write = &mut shape as *mut Box<Shape>;
        let node_id = self.arena.new_node(shape);
        unsafe { (*write).set_tag(node_id) }
        node_id
    }

    pub fn get<'a>(&self, id: NodeId) -> &Box<Shape> {
        &self.arena.get(id).expect("can't find node in registry").data
    }

    pub fn add(&mut self, group_id: NodeId, child_id: NodeId) {
        group_id.append(child_id, &mut self.arena).unwrap();
    }

    pub fn parent(&self, child_id: NodeId) -> Option<&Box<Shape>> {
        let mut ancestors = child_id.ancestors(&self.arena);
        ancestors.next().unwrap(); //skip self
        let parent_id = ancestors.next()?;
        let parent = self.arena.get(parent_id)?;
        Some(&parent.data)
    }

    pub fn children(&self, group_id: NodeId) -> Vec<&Box<Shape>> {
        let ids = group_id.children(&self.arena);
        let mut shapes = vec![];
        for id in ids {
            shapes.push(&self.arena.get(id).unwrap().data);
        }
        shapes
    }

    pub fn all(&self) -> Vec<&Box<Shape>> {
        self.arena.iter().map(|x| &x.data).collect()
    }

    pub fn all_ids(&self) -> Vec<NodeId> {
        let count = self.arena.count();
        let mut out: Vec<NodeId> = vec![];
        for i in 1..=count { // 1-based indexing (dont ask)
            out.push(NodeId::from_non_zero_usize(NonZeroUsize::new(i).unwrap()))
        }
        out
    }

    pub fn get_mut(&mut self, id: NodeId) -> &mut Box<Shape> {
        &mut self.arena.get_mut(id).expect("can't find node by id").data
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shapes::{Group, Sphere};

    #[test]
    fn writing_and_reading() {
        let mut reg = Registry::new();
        let s: Box<Shape> = Box::from(Sphere::new());
        reg.register(s);
    }

    #[test]
    fn get_parent() {
        let mut reg = Registry::new();
        let parent = Group::new();
        let parent_id = parent.id();
        let s: Box<Shape> = Box::from(Sphere::new());
        let g: Box<Shape> = Box::from(parent);
        let child_id = reg.register(s);
        let group_id = reg.register(g);
        reg.add(group_id, child_id);
        assert_eq!(reg.parent(child_id).unwrap().id(), parent_id);
    }

    #[test]
    fn get_children() {
        let mut reg = Registry::new();
        let parent = Group::new();
        let sphere = Sphere::new();
        let sphere_id = sphere.id();
        let s: Box<Shape> = Box::from(sphere);
        let g: Box<Shape> = Box::from(parent);
        let child_id = reg.register(s);
        let group_id = reg.register(g);
        reg.add(group_id, child_id);
        let parent = reg.get(group_id);
        assert_eq!(reg.children(group_id).first().unwrap().id(), sphere_id);
    }
}