use indextree::{Arena, NodeId};
use crate::shapes::Shape;

pub struct Registry {
    arena: Arena<Box<Shape>>
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

    pub fn get(&self, id: NodeId) -> Option<&Box<Shape>> {
        let node = self.arena.get(id)?;
        Some(&node.data)
    }

    pub fn add(&mut self, group_id: NodeId, child_id: NodeId) {
        group_id.append(child_id, &mut self.arena).unwrap();
    }

    pub fn parent(&self, child: &Box<Shape>) -> Option<&Box<Shape>> {
        let child_id = child.tag().unwrap();
        let mut ancestors = child_id.ancestors(&self.arena);
        ancestors.next().unwrap(); // skip self
        let parent_id = ancestors.next()?;
        let parent = self.arena.get(parent_id)?;
        Some(&parent.data)
    }

    pub fn children(&self, group: &Box<Shape>) -> Vec<&Box<Shape>> {
        let parent_id = group.tag().unwrap();
        parent_id.children(&self.arena).map(|x| &self.arena.get(x).unwrap().data).collect()
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
        let child = reg.get(child_id).unwrap();
        assert_eq!(reg.parent(child).unwrap().id(), parent_id);
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
        let parent = reg.get(group_id).unwrap();
        assert_eq!(reg.children(parent).first().unwrap().id(), sphere_id);
    }
}