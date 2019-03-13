use crate::rays::Sphere;

pub trait Object
where
    Self: std::fmt::Debug + PartialEq,
{
}

impl Object for Sphere {}
