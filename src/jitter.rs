use rand::Rng;
use rand::rngs::ThreadRng;

#[derive(Debug)]
pub struct RandomJitter {
    rng: ThreadRng
}

impl RandomJitter {
    pub fn new() -> RandomJitter {
        RandomJitter {
            rng: rand::thread_rng()
        }
    }
}

impl Iterator for RandomJitter {
    type Item = f32;

    fn next(&mut self) -> Option<f32> {
        Some(self.rng.gen())
    }
}

pub struct ConstantJitter {
    value: f32
}

impl ConstantJitter {
    pub fn new(value: f32) -> ConstantJitter {
        ConstantJitter { value }
    }
}

impl Iterator for ConstantJitter {
    type Item = f32;

    fn next(&mut self) -> Option<f32> {
        Some(self.value)
    }
}