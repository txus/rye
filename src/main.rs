extern crate yaml_rust;
extern crate structopt;
use std::path::PathBuf;
use structopt::StructOpt;
use std::time::Instant;

mod camera;
mod canvas;
mod color;
mod light;
mod linear;
mod materials;
mod output;
mod patterns;
mod rays;
mod shapes;
mod world;
mod parser;
mod obj_parser;
mod registry;
mod jitter;

/// A rye (tracer)
#[derive(StructOpt, Debug)]
#[structopt(name = "rye")]
struct Opt {
    /// path to the scene YAML file
    #[structopt(short = "s", long = "scene", parse(from_os_str))]
    scene: PathBuf,

    /// desired path to the output PPM file
    #[structopt(short = "o", long = "output", parse(from_os_str))]
    output: PathBuf,

    /// width in pixels
    #[structopt(short = "w", long = "width")]
    width: u32,
 
    /// height in pixels
    #[structopt(short = "h", long = "height")]
    height: u32,

    /// samples per pixel = supersampling^2
    #[structopt(short = "s", long = "supersampling")]
    supersampling: usize
}

fn main() {
    let opt = Opt::from_args();
    let (world, camera) = parser::read_filename(&opt.scene.to_str().unwrap(), opt.width, opt.height, opt.supersampling).unwrap();
    let instant = Instant::now();
    let canvas = camera.render(&world);
    output::render(&canvas, output::ppm::PPM {}, opt.output.to_str().unwrap());
    println!("Rendered {} to {} in {:?}", opt.scene.to_str().unwrap(), opt.output.to_str().unwrap(), instant.elapsed());
}
