# rye (tracer)

My attempt at learning raytracing (and Rust!) with the fantastic [Ray Tracer Challenge](http://raytracerchallenge.com) book.

As an example, this scene (at 1200 * 800 resolution) took 3.3 min to render on my Macbook Pro:

![Example](scenes/basic.png?raw=true "Example scene")

The scene description is at [scenes/basic.yml](scenes/basic.yml).

## Running tests

```
cargo test
```

## Rendering a scene

To render a scene from the example [scenes](scenes) folder:

```
cargo build --release
./target/release/rye  --height 150 --width 200 --scene scenes/basic.yml --output render.ppm --supersampling 4

```
