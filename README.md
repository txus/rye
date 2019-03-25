# rye (tracer)

My attempt at learning raytracing (and Rust!) with the fantastic [Ray Tracer Challenge](http://raytracerchallenge.com) book.

```
cargo test
```

To render whatever scene is described in `src/main.rs`:

```
cargo build --release && ./target/release/rye
```

Or if you don't mind waiting longer:

```
cargo run
```

The output will be at `/tmp/refraction.ppm`.
