# `biquad`

[![Build Status](https://www.travis-ci.org/korken89/biquad-rs.svg?branch=master)](https://www.travis-ci.org/korken89/biquad-rs)

`biquad` is a `#![no_std]` library for creating first and second order IIR
filters for signal processing based on
[Biquads](https://en.wikipedia.org/wiki/Digital_biquad_filter). Both a
Direct Form 1 (DF1) and Direct Form 2 Transposed (DF2T) implementation is
available, where the DF1 is better used when the filter needs retuning
online, as it has the property to introduce minimal artifacts under retuning,
while the DF2T is best used for static filters as it has the least
computational complexity and best numerical stability.

This crate implements the biquads for `f32` and `f64`.

## Example

```rust
fn main() {
    use biquad::*;

    // Cutoff and sampling frequencies
    let f0 = 10.hz();
    let fs = 1.khz();

    // Create coefficients for the biquads
    let coeffs = Coefficients::<f32>::from_params(Type::LowPass, fs, f0, Q_BUTTERWORTH_F32).unwrap();

    // Create two different biquads
    let mut biquad1 = DirectForm1::<f32>::new(coeffs);
    let mut biquad2 = DirectForm2Transposed::<f32>::new(coeffs);

    let input_vec = vec![0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
    let mut output_vec1 = Vec::new();
    let mut output_vec2 = Vec::new();

    // Run for all the inputs
    for elem in input_vec {
        output_vec1.push(biquad1.run(elem));
        output_vec2.push(biquad2.run(elem));
    }
}
```

# [Documentation](https://docs.rs/biquad)

# License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or
  http://www.apache.org/licenses/LICENSE-2.0)

- MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## Contribution

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, as defined in the Apache-2.0 license, shall be
dual licensed as above, without any additional terms or conditions.
