//! # biquad
//!
//! A library for creating second order IIR filters for signal processing based on
//! [Biquads](https://en.wikipedia.org/wiki/Digital_biquad_filter). Both a
//! Direct Form 1 (DF1) and Direct Form 2 Transposed (DF2T) implementation is
//! available, where the DF1 is better used when the filter needs retuning
//! online, as it has the property to introduce minimal artifacts under retuning,
//! while the DF2T is best used for static filters as it has the least
//! computational complexity and best numerical stability.
//!
//! # Usage example
//!
//! ```
//! fn main() {
//!     use biquad::*;
//!
//!     // Cutoff and sampling frequencies
//!     let f0 = 10.hz();
//!     let fs = 1.khz();
//!
//!     // Create coefficients for the biquads
//!     let coeffs = Coefficients::new(Type::LowPass, fs, f0, Q_BUTTERWORTH).unwrap();
//!
//!     // Create two different biquads
//!     let mut biquad1 = DirectForm1::new(coeffs);
//!     let mut biquad2 = DirectForm2Transposed::new(coeffs);
//!
//!     let input_vec = vec![0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
//!     let mut output_vec1 = Vec::new();
//!     let mut output_vec2 = Vec::new();
//!
//!     // Run for all the inputs
//!     for elem in input_vec {
//!         output_vec1.push(biquad1.run(elem));
//!         output_vec2.push(biquad2.run(elem));
//!     }
//! }
//! ```

#![no_std]

pub mod coefficients;
pub mod frequency;

pub use coefficients::*;
pub use frequency::*;

/// The required functions of a biquad implementation
pub trait Biquad {
    /// A single iteration of a biquad, applying the filtering on the input
    fn run(&mut self, input: f32) -> f32;

    /// Updating of coefficients
    fn update_coefficients(&mut self, new_coefficients: Coefficients);
}

/// Internal states and coefficients of the Direct Form 1 form
#[derive(Clone, Debug)]
pub struct DirectForm1 {
    y1: f32,
    y2: f32,
    x1: f32,
    x2: f32,
    coeffs: Coefficients,
}

/// Internal states and coefficients of the Direct Form 2 Transposed form
#[derive(Clone, Debug)]
pub struct DirectForm2Transposed {
    pub s1: f32,
    pub s2: f32,
    coeffs: Coefficients,
}

impl DirectForm1 {
    /// Creates a Direct Form 1 biquad from a set of filter coefficients
    pub fn new(coefficients: Coefficients) -> DirectForm1 {
        DirectForm1 {
            y1: 0.0,
            y2: 0.0,
            x1: 0.0,
            x2: 0.0,
            coeffs: coefficients,
        }
    }
}

impl Biquad for DirectForm1 {
    fn run(&mut self, input: f32) -> f32 {
        let out = self.coeffs.b0 * input + self.coeffs.b1 * self.x1 + self.coeffs.b2 * self.x2
            - self.coeffs.a1 * self.y1
            - self.coeffs.a2 * self.y2;

        self.x2 = self.x1;
        self.x1 = input;
        self.y2 = self.y1;
        self.y1 = out;

        out
    }

    fn update_coefficients(&mut self, new_coefficients: Coefficients) {
        self.coeffs = new_coefficients;
    }
}

impl DirectForm2Transposed {
    /// Creates a Direct Form 2 Transposed biquad from a set of filter coefficients
    pub fn new(coefficients: Coefficients) -> DirectForm2Transposed {
        DirectForm2Transposed {
            s1: 0.0,
            s2: 0.0,
            coeffs: coefficients,
        }
    }
}

impl Biquad for DirectForm2Transposed {
    fn run(&mut self, input: f32) -> f32 {
        let out = self.s1 + self.coeffs.b0 * input;
        self.s1 = self.s2 + self.coeffs.b1 * input - self.coeffs.a1 * out;
        self.s2 = self.coeffs.b2 * input - self.coeffs.a2 * out;

        out
    }

    fn update_coefficients(&mut self, new_coefficients: Coefficients) {
        self.coeffs = new_coefficients;
    }
}

#[cfg(test)]
#[macro_use]
extern crate std;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_frequency() {
        let f1 = 10.hz();
        let f2 = 10.khz();
        let f3 = 10.mhz();
        let f4 = 10.dt();

        assert_eq!(f1, Hertz::new(10.));
        assert_eq!(f2, Hertz::new(10000.));
        assert_eq!(f3, Hertz::new(10000000.));
        assert_eq!(f4, Hertz::new(0.1));

        assert!(f1 < f2);
        assert!(f3 > f2);
        assert!(f1 == f1);
        assert!(f1 != f2);
    }

    #[test]
    fn test_coefficients_normal() {
        let f0 = 10.hz();
        let fs = 1.khz();

        let coeffs = Coefficients::new(Type::LowPass, fs, f0, Q_BUTTERWORTH);

        match coeffs {
            Ok(_) => {}
            Err(s) => {
                panic!("Coefficients creation failed! Error: {}", s);
            }
        }
    }

    #[test]
    #[should_panic]
    fn test_coefficients_fail_flipped_frequencies() {
        let f0 = 10.hz();
        let fs = 1.khz();

        let coeffs = Coefficients::new(Type::LowPass, f0, fs, Q_BUTTERWORTH);

        match coeffs {
            Ok(_) => {}
            Err(s) => {
                panic!("Coefficients creation failed! Error: {}", s);
            }
        }
    }

    #[test]
    #[should_panic]
    fn test_coefficients_fail_negative_q() {
        let f0 = 10.hz();
        let fs = 1.khz();

        let coeffs = Coefficients::new(Type::LowPass, fs, f0, -1.0);

        match coeffs {
            Ok(_) => {}
            Err(s) => {
                panic!("Coefficients creation failed! Error: {}", s);
            }
        }
    }

    #[test]
    fn test_biquad_zeros() {
        use tests::std::vec::Vec;
        //use std::vec::Vec::*;

        let f0 = 10.hz();
        let fs = 1.khz();

        let coeffs = Coefficients::new(Type::LowPass, fs, f0, Q_BUTTERWORTH).unwrap();

        let mut biquad1 = DirectForm1::new(coeffs);
        let mut biquad2 = DirectForm2Transposed::new(coeffs);

        let input_vec = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut output_vec1 = Vec::new();
        let mut output_vec2 = Vec::new();

        for elem in input_vec {
            output_vec1.push(biquad1.run(elem));
            output_vec2.push(biquad2.run(elem));
        }
    }
}
