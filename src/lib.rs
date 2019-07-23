//! # biquad
//!
//! `biquad` is a library for creating second order IIR filters for signal processing based on
//! [Biquads](https://en.wikipedia.org/wiki/Digital_biquad_filter). Both a
//! Direct Form 1 (DF1) and Direct Form 2 Transposed (DF2T) implementation is
//! available, where the DF1 is better used when the filter needs retuning
//! online, as it has the property to introduce minimal artifacts under retuning,
//! while the DF2T is best used for static filters as it has the least
//! computational complexity and best numerical stability.
//!
//! # Examples
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
//!     let coeffs = Coefficients::<f32>::from_params(Type::LowPass, fs, f0, Q_BUTTERWORTH_F32).unwrap();
//!
//!     // Create two different biquads
//!     let mut biquad1 = DirectForm1::<f32>::new(coeffs);
//!     let mut biquad2 = DirectForm2Transposed::<f32>::new(coeffs);
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
//!
//! # Errors
//!
//! `Coefficients::from_params(...)` can error if the cutoff frequency does not adhere to the
//! [Nyquist Frequency](https://en.wikipedia.org/wiki/Nyquist_frequency), or if the Q value is
//! negative.
//!
//! `Hertz::from_hz(...)` and `Hertz::from_dt(...)` will error if the frequency is negative.
//!
//! # Panics
//!
//! `x.hz()`, `x.khz()`, `x.mhz()`, `x.dt()` will panic for `f32`/`f64` if they are negative.
//!

#![no_std]

pub mod coefficients;
pub mod frequency;

pub use crate::coefficients::*;
pub use crate::frequency::*;

/// The required functions of a biquad implementation
pub trait Biquad<T> {
    /// A single iteration of a biquad, applying the filtering on the input
    fn run(&mut self, input: T) -> T;

    /// Updating of coefficients
    fn update_coefficients(&mut self, new_coefficients: Coefficients<T>);
}

/// Possible errors
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Errors {
    OutsideNyquist,
    NegativeQ,
    NegativeFrequency,
}

/// Internal states and coefficients of the Direct Form 1 form
#[derive(Copy, Clone, Debug)]
pub struct DirectForm1<T> {
    y1: T,
    y2: T,
    x1: T,
    x2: T,
    coeffs: Coefficients<T>,
}

/// Internal states and coefficients of the Direct Form 2 Transposed form
#[derive(Copy, Clone, Debug)]
pub struct DirectForm2Transposed<T> {
    pub s1: T,
    pub s2: T,
    coeffs: Coefficients<T>,
}

impl DirectForm1<f32> {
    /// Creates a Direct Form 1 biquad from a set of filter coefficients
    pub fn new(coefficients: Coefficients<f32>) -> Self {
        DirectForm1 {
            y1: 0.0_f32,
            y2: 0.0_f32,
            x1: 0.0_f32,
            x2: 0.0_f32,
            coeffs: coefficients,
        }
    }
}

impl Biquad<f32> for DirectForm1<f32> {
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

    fn update_coefficients(&mut self, new_coefficients: Coefficients<f32>) {
        self.coeffs = new_coefficients;
    }
}

impl DirectForm1<f64> {
    /// Creates a Direct Form 1 biquad from a set of filter coefficients
    pub fn new(coefficients: Coefficients<f64>) -> Self {
        DirectForm1 {
            y1: 0.0_f64,
            y2: 0.0_f64,
            x1: 0.0_f64,
            x2: 0.0_f64,
            coeffs: coefficients,
        }
    }
}

impl Biquad<f64> for DirectForm1<f64> {
    fn run(&mut self, input: f64) -> f64 {
        let out = self.coeffs.b0 * input + self.coeffs.b1 * self.x1 + self.coeffs.b2 * self.x2
            - self.coeffs.a1 * self.y1
            - self.coeffs.a2 * self.y2;

        self.x2 = self.x1;
        self.x1 = input;
        self.y2 = self.y1;
        self.y1 = out;

        out
    }

    fn update_coefficients(&mut self, new_coefficients: Coefficients<f64>) {
        self.coeffs = new_coefficients;
    }
}

impl DirectForm2Transposed<f32> {
    /// Creates a Direct Form 2 Transposed biquad from a set of filter coefficients
    pub fn new(coefficients: Coefficients<f32>) -> Self {
        DirectForm2Transposed {
            s1: 0.0_f32,
            s2: 0.0_f32,
            coeffs: coefficients,
        }
    }
}

impl Biquad<f32> for DirectForm2Transposed<f32> {
    fn run(&mut self, input: f32) -> f32 {
        let out = self.s1 + self.coeffs.b0 * input;
        self.s1 = self.s2 + self.coeffs.b1 * input - self.coeffs.a1 * out;
        self.s2 = self.coeffs.b2 * input - self.coeffs.a2 * out;

        out
    }

    fn update_coefficients(&mut self, new_coefficients: Coefficients<f32>) {
        self.coeffs = new_coefficients;
    }
}

impl DirectForm2Transposed<f64> {
    /// Creates a Direct Form 2 Transposed biquad from a set of filter coefficients
    pub fn new(coefficients: Coefficients<f64>) -> Self {
        DirectForm2Transposed {
            s1: 0.0_f64,
            s2: 0.0_f64,
            coeffs: coefficients,
        }
    }
}

impl Biquad<f64> for DirectForm2Transposed<f64> {
    fn run(&mut self, input: f64) -> f64 {
        let out = self.s1 + self.coeffs.b0 * input;
        self.s1 = self.s2 + self.coeffs.b1 * input - self.coeffs.a1 * out;
        self.s2 = self.coeffs.b2 * input - self.coeffs.a2 * out;

        out
    }

    fn update_coefficients(&mut self, new_coefficients: Coefficients<f64>) {
        self.coeffs = new_coefficients;
    }
}

#[cfg(test)]
#[macro_use]
extern crate std;

#[cfg(test)]
mod tests {
    use crate::*;

    #[test]
    fn test_frequency_f32() {
        let f1 = 10.hz();
        let f2 = 10.khz();
        let f3 = 10.mhz();
        let f4 = 10.dt();

        assert_eq!(f1, Hertz::<f32>::from_hz(10.).unwrap());
        assert_eq!(f2, Hertz::<f32>::from_hz(10000.).unwrap());
        assert_eq!(f3, Hertz::<f32>::from_hz(10000000.).unwrap());
        assert_eq!(f4, Hertz::<f32>::from_hz(0.1).unwrap());

        assert!(f1 < f2);
        assert!(f3 > f2);
        assert!(f1 == f1);
        assert!(f1 != f2);
    }

    #[test]
    fn test_frequency_f64() {
        let f1 = 10.hz();
        let f2 = 10.khz();
        let f3 = 10.mhz();
        let f4 = 10.dt();

        assert_eq!(f1, Hertz::<f64>::from_hz(10.).unwrap());
        assert_eq!(f2, Hertz::<f64>::from_hz(10000.).unwrap());
        assert_eq!(f3, Hertz::<f64>::from_hz(10000000.).unwrap());
        assert_eq!(f4, Hertz::<f64>::from_hz(0.1).unwrap());

        assert!(f1 < f2);
        assert!(f3 > f2);
        assert!(f1 == f1);
        assert!(f1 != f2);
    }

    #[test]
    #[should_panic]
    fn test_frequency_panic() {
        let _f1 = (-10.0).hz();
    }

    #[test]
    fn test_hertz_from_f32() {
        assert_eq!(
            Hertz::<f32>::from_dt(1.0).unwrap(),
            Hertz::<f32>::from_hz(1.0).unwrap()
        );
    }

    #[test]
    fn test_hertz_from_f64() {
        assert_eq!(
            Hertz::<f64>::from_dt(1.0).unwrap(),
            Hertz::<f64>::from_hz(1.0).unwrap()
        );
    }

    #[test]
    fn test_coefficients_normal_f32() {
        let f0 = 10.hz();
        let fs = 1.khz();

        let coeffs = Coefficients::<f32>::from_params(Type::LowPass, fs, f0, Q_BUTTERWORTH_F32);

        match coeffs {
            Ok(_) => {}
            Err(_) => {
                panic!("Coefficients creation failed!");
            }
        }
    }

    #[test]
    fn test_coefficients_normal_f64() {
        let f0 = 10.hz();
        let fs = 1.khz();

        let coeffs = Coefficients::<f64>::from_params(Type::LowPass, fs, f0, Q_BUTTERWORTH_F64);

        match coeffs {
            Ok(_) => {}
            Err(_) => {
                panic!("Coefficients creation failed!");
            }
        }
    }

    #[test]
    fn test_coefficients_fail_flipped_frequencies_f32() {
        let f0 = 10.hz();
        let fs = 1.khz();

        let coeffs = Coefficients::<f32>::from_params(Type::LowPass, f0, fs, Q_BUTTERWORTH_F32);

        match coeffs {
            Ok(_) => {
                panic!("Should not come here");
            }
            Err(e) => {
                assert_eq!(e, Errors::OutsideNyquist);
            }
        }
    }

    #[test]
    fn test_coefficients_fail_flipped_frequencies_f64() {
        let f0 = 10.hz();
        let fs = 1.khz();

        let coeffs = Coefficients::<f64>::from_params(Type::LowPass, f0, fs, Q_BUTTERWORTH_F64);

        match coeffs {
            Ok(_) => {
                panic!("Should not come here");
            }
            Err(e) => {
                assert_eq!(e, Errors::OutsideNyquist);
            }
        }
    }

    #[test]
    fn test_coefficients_fail_negative_q_f32() {
        let f0 = 10.hz();
        let fs = 1.khz();

        let coeffs = Coefficients::<f32>::from_params(Type::LowPass, fs, f0, -1.0);

        match coeffs {
            Ok(_) => {
                panic!("Should not come here");
            }
            Err(e) => {
                assert_eq!(e, Errors::NegativeQ);
            }
        }
    }

    #[test]
    fn test_coefficients_fail_negative_q_f64() {
        let f0 = 10.hz();
        let fs = 1.khz();

        let coeffs = Coefficients::<f64>::from_params(Type::LowPass, fs, f0, -1.0);

        match coeffs {
            Ok(_) => {
                panic!("Should not come here");
            }
            Err(e) => {
                assert_eq!(e, Errors::NegativeQ);
            }
        }
    }

    #[test]
    fn test_biquad_zeros_f32() {
        use std::vec::Vec;
        //use std::vec::Vec::*;

        let f0 = 10.hz();
        let fs = 1.khz();

        let coeffs =
            Coefficients::<f32>::from_params(Type::LowPass, fs, f0, Q_BUTTERWORTH_F32).unwrap();

        let mut biquad1 = DirectForm1::<f32>::new(coeffs);
        let mut biquad2 = DirectForm2Transposed::<f32>::new(coeffs);

        let input_vec = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut output_vec1 = Vec::new();
        let mut output_vec2 = Vec::new();

        for elem in input_vec {
            output_vec1.push(biquad1.run(elem));
            output_vec2.push(biquad2.run(elem));
        }
    }

    #[test]
    fn test_biquad_zeros_f64() {
        use std::vec::Vec;
        //use std::vec::Vec::*;

        let f0 = 10.hz();
        let fs = 1.khz();

        let coeffs =
            Coefficients::<f64>::from_params(Type::LowPass, fs, f0, Q_BUTTERWORTH_F64).unwrap();

        let mut biquad1 = DirectForm1::<f64>::new(coeffs);
        let mut biquad2 = DirectForm2Transposed::<f64>::new(coeffs);

        let input_vec = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut output_vec1 = Vec::new();
        let mut output_vec2 = Vec::new();

        for elem in input_vec {
            output_vec1.push(biquad1.run(elem));
            output_vec2.push(biquad2.run(elem));
        }
    }
}
