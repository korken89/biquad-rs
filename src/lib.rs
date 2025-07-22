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

use core::ops::{Add, Mul, Sub};

use num_traits::Zero;

pub use crate::coefficients::*;
pub use crate::frequency::*;

/// The required functions of a biquad implementation
pub trait Biquad<C, T = C> {
    /// A single iteration of a biquad, applying the filtering on the input

    fn run(&mut self, input: T) -> T;

    /// Updating of coefficients
    fn update_coefficients(&mut self, new_coefficients: Coefficients<C>);

    /// Updating coefficients and returning the old ones. This is useful to avoid deallocating on the audio thread, since
    /// the `Coefficients` can then be sent to another thread for deallocation.
    fn replace_coefficients(&mut self, new_coefficients: Coefficients<C>) -> Coefficients<C>;

    /// Set the internal state of the biquad to 0 without allocation.
    fn reset_state(&mut self);
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
pub struct DirectForm1<C, T = C> {
    y1: T,
    y2: T,
    x1: T,
    x2: T,
    coeffs: Coefficients<C>,
}

/// Internal states and coefficients of the Direct Form 2 Transposed form
#[derive(Copy, Clone, Debug)]
pub struct DirectForm2Transposed<C, T = C> {
    pub s1: T,
    pub s2: T,
    coeffs: Coefficients<C>,
}

impl<C, T> DirectForm1<C, T>
where
    T: Zero,
{
    /// Creates a Direct Form 1 biquad from a set of filter coefficients
    pub fn new(coefficients: Coefficients<C>) -> Self {
        DirectForm1 {
            y1: T::zero(),
            y2: T::zero(),
            x1: T::zero(),
            x2: T::zero(),
            coeffs: coefficients,
        }
    }
}

impl<C, T> Biquad<C, T> for DirectForm1<C, T>
where
    T: Copy + Add<T, Output = T> + Sub<T, Output = T> + Zero,
    C: Copy + Mul<T, Output = T>,
{
    fn run(&mut self, input: T) -> T {
        let out = self.coeffs.b0 * input + self.coeffs.b1 * self.x1 + self.coeffs.b2 * self.x2
            - self.coeffs.a1 * self.y1
            - self.coeffs.a2 * self.y2;

        self.x2 = self.x1;
        self.x1 = input;
        self.y2 = self.y1;
        self.y1 = out;

        out
    }

    fn update_coefficients(&mut self, new_coefficients: Coefficients<C>) {
        self.coeffs = new_coefficients;
    }

    fn replace_coefficients(&mut self, new_coefficients: Coefficients<C>) -> Coefficients<C> {
        core::mem::replace(&mut self.coeffs, new_coefficients)
    }

    fn reset_state(&mut self) {
        self.x1 = T::zero();
        self.x2 = T::zero();
        self.y1 = T::zero();
        self.y2 = T::zero();
    }
}

impl<C, T> DirectForm2Transposed<C, T> where T: Zero {
    /// Creates a Direct Form 2 Transposed biquad from a set of filter coefficients
    pub fn new(coefficients: Coefficients<C>) -> Self {
        DirectForm2Transposed {
            s1: T::zero(),
            s2: T::zero(),
            coeffs: coefficients,
        }
    }
}

impl<C, T> Biquad<C, T> for DirectForm2Transposed<C, T>
where
    T: Copy + Add<T, Output = T> + Sub<T, Output = T> + Zero,
    C: Copy + Mul<T, Output = T>,
{
    fn run(&mut self, input: T) -> T {
        let out = self.s1 + self.coeffs.b0 * input;
        self.s1 = self.s2 + self.coeffs.b1 * input - self.coeffs.a1 * out;
        self.s2 = self.coeffs.b2 * input - self.coeffs.a2 * out;

        out
    }

    fn update_coefficients(&mut self, new_coefficients: Coefficients<C>) {
        self.coeffs = new_coefficients;
    }

    fn replace_coefficients(&mut self, new_coefficients: Coefficients<C>) -> Coefficients<C> {
        core::mem::replace(&mut self.coeffs, new_coefficients)
    }

    fn reset_state(&mut self) {
        self.s1 = T::zero();
        self.s2 = T::zero();
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
        let _f1: Hertz<f32> = (-10.0f32).hz();
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
        let f0 = 100.hz();
        let fs = 1.khz();

        let coeffs = Coefficients::<f32>::from_params(Type::LowPass, fs, f0, Q_BUTTERWORTH_F32);
        println!("{:?}", coeffs.unwrap());
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

    #[test]
    fn test_biquad_lowpass_vs_scipy() {
        use rand::prelude::*;
        use std::process::Command;
        use std::string::String;
        use std::vec::Vec;
        let normalize_f0 = 0.2;
        let mut rng = rand::thread_rng();

        let coeffs: Coefficients<f32> = Coefficients::<f32>::from_normalized_params(
            Type::LowPass,
            normalize_f0,
            Q_BUTTERWORTH_F32,
        )
        .unwrap();

        let mut biquad = DirectForm1::<f32>::new(coeffs);
        let mut in_vec = Vec::<f32>::with_capacity(1000);
        let mut out_vec = Vec::<f32>::with_capacity(1000);
        for i in 0..in_vec.capacity() {
            let x = ((i as f32) * 0.1).sin() + rng.gen::<f32>();
            in_vec.push(x);
            out_vec.push(biquad.run(x));
        }
        let in_vec_str = format!("{:?}", in_vec);
        let out_vec_str = format!("{:?}", out_vec);

        let cmd_output = Command::new("python")
            .args([
                "test/comparaison_lowpass_to_scipy.py",
                format!("{normalize_f0}").as_str(),
                in_vec_str.as_str(),
                out_vec_str.as_str(),
            ])
            .output()
            .expect("You need python and scipy and numpy to run this one.");
        let py_res_stderr = String::from_utf8(cmd_output.stderr).unwrap();
        if !py_res_stderr.is_empty() {
            println!("{:?}", py_res_stderr);
        }
        assert!(py_res_stderr.is_empty());
        let py_res = String::from_utf8(cmd_output.stdout).unwrap();
        let py_res_parsed_vec: Vec<_> = py_res
            .split(",")
            .map_while(|x| {
                let a = x.trim().parse::<f32>();
                if let Ok(fl) = a {
                    return Some(fl);
                }
                return None;
            })
            .collect();
        println!("{}", py_res_parsed_vec.len());
        let sum_err: f32 = py_res_parsed_vec
            .iter()
            .zip(out_vec)
            .map(|(x, y)| (x - y).abs())
            .sum();
        let avg_error = sum_err / (in_vec.len() as f32);
        println!("Sum of absolute error = {sum_err}");

        println!("Average absolute error = {avg_error}");

        assert!(avg_error <= 1e-6);
    }
    #[test]
    fn test_biquad_notch_vs_scipy() {
        use core::f32::consts::PI;
        use rand::prelude::*;
        use std::process::Command;
        use std::string::String;
        use std::vec::Vec;
        let mut rng = rand::thread_rng();
        let f0 = 0.2;
        let w0 = PI * f0;
        let normalize_f01 = f0 / 2.;
        let normalize_f02 = f0 * 2.;
        let coeffs: Coefficients<f32> = Coefficients::<f32>::band_0db_from_cutting_frequencies(
            Type::BandPass,
            normalize_f01,
            normalize_f02,
        )
        .unwrap();
        let mut biquad = DirectForm1::<f32>::new(coeffs);
        let vec_capacity = 1000;
        let mut in_vec = Vec::<f32>::with_capacity(vec_capacity);
        let mut out_vec = Vec::<f32>::with_capacity(vec_capacity);
        for i in 0..vec_capacity {
            let x = ((i as f32) * w0).sin() + rng.gen::<f32>();
            in_vec.push(x);
            out_vec.push(biquad.run(x));
        }
        let in_vec_str = format!("{:?}", in_vec);
        let out_vec_str = format!("{:?}", out_vec);

        let cmd_output = Command::new("python")
            .args([
                "test/comparaison_notch_to_scipy.py",
                format!("{normalize_f01}").as_str(),
                format!("{normalize_f02}").as_str(),
                in_vec_str.as_str(),
                out_vec_str.as_str(),
            ])
            .output()
            .expect("You need python and scipy and numpy to run this one.");
        let py_res_stderr = String::from_utf8(cmd_output.stderr).unwrap();
        if !py_res_stderr.is_empty() {
            println!("{:?}", py_res_stderr);
        }
        let py_res = String::from_utf8(cmd_output.stdout).unwrap();
        let py_res_parsed_vec: Vec<_> = py_res
            .split(",")
            .map_while(|x| {
                let a = x.trim().parse::<f32>();
                if let Ok(fl) = a {
                    return Some(fl);
                }
                return None;
            })
            .collect();
        let sum_err: f32 = py_res_parsed_vec
            .iter()
            .zip(out_vec)
            .map(|(x, y)| (x - y).abs())
            .sum();
        let avg_error = sum_err / (in_vec.len() as f32);
        println!("Sum of absolute error = {sum_err}");
        println!("Average absolute error = {avg_error}");

        assert!(avg_error <= 1e-1);
        assert!(py_res_stderr.is_empty());
    }
}
