//! # coefficients
//!
//! Module for generating filter coefficients for second order IIR biquads, where the coefficients
//! form the following Z-domain transfer function:
//! ```text
//!         b0 + b1 * z^-1 + b2 * z^-2
//! H(z) =  --------------------------
//!          1 + a1 * z^-1 + a2 * z^-2
//! ```
//!
//! The second orders filter are based on the
//! [Audio EQ Cookbook](http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt), while the first order
//! low pass filter is based on the following
//! [Wikipedia article](https://en.wikipedia.org/wiki/Low-pass_filter#Discrete-time_realization).
//!
//!
//! # Examples
//!
//! ```
//! fn main() {
//!     use biquad::*;
//!
//!     // Cutoff frequency
//!     let f0 = 10.hz();
//!
//!     // Sampling frequency
//!     let fs = 1.khz();
//!
//!     // Create coefficients
//!     let coeffs = Coefficients::<f32>::from_params(Type::LowPass, fs, f0, Q_BUTTERWORTH_F32);
//! }
//! ```
//!
//! # Errors
//!
//! `Coefficients::from_params(...)` can error if the cutoff frequency does not adhere to the
//! [Nyquist Frequency](https://en.wikipedia.org/wiki/Nyquist_frequency), or if the Q value is
//! negative.

use crate::{frequency::Hertz, Errors};

// For some reason this is not detected properly
#[allow(unused_imports)]
use libm::{F32Ext, F64Ext};

/// Common Q value of the Butterworth low-pass filter
pub const Q_BUTTERWORTH_F32: f32 = core::f32::consts::FRAC_1_SQRT_2;
pub const Q_BUTTERWORTH_F64: f64 = core::f64::consts::FRAC_1_SQRT_2;

/// The supported types of biquad coefficients. Note that single pole low pass filters are faster to
/// retune, as all other filter types require evaluations of sin/cos functions
#[derive(Clone, Copy, Debug)]
pub enum Type {
    SinglePoleLowPass,
    LowPass,
    HighPass,
    BandPass,
    Notch,
}

/// Holder of the biquad coefficients, utilizes normalized form
#[derive(Clone, Copy, Debug)]
pub struct Coefficients<T> {
    // Denominator coefficients
    pub a1: T,
    pub a2: T,

    // Nominator coefficients
    pub b0: T,
    pub b1: T,
    pub b2: T,
}

impl Coefficients<f32> {
    /// Creates coefficients based on the biquad filter type, sampling and cutoff frequency, and Q
    /// value. Note that the cutoff frequency must be smaller than half the sampling frequency and
    /// that Q may not be negative, this will result in an `Err()`.
    pub fn from_params(
        filter: Type,
        fs: Hertz<f32>,
        f0: Hertz<f32>,
        q_value: f32,
    ) -> Result<Coefficients<f32>, Errors> {
        if 2.0 * f0.hz() > fs.hz() {
            return Err(Errors::OutsideNyquist);
        }

        if q_value < 0.0 {
            return Err(Errors::NegativeQ);
        }

        let omega = 2.0 * core::f32::consts::PI * f0.hz() / fs.hz();

        match filter {
            Type::SinglePoleLowPass => {
                let alpha = omega / (omega + 1.0);

                Ok(Coefficients {
                    a1: alpha - 1.0,
                    a2: 0.0,
                    b0: alpha,
                    b1: 0.0,
                    b2: 0.0,
                })
            }
            Type::LowPass => {
                // The code for omega_s/c and alpha is currently duplicated due to the single pole
                // low pass filter not needing it and when creating coefficients are commonly
                // assumed to be of low computational complexity.
                let omega_s = omega.sin();
                let omega_c = omega.cos();
                let alpha = omega_s / (2.0 * q_value);

                let b0 = (1.0 - omega_c) * 0.5;
                let b1 = 1.0 - omega_c;
                let b2 = (1.0 - omega_c) * 0.5;
                let a0 = 1.0 + alpha;
                let a1 = -2.0 * omega_c;
                let a2 = 1.0 - alpha;

                Ok(Coefficients {
                    a1: a1 / a0,
                    a2: a2 / a0,
                    b0: b0 / a0,
                    b1: b1 / a0,
                    b2: b2 / a0,
                })
            }
            Type::HighPass => {
                let omega_s = omega.sin();
                let omega_c = omega.cos();
                let alpha = omega_s / (2.0 * q_value);

                let b0 = (1.0 + omega_c) * 0.5;
                let b1 = -(1.0 + omega_c);
                let b2 = (1.0 + omega_c) * 0.5;
                let a0 = 1.0 + alpha;
                let a1 = -2.0 * omega_c;
                let a2 = 1.0 - alpha;

                Ok(Coefficients {
                    a1: a1 / a0,
                    a2: a2 / a0,
                    b0: b0 / a0,
                    b1: b1 / a0,
                    b2: b2 / a0,
                })
            }
            Type::BandPass => {
                let omega_s = omega.sin();
                let omega_c = omega.cos();
                let alpha = omega_s / (2.0 * q_value);

                let b0 = omega_s / 2.0;
                let b1 = 0.;
                let b2 = -(omega_s / 2.0);
                let a0 = 1.0 + alpha;
                let a1 = -2.0 * omega_c;
                let a2 = 1.0 - alpha;

                let div = 1.0 / a0;

                Ok(Coefficients {
                    a1: a1 * div,
                    a2: a2 * div,
                    b0: b0 * div,
                    b1: b1 * div,
                    b2: b2 * div,
                })
            }
            Type::Notch => {
                let omega_s = omega.sin();
                let omega_c = omega.cos();
                let alpha = omega_s / (2.0 * q_value);

                let b0 = 1.0;
                let b1 = -2.0 * omega_c;
                let b2 = 1.0;
                let a0 = 1.0 + alpha;
                let a1 = -2.0 * omega_c;
                let a2 = 1.0 - alpha;

                Ok(Coefficients {
                    a1: a1 / a0,
                    a2: a2 / a0,
                    b0: b0 / a0,
                    b1: b1 / a0,
                    b2: b2 / a0,
                })
            }
        }
    }
}

impl Coefficients<f64> {
    /// Creates coefficients based on the biquad filter type, sampling and cutoff frequency, and Q
    /// value. Note that the cutoff frequency must be smaller than half the sampling frequency and
    /// that Q may not be negative, this will result in an `Err()`.
    pub fn from_params(
        filter: Type,
        fs: Hertz<f64>,
        f0: Hertz<f64>,
        q_value: f64,
    ) -> Result<Coefficients<f64>, Errors> {
        if 2.0 * f0.hz() > fs.hz() {
            return Err(Errors::OutsideNyquist);
        }

        if q_value < 0.0 {
            return Err(Errors::NegativeQ);
        }

        let omega = 2.0 * core::f64::consts::PI * f0.hz() / fs.hz();

        match filter {
            Type::SinglePoleLowPass => {
                let alpha = omega / (omega + 1.0);

                Ok(Coefficients {
                    a1: alpha - 1.0,
                    a2: 0.0,
                    b0: alpha,
                    b1: 0.0,
                    b2: 0.0,
                })
            }
            Type::LowPass => {
                // The code for omega_s/c and alpha is currently duplicated due to the single pole
                // low pass filter not needing it and when creating coefficients are commonly
                // assumed to be of low computational complexity.
                let omega_s = omega.sin();
                let omega_c = omega.cos();
                let alpha = omega_s / (2.0 * q_value);

                let b0 = (1.0 - omega_c) * 0.5;
                let b1 = 1.0 - omega_c;
                let b2 = (1.0 - omega_c) * 0.5;
                let a0 = 1.0 + alpha;
                let a1 = -2.0 * omega_c;
                let a2 = 1.0 - alpha;

                let div = 1.0 / a0;

                Ok(Coefficients {
                    a1: a1 * div,
                    a2: a2 * div,
                    b0: b0 * div,
                    b1: b1 * div,
                    b2: b2 * div,
                })
            }
            Type::HighPass => {
                let omega_s = omega.sin();
                let omega_c = omega.cos();
                let alpha = omega_s / (2.0 * q_value);

                let b0 = (1.0 + omega_c) * 0.5;
                let b1 = -(1.0 + omega_c);
                let b2 = (1.0 + omega_c) * 0.5;
                let a0 = 1.0 + alpha;
                let a1 = -2.0 * omega_c;
                let a2 = 1.0 - alpha;

                let div = 1.0 / a0;

                Ok(Coefficients {
                    a1: a1 * div,
                    a2: a2 * div,
                    b0: b0 * div,
                    b1: b1 * div,
                    b2: b2 * div,
                })
            }
            Type::Notch => {
                let omega_s = omega.sin();
                let omega_c = omega.cos();
                let alpha = omega_s / (2.0 * q_value);

                let b0 = 1.0;
                let b1 = -2.0 * omega_c;
                let b2 = 1.0;
                let a0 = 1.0 + alpha;
                let a1 = -2.0 * omega_c;
                let a2 = 1.0 - alpha;

                let div = 1.0 / a0;

                Ok(Coefficients {
                    a1: a1 * div,
                    a2: a2 * div,
                    b0: b0 * div,
                    b1: b1 * div,
                    b2: b2 * div,
                })
            }
            Type::BandPass => {
                let omega_s = omega.sin();
                let omega_c = omega.cos();
                let alpha = omega_s / (2.0 * q_value);

                let b0 = omega_s / 2.0;
                let b1 = 0.;
                let b2 = -(omega_s / 2.0);
                let a0 = 1.0 + alpha;
                let a1 = -2.0 * omega_c;
                let a2 = 1.0 - alpha;

                let div = 1.0 / a0;

                Ok(Coefficients {
                    a1: a1 * div,
                    a2: a2 * div,
                    b0: b0 * div,
                    b1: b1 * div,
                    b2: b2 * div,
                })
            }
        }
    }
}
