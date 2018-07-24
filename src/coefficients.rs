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
//! # Usage
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
//!     let coeffs = Coefficients::new(Type::LowPass, fs, f0, Q_BUTTERWORTH);
//! }
//! ```

use core::f32::consts::FRAC_1_SQRT_2;
use core::f32::consts::FRAC_PI_2;
use core::f32::consts::PI;
use frequency::Hertz;

/// Common Q value of the Butterworth low-pass filter
pub const Q_BUTTERWORTH: f32 = FRAC_1_SQRT_2;

/// The supported types of biquad coefficients. Note that single pole low pass filters are faster to
/// retune, as all other filter types require evaluations of sin/cos functions
#[derive(Clone, Copy, Debug)]
pub enum Type {
    SinglePoleLowPass,
    LowPass,
    HighPass,
    Notch,
}

/// Holder of the biquad coefficients, utilizes normalized form
#[derive(Clone, Copy, Debug)]
pub struct Coefficients {
    // Denominator coefficients
    pub a1: f32,
    pub a2: f32,

    // Nominator coefficients
    pub b0: f32,
    pub b1: f32,
    pub b2: f32,
}

/// An [accurate and simple sin function](http://mooooo.ooo/chebyshev-sine-approximation/)
/// over the input in [-pi, pi] which is accurate to 4.58 ULPs. The sin/cos calculations used in the
/// implementation of coefficients are bounded to [0, pi] which allows for the use of simplified
/// sin/cos implementations.
fn sin(x: f32) -> f32 {
    let coeffs = [
        -0.10132118f32,         // x
        0.0066208798f32,        // x^3
        -0.00017350505f32,      // x^5
        0.0000025222919f32,     // x^7
        -0.000000023317787f32,  // x^9
        0.00000000013291342f32, // x^11
    ];
    let pi_major = 3.1415927f32;
    let pi_minor = -0.00000008742278f32;

    // Horner's rule
    let x2 = x * x;
    let p11 = coeffs[5];
    let p9 = p11 * x2 + coeffs[4];
    let p7 = p9 * x2 + coeffs[3];
    let p5 = p7 * x2 + coeffs[2];
    let p3 = p5 * x2 + coeffs[1];
    let p1 = p3 * x2 + coeffs[0];

    // Remove scaling function
    (x - pi_major - pi_minor) * (x + pi_major + pi_minor) * p1 * x
}

/// Used to convert the argument of `cos` into a `sin`, with the arguments being in [0, pi], the
/// transformation here will still be valid for `sin`.
fn cos(x: f32) -> f32 {
    sin(FRAC_PI_2 - x)
}

impl Coefficients {
    /// Creates coefficients based on the biquad filter type, sampling and cutoff frequency, and Q
    /// value. Note that the cutoff frequency must be smaller than half the sampling frequency and
    /// that Q may not be negative, this will result in an `Err(&str)`.
    pub fn new(
        filter: Type,
        fs: Hertz,
        f0: Hertz,
        q_value: f32,
    ) -> Result<Coefficients, &'static str> {
        if f0.hz() > fs.hz() / 2.0 {
            return Err("Center frequency outside half of the sampling frequency");
        }

        if q_value < 0.0 {
            return Err("Negative Q value");
        }

        let omega = 2.0 * PI * f0.hz() / fs.hz();

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
                let omega_s = sin(omega);
                let omega_c = cos(omega);
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
                let omega_s = sin(omega);
                let omega_c = cos(omega);
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
            Type::Notch => {
                let omega_s = sin(omega);
                let omega_c = cos(omega);
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
