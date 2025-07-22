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
//! [Audio EQ Cookbook](https://webaudio.github.io/Audio-EQ-Cookbook/audio-eq-cookbook.html), while the first order
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

use num_traits::Float;
/// Common Q value of the Butterworth low-pass filter
pub const Q_BUTTERWORTH_F32: f32 = core::f32::consts::FRAC_1_SQRT_2;
pub const Q_BUTTERWORTH_F64: f64 = core::f64::consts::FRAC_1_SQRT_2;

/// The supported types of biquad coefficients. Note that single pole low pass filters are faster to
/// retune, as all other filter types require evaluations of sin/cos functions
/// The `LowShelf`, `HighShelf`, and `PeakingEQ` all have a gain value for its
/// field, and represents the gain, in decibels, that the filter provides.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Type<DBGain> {
    SinglePoleLowPassApprox,
    SinglePoleLowPass,
    LowPass,
    HighPass,
    BandPass,
    Notch,
    AllPass,
    LowShelf(DBGain),
    HighShelf(DBGain),
    PeakingEQ(DBGain),
}

/// Holder of the biquad coefficients, utilizes normalized form
#[derive(Clone, Copy, Debug)]
pub struct Coefficients<C> {
    // Denominator coefficients
    pub a1: C,
    pub a2: C,

    // Nominator coefficients
    pub b0: C,
    pub b1: C,
    pub b2: C,
}

impl<C> Coefficients<C> where C: Float {
    /// Creates coefficients based on the biquad filter type, normalized cutoff frequency, and Q
    /// value. Note that the cutoff frequency must be smaller than 1 and that Q may not be negative,
    ///  this will result in an `Err()`.
    /// * `filter_type` - Type of filter desired
    /// * `normalized_f0` - Cut frequency devided by two times sampling frequency (0<F_cut<1 respects the shanon condition)
    /// * `q_value` - Text about bar.
    pub fn from_normalized_params(
        filter_type: Type<C>,
        normalized_f0: C,
        q_value: C,
    ) -> Result<Coefficients<C>, Errors> {
        #[allow(non_snake_case)]
        let TWO: C = C::from(2).unwrap();
        #[allow(non_snake_case)]
        let PI: C = C::from(core::f64::consts::PI).unwrap();
        #[allow(non_snake_case)]
        let FORTY: C = C::from(40).unwrap();
        if normalized_f0 >= C::one() || normalized_f0 < C::zero() {
            return Err(Errors::OutsideNyquist);
        }
        if q_value <= C::zero() {
            return Err(Errors::NegativeQ);
        }
        let omega = PI * normalized_f0;

        // The code for omega_s/c and alpha is hidden behind a condition due to the single pole
        // low pass filter not needing it and when creating coefficients are commonly
        // assumed to be of low computational complexity.
        let (omega_s, omega_c, alpha) = match filter_type {
            Type::SinglePoleLowPassApprox | Type::SinglePoleLowPass => {
                (C::nan(), C::nan(), C::nan())
            }
            _ => {
                let omega_s = omega.sin();
                (omega_s, omega.cos(), omega_s / (TWO * q_value))
            }
        };
        match filter_type {
            Type::SinglePoleLowPassApprox => {
                let alpha = omega / (omega + C::one());

                Ok(Coefficients {
                    a1: alpha - C::one(),
                    a2: C::zero(),
                    b0: alpha,
                    b1: C::zero(),
                    b2: C::zero(),
                })
            }
            Type::SinglePoleLowPass => {
                let omega_t = (omega / TWO).tan();
                let a0 = C::one() + omega_t;

                Ok(Coefficients {
                    a1: (omega_t - C::one()) / a0,
                    a2: C::zero(),
                    b0: omega_t / a0,
                    b1: omega_t / a0,
                    b2: C::zero(),
                })
            }
            Type::LowPass => {
                let b0 = (C::one() - omega_c) / TWO;
                let b1 = C::one() - omega_c;
                let b2 = (C::one() - omega_c) / TWO;
                let a0 = C::one() + alpha;
                let a1 = -TWO * omega_c;
                let a2 = C::one() - alpha;

                let div = C::one() / a0;

                Ok(Coefficients {
                    a1: a1 * div,
                    a2: a2 * div,
                    b0: b0 * div,
                    b1: b1 * div,
                    b2: b2 * div,
                })
            }
            Type::HighPass => {
                let b0 = (C::one() + omega_c) / TWO;
                let b1 = -(C::one() + omega_c);
                let b2 = (C::one() + omega_c) / TWO;
                let a0 = C::one() + alpha;
                let a1 = -TWO * omega_c;
                let a2 = C::one() - alpha;

                let div = C::one() / a0;

                Ok(Coefficients {
                    a1: a1 * div,
                    a2: a2 * div,
                    b0: b0 * div,
                    b1: b1 * div,
                    b2: b2 * div,
                })
            }
            Type::BandPass => {
                let b0 = omega_s / TWO;
                let b1 = C::zero();
                let b2 = -(omega_s / TWO);
                let a0 = C::one() + alpha;
                let a1 = -TWO * omega_c;
                let a2 = C::one() - alpha;

                let div = C::one() / a0;

                Ok(Coefficients {
                    a1: a1 * div,
                    a2: a2 * div,
                    b0: b0 * div,
                    b1: b1 * div,
                    b2: b2 * div,
                })
            }
            Type::Notch => {
                let b0 = C::one();
                let b1 = -TWO * omega_c;
                let b2 = C::one();
                let a0 = C::one() + alpha;
                let a1 = -TWO * omega_c;
                let a2 = C::one() - alpha;

                let div = C::one() / a0;

                Ok(Coefficients {
                    a1: a1 * div,
                    a2: a2 * div,
                    b0: b0 * div,
                    b1: b1 * div,
                    b2: b2 * div,
                })
            }
            Type::AllPass => {
                let b0 = C::one() - alpha;
                let b1 = -TWO * omega_c;
                let b2 = C::one() + alpha;
                let a0 = C::one() + alpha;
                let a1 = -TWO * omega_c;
                let a2 = C::one() - alpha;

                let div = C::one() / a0;

                Ok(Coefficients {
                    a1: a1 * div,
                    a2: a2 * div,
                    b0: b0 * div,
                    b1: b1 * div,
                    b2: b2 * div,
                })
            }
            Type::LowShelf(db_gain) => {
                let a = C::from(C::from(10).unwrap().powf(db_gain / FORTY)).unwrap();
                let b0 = a * ((a + C::one()) - (a - C::one()) * omega_c + TWO * alpha * a.sqrt());
                let b1 = TWO * a * ((a - C::one()) - (a + C::one()) * omega_c);
                let b2 = a * ((a + C::one()) - (a - C::one()) * omega_c - TWO * alpha * a.sqrt());
                let a0 = (a + C::one()) + (a - C::one()) * omega_c + TWO * alpha * a.sqrt();
                let a1 = -TWO * ((a - C::one()) + (a + C::one()) * omega_c);
                let a2 = (a + C::one()) + (a - C::one()) * omega_c - TWO * alpha * a.sqrt();

                let div = C::one() / a0;

                Ok(Coefficients {
                    a1: a1 * div,
                    a2: a2 * div,
                    b0: b0 * div,
                    b1: b1 * div,
                    b2: b2 * div,
                })
            }
            Type::HighShelf(db_gain) => {
                let a = C::from(C::from(10).unwrap().powf(db_gain / FORTY)).unwrap();

                let b0 = a * ((a + C::one()) + (a - C::one()) * omega_c + TWO * alpha * a.sqrt());
                let b1 = -TWO * a * ((a - C::one()) + (a + C::one()) * omega_c);
                let b2 = a * ((a + C::one()) + (a - C::one()) * omega_c - TWO * alpha * a.sqrt());
                let a0 = (a + C::one()) - (a - C::one()) * omega_c + TWO * alpha * a.sqrt();
                let a1 = TWO * ((a - C::one()) - (a + C::one()) * omega_c);
                let a2 = (a + C::one()) - (a - C::one()) * omega_c - TWO * alpha * a.sqrt();

                let div = C::one() / a0;

                Ok(Coefficients {
                    a1: a1 * div,
                    a2: a2 * div,
                    b0: b0 * div,
                    b1: b1 * div,
                    b2: b2 * div,
                })
            }
            Type::PeakingEQ(db_gain) => {
                let a = C::from(C::from(10).unwrap().powf(db_gain / FORTY)).unwrap();

                let b0 = C::one() + alpha * a;
                let b1 = -TWO * omega_c;
                let b2 = C::one() - alpha * a;
                let a0 = C::one() + alpha / a;
                let a1 = -TWO * omega_c;
                let a2 = C::one() - alpha / a;

                let div = C::one() / a0;

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
    /// Creates coefficients based on the biquad filter type, sampling and cutoff frequency, and Q
    /// value. Note that the cutoff frequency must be smaller than half the sampling frequency and
    /// that Q may not be negative, this will result in an `Err()`.
    pub fn from_params(
        filter: Type<C>,
        fs: Hertz<C>,
        f0: Hertz<C>,
        q_value: C,
    ) -> Result<Coefficients<C>, Errors> {
        let normalized_f0 = f0.hz() / (fs.hz() / C::from(2).unwrap());
        Self::from_normalized_params(filter, normalized_f0, q_value)
    }

    //
    pub fn band_0db_from_cutting_frequencies(
        filter_type: Type<C>,
        normalized_f01: C,
        normalized_f02: C,
    ) -> Result<Coefficients<C>, Errors> {
        if normalized_f01 < C::zero()
            || normalized_f02 < C::zero()
            || normalized_f01 > C::one()
            || normalized_f02 > C::one()
            || normalized_f01 > normalized_f02
        {
            return Err(Errors::OutsideNyquist);
        }
        #[allow(non_snake_case)]
        let TWO: C = C::from(2).unwrap();
        assert!(filter_type == Type::<_>::BandPass || filter_type == Type::<_>::Notch);
        let frac_freq: C = normalized_f02 / normalized_f01;
        let bw_in_octaves = frac_freq.log2();
        let normalized_f0 = TWO.powf(normalized_f01.log2() + bw_in_octaves / TWO);
        // #[allow(non_snake_case)]
        // let PI: T = T::from(core::f64::consts::PI).unwrap();
        // let omega_0 = PI * normalized_f0;
        //
        // let q_value =
        //     T::one() / (TWO * (TWO.ln() * bw_in_octaves * (omega_0 / omega_0.sin()) / TWO).sinh());
        Self::from_normalized_params(filter_type, normalized_f0, C::one())
    }
}
