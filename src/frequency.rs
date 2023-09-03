//! # frequency
//!
//! A helper module for creating type-safe frequencies, while also allowing to create frequencies
//! from float and integer literals.
//!
//! # Examples
//!
//! ```ignore
//! fn main() {
//!     use biquad::frequency::*;
//!
//!     // Integer literals
//!     let one_hz = 1.hz();
//!     let one_khz = 1.khz();
//!     let one_mhz = 1.mhz();
//!
//!     // Float literals
//!     let one_hz = 1.0.dt();
//!     let ten_hz = 0.1.dt();
//! }
//! ```
//!
//! # Errors
//!
//! `Hertz::from_hz(...)` will error if the frequency is negative.
//!
//! # Panics
//!
//! `x.hz()`, `x.khz()`, `x.mhz()`, `x.dt()` will panic for `f32` if they are negative.
//!

use num_traits::Float;

use crate::Errors;

/// Base type for frequency, everything is based on Hertz
#[derive(PartialOrd, PartialEq, Debug, Copy, Clone)]
pub struct Hertz<T: Float>(T);

/// Used to implement conversions to the Hertz struct
pub trait ToHertz<T: Float> {
    /// From hertz
    fn hz(self) -> Hertz<T>;

    /// From kilohertz
    fn khz(self) -> Hertz<T>;

    /// From megahertz
    fn mhz(self) -> Hertz<T>;

    /// From delta time (in seconds)
    fn dt(self) -> Hertz<T>;
}

// -----------------------------------------------
// generic implementation
// -----------------------------------------------
use num_traits::NumCast;
impl<T: Float, X: NumCast> ToHertz<T> for X {
    fn hz(self) -> Hertz<T> {
        Hertz::<T>::from_hz(T::from(self).unwrap()).unwrap()
    }

    fn khz(self) -> Hertz<T> {
        Hertz::<T>::from_hz(T::from(self).unwrap() * T::from(1_000).unwrap()).unwrap()
    }

    fn mhz(self) -> Hertz<T> {
        Hertz::<T>::from_hz(T::from(self).unwrap() * T::from(1_000_000).unwrap()).unwrap()
    }

    fn dt(self) -> Hertz<T> {
        Hertz::<T>::from_hz(T::one() / T::from(self).unwrap()).unwrap()
    }
}

impl<T: Float> Hertz<T> {
    pub fn from_hz(hz: T) -> Result<Self, Errors> {
        if hz > T::zero() {
            Ok(Hertz(hz))
        } else {
            Err(Errors::NegativeFrequency)
        }
    }

    pub fn from_dt(dt: T) -> Result<Self, Errors> {
        if dt > T::zero() {
            Ok(Hertz(T::one() / dt))
        } else {
            Err(Errors::NegativeFrequency)
        }
    }

    pub fn hz(self) -> T {
        self.0
    }
}
