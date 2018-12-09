//! # frequency
//!
//! A helper module for creating type-safe frequencies, while also allowing to create frequencies
//! from float and integer literals.
//!
//! # Examples
//!
//! ```
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

use crate::Errors;

/// Base type for frequency, everything is based on Hertz
#[derive(PartialOrd, PartialEq, Debug, Copy, Clone)]
pub struct Hertz(f32);

/// Used to implement conversions to the Hertz struct
pub trait ToHertz {
    /// From hertz
    fn hz(self) -> Hertz;

    /// From kilohertz
    fn khz(self) -> Hertz;

    /// From megahertz
    fn mhz(self) -> Hertz;

    /// From delta time (in seconds)
    fn dt(self) -> Hertz;
}

impl ToHertz for f32 {
    fn hz(self) -> Hertz {
        Hertz::from_hz(self).unwrap()
    }

    fn khz(self) -> Hertz {
        Hertz::from_hz(self * 1_000.0).unwrap()
    }

    fn mhz(self) -> Hertz {
        Hertz::from_hz(self * 1_000_000.0).unwrap()
    }

    fn dt(self) -> Hertz {
        Hertz::from_hz(1.0 / self).unwrap()
    }
}

impl ToHertz for u32 {
    fn hz(self) -> Hertz {
        Hertz::from_hz(self as f32).unwrap()
    }

    fn khz(self) -> Hertz {
        Hertz::from_hz((self * 1_000) as f32).unwrap()
    }

    fn mhz(self) -> Hertz {
        Hertz::from_hz((self * 1_000_000) as f32).unwrap()
    }

    fn dt(self) -> Hertz {
        Hertz::from_hz(1 as f32 / self as f32).unwrap()
    }
}

impl Hertz {
    pub fn from_hz(hz: f32) -> Result<Self, Errors> {
        if hz > 0.0 {
            Ok(Hertz(hz))
        } else {
            Err(Errors::NegativeFrequency)
        }
    }

    pub fn from_dt(dt: f32) -> Result<Self, Errors> {
        if dt > 0.0 {
            Ok(Hertz(1.0 / dt))
        } else {
            Err(Errors::NegativeFrequency)
        }
    }

    pub fn hz(self) -> f32 {
        self.0
    }
}
