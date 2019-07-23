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

use crate::Errors;

/// Base type for frequency, everything is based on Hertz
#[derive(PartialOrd, PartialEq, Debug, Copy, Clone)]
pub struct Hertz<T>(T);

/// Used to implement conversions to the Hertz struct
pub trait ToHertz<T> {
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
// f32 implementation
// -----------------------------------------------

impl ToHertz<f32> for f32 {
    fn hz(self) -> Hertz<f32> {
        Hertz::<f32>::from_hz(self).unwrap()
    }

    fn khz(self) -> Hertz<f32> {
        Hertz::<f32>::from_hz(self * 1_000.0).unwrap()
    }

    fn mhz(self) -> Hertz<f32> {
        Hertz::<f32>::from_hz(self * 1_000_000.0).unwrap()
    }

    fn dt(self) -> Hertz<f32> {
        Hertz::<f32>::from_hz(1.0 / self).unwrap()
    }
}

impl ToHertz<f32> for u32 {
    fn hz(self) -> Hertz<f32> {
        Hertz::<f32>::from_hz(self as f32).unwrap()
    }

    fn khz(self) -> Hertz<f32> {
        Hertz::<f32>::from_hz((self * 1_000) as f32).unwrap()
    }

    fn mhz(self) -> Hertz<f32> {
        Hertz::<f32>::from_hz((self * 1_000_000) as f32).unwrap()
    }

    fn dt(self) -> Hertz<f32> {
        Hertz::<f32>::from_hz(1 as f32 / self as f32).unwrap()
    }
}

impl ToHertz<f32> for i32 {
    fn hz(self) -> Hertz<f32> {
        Hertz::<f32>::from_hz(self as f32).unwrap()
    }

    fn khz(self) -> Hertz<f32> {
        Hertz::<f32>::from_hz((self * 1_000) as f32).unwrap()
    }

    fn mhz(self) -> Hertz<f32> {
        Hertz::<f32>::from_hz((self * 1_000_000) as f32).unwrap()
    }

    fn dt(self) -> Hertz<f32> {
        Hertz::<f32>::from_hz(1 as f32 / self as f32).unwrap()
    }
}

impl Hertz<f32> {
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

// -----------------------------------------------
// f64 implementation
// -----------------------------------------------

impl ToHertz<f64> for f64 {
    fn hz(self) -> Hertz<f64> {
        Hertz::<f64>::from_hz(self).unwrap()
    }

    fn khz(self) -> Hertz<f64> {
        Hertz::<f64>::from_hz(self * 1_000.0).unwrap()
    }

    fn mhz(self) -> Hertz<f64> {
        Hertz::<f64>::from_hz(self * 1_000_000.0).unwrap()
    }

    fn dt(self) -> Hertz<f64> {
        Hertz::<f64>::from_hz(1.0 / self).unwrap()
    }
}

impl ToHertz<f64> for f32 {
    fn hz(self) -> Hertz<f64> {
        Hertz::<f64>::from_hz(self as f64).unwrap()
    }

    fn khz(self) -> Hertz<f64> {
        Hertz::<f64>::from_hz(self as f64 * 1_000.0).unwrap()
    }

    fn mhz(self) -> Hertz<f64> {
        Hertz::<f64>::from_hz(self as f64 * 1_000_000.0).unwrap()
    }

    fn dt(self) -> Hertz<f64> {
        Hertz::<f64>::from_hz(1.0 / self as f64).unwrap()
    }
}

impl ToHertz<f64> for u32 {
    fn hz(self) -> Hertz<f64> {
        Hertz::<f64>::from_hz(self as f64).unwrap()
    }

    fn khz(self) -> Hertz<f64> {
        Hertz::<f64>::from_hz((self * 1_000) as f64).unwrap()
    }

    fn mhz(self) -> Hertz<f64> {
        Hertz::<f64>::from_hz((self * 1_000_000) as f64).unwrap()
    }

    fn dt(self) -> Hertz<f64> {
        Hertz::<f64>::from_hz(1 as f64 / self as f64).unwrap()
    }
}

impl ToHertz<f64> for i32 {
    fn hz(self) -> Hertz<f64> {
        Hertz::<f64>::from_hz(self as f64).unwrap()
    }

    fn khz(self) -> Hertz<f64> {
        Hertz::<f64>::from_hz((self * 1_000) as f64).unwrap()
    }

    fn mhz(self) -> Hertz<f64> {
        Hertz::<f64>::from_hz((self * 1_000_000) as f64).unwrap()
    }

    fn dt(self) -> Hertz<f64> {
        Hertz::<f64>::from_hz(1 as f64 / self as f64).unwrap()
    }
}

impl Hertz<f64> {
    pub fn from_hz(hz: f64) -> Result<Self, Errors> {
        if hz > 0.0 {
            Ok(Hertz(hz))
        } else {
            Err(Errors::NegativeFrequency)
        }
    }

    pub fn from_dt(dt: f64) -> Result<Self, Errors> {
        if dt > 0.0 {
            Ok(Hertz(1.0 / dt))
        } else {
            Err(Errors::NegativeFrequency)
        }
    }

    pub fn hz(self) -> f64 {
        self.0
    }
}
