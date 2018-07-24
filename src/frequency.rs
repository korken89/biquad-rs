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
//! `Hertz::new(...)` can error if the frequency is negative.

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
        Hertz::new(self).unwrap()
    }

    fn khz(self) -> Hertz {
        Hertz::new(self * 1_000.0).unwrap()
    }

    fn mhz(self) -> Hertz {
        Hertz::new(self * 1_000_000.0).unwrap()
    }

    fn dt(self) -> Hertz {
        Hertz::new(1.0 / self).unwrap()
    }
}

impl ToHertz for u32 {
    fn hz(self) -> Hertz {
        Hertz::new(self as f32).unwrap()
    }

    fn khz(self) -> Hertz {
        Hertz::new((self * 1_000) as f32).unwrap()
    }

    fn mhz(self) -> Hertz {
        Hertz::new((self * 1_000_000) as f32).unwrap()
    }

    fn dt(self) -> Hertz {
        Hertz::new(1 as f32 / self as f32).unwrap()
    }
}

impl Hertz {
    pub fn new(hz: f32) -> Result<Hertz, &'static str> {
        if hz < 0.0 {
            return Err("Negative frequency");
        }

        Ok(Hertz(hz))
    }

    pub fn hz(self) -> f32 {
        self.0
    }
}
