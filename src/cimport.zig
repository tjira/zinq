//! Imports external C interfaces for numerical computations, Fourier transforms, and quantum chemical simulations.

pub const cblas = @import("cblas");
pub const exprtk = @import("exprtk");
pub const fftw = @import("fftw");
pub const lapacke = @import("lapacke");
pub const libint = @import("libint");
pub const libxc = @import("libxc");
