//! Implmentation of DFT functionals.

const std = @import("std");

/// Enumeration of available DFT functionals.
pub fn DFTFunctional(comptime T: type) type {
    return union(enum) {
        lda: LDA(T), LDA: LDA(T)
    };
}

/// Return type of the DFT functional.
pub fn DFTResult(comptime T: type) type {
    return struct {eps_xc: T, v_xc: T};
}

/// Compute the exchange-correlation energy density and potential for a given electron density `rho` using the specified DFT functional.
pub fn computeExchangeCorrelation(comptime T: type, functional: DFTFunctional(T), rho: T) DFTResult(T) {
    switch (functional) {
        inline else => |dft| return dft.evaluate(rho)
    }
}

/// The Local Density Approximation (LDA) functional for exchange energy.
pub fn LDA(comptime T: type) type {
    return struct {

        /// Evaluate function.
        pub fn evaluate(_: @This(), rho: T) DFTResult(T) {
            if (rho <= 1e-12) return .{.eps_xc = 0, .v_xc = 0};

            const cx_energy = -0.75 * std.math.pow(T, 3.0 / std.math.pi, 1.0 / 3.0);
            const cx_potential = -std.math.pow(T, 3.0 / std.math.pi, 1.0 / 3.0);

            const rho_third = std.math.pow(T, rho, 1.0 / 3.0);

            return .{
                .eps_xc = cx_energy * rho_third,
                .v_xc = cx_potential * rho_third,
            };
        }
    };
}
