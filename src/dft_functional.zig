//! Implmentation of DFT functionals.

const std = @import("std");

/// Enumeration of available correlation functionals.
pub fn CorrelationFunctional(comptime T: type) type {
    return union(enum) {
        chachiyo: ChachiyoCorrelation(T), Chachiyo: ChachiyoCorrelation(T)
    };
}

/// Enumeration of available exchange functionals.
pub fn ExchangeFunctional(comptime T: type) type {
    return union(enum) {
        slater: SlaterExchange(T), Slater: SlaterExchange(T)
    };
}

/// Compute the exchange-correlation energy density and potential for a given electron density `rho` using the specified DFT functional.
pub fn computeExchangeCorrelation(comptime T: type, exchange: ExchangeFunctional(T), correlation: CorrelationFunctional(T), rho: T) struct {eps_xc: T, v_xc: T} {
    const res_exchange = switch (exchange) {
        inline else => |x| x.evaluate(rho)
    };

    const res_correlation = switch (correlation) {
        inline else => |c| c.evaluate(rho)
    };

    return .{
        .eps_xc = res_exchange.eps_x + res_correlation.eps_c,
        .v_xc = res_exchange.v_x + res_correlation.v_c
    };
}

/// The Local Density Approximation (LDA) functional for exchange-correlation energy.
pub fn SlaterExchange(comptime T: type) type {
    return struct {

        /// Evaluate function.
        pub fn evaluate(_: @This(), rho: T) struct {eps_x: T, v_x: T} {
            if (rho <= 1e-12) return .{.eps_x = 0, .v_x = 0};

            const cx_energy = -0.75 * std.math.pow(T, 3.0 / std.math.pi, 1.0 / 3.0);
            const cx_potential = -std.math.pow(T, 3.0 / std.math.pi, 1.0 / 3.0);

            const rho_third = std.math.pow(T, rho, 1.0 / 3.0);

            return .{
                .eps_x = cx_energy * rho_third,
                .v_x = cx_potential * rho_third
            };
        }
    };
}

/// The Local Density Approximation (LDA) functional for correlation energy by Chachiyo. https://doi.org/10.1063/1.4958669
pub fn ChachiyoCorrelation(comptime T: type) type {
    return struct {

        /// Evaluate function.
        pub fn evaluate(_: @This(), rho: T) struct {eps_c: T, v_c: T} {
            if (rho <= 1e-12) return .{.eps_c = 0, .v_c = 0};

            const a: T = (std.math.ln2 - 1) / (2 * std.math.pi * std.math.pi); const b: T = 20.4562557;

            const rs = std.math.pow(T, 3.0 / (4 * rho * std.math.pi), 1.0 / 3.0);

            const eps_c = a * std.math.log(T, std.math.e, 1 + b / rs + b / (rs * rs));
            const v_c = eps_c + a * b * (rs + 2) / (3 * (rs * rs + b * rs + b));

            return .{
                .eps_c = eps_c,
                .v_c = v_c
            };
        }
    };
}
