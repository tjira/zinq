//! Implmentation of DFT functionals.

const std = @import("std");

/// Enumeration of available correlation functionals.
pub fn CorrelationFunctional(comptime T: type) type {
    return union(enum) {
        chachiyo: ChachiyoCorrelation(T), Chachiyo: ChachiyoCorrelation(T),
        vwn5: VWN5Correlation(T), VWN5: VWN5Correlation(T)
    };
}

/// Enumeration of available exchange functionals.
pub fn ExchangeFunctional(comptime T: type) type {
    return union(enum) {
        slater: SlaterExchange(T), Slater: SlaterExchange(T)
    };
}

/// Compute the exchange-correlation energy density and potential for a given electron density `rho` using the specified DFT functional.
pub fn computeExchangeCorrelation(comptime T: type, exchange: ?ExchangeFunctional(T), correlation: ?CorrelationFunctional(T), rho: T) struct {T, T} {
    const eps_x, const v_x = if (exchange) |ex| switch (ex) {
        inline else => |x| x.evaluate(rho)
    } else .{0, 0};

    const eps_c, const v_c = if (correlation) |corr| switch (corr) {
        inline else => |c| c.evaluate(rho)
    } else .{0, 0};

    return .{eps_x + eps_c, v_x + v_c};
}

/// The Local Density Approximation (LDA) functional for exchange-correlation energy.
pub fn SlaterExchange(comptime T: type) type {
    return struct {

        /// Evaluate function.
        pub fn evaluate(_: @This(), rho: T) struct {T, T} {
            if (rho <= 1e-12) return .{0, 0};

            const cx_energy = -0.75 * std.math.pow(T, 3.0 / std.math.pi, 1.0 / 3.0);
            const cx_potential = -std.math.pow(T, 3.0 / std.math.pi, 1.0 / 3.0);

            const rho_third = std.math.pow(T, rho, 1.0 / 3.0);

            return .{
                cx_energy * rho_third,
                cx_potential * rho_third
            };
        }
    };
}

/// The Local Density Approximation (LDA) functional for correlation energy by Chachiyo. https://doi.org/10.1063/1.4958669
pub fn ChachiyoCorrelation(comptime T: type) type {
    return struct {

        /// Evaluate function.
        pub fn evaluate(_: @This(), rho: T) struct {T, T} {
            if (rho <= 1e-12) return .{0, 0};

            const a: T = (std.math.ln2 - 1) / (2 * std.math.pi * std.math.pi); const b: T = 20.4562557;

            const rs = std.math.pow(T, 3 / (4 * rho * std.math.pi), 1.0 / 3.0);

            const eps = a * std.math.log(T, std.math.e, 1 + b / rs + b / (rs * rs));
            const v = eps + a * b * (rs + 2) / (3 * (rs * rs + b * rs + b));

            return .{eps, v};
        }
    };
}

/// The VWN5 Local Density Approximation (LDA) functional for correlation energy. https://doi.org/10.1139/p80-159
pub fn VWN5Correlation(comptime T: type) type {
    return struct {
        /// Evaluate function.
        pub fn evaluate(_: @This(), rho: T) struct {T, T} {
            if (rho <= 1e-12) return .{0, 0};

            const A: T = 0.0310907; const x0: T = -0.10498; const b: T = 3.72744; const c: T = 12.9352;

            const rs = std.math.pow(T, 3 / (4 * std.math.pi * rho), 1.0 / 3.0);

            const x = std.math.sqrt(rs); const Q = std.math.sqrt(4 * c - b * b);

            const Xx = x * x + b * x + c; const Xx0 = x0 * x0 + b * x0 + c;

            const term1 = std.math.log(T, std.math.e, x * x / Xx);
            const term2 = 2 * b / Q * std.math.atan(Q / (2 * x + b));
            const term3 = std.math.log(T, std.math.e, (x - x0) * (x - x0) / Xx);
            const term4 = 2 * (b + 2 * x0) / Q * std.math.atan(Q / (2 * x + b));
            
            const eps = A * (term1 + term2 - b * x0 / Xx0 * (term3 + term4));

            const v = eps - A * x * (2 / x - 2 * (x + b) / Xx - b * x0 / Xx0 * (2 / (x - x0) - 2 * (x + b + x0) / Xx)) / 6;

            return .{eps, v};
        }
    };
}
