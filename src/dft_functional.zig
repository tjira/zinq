//! Implmentation of DFT functionals.

const std = @import("std");

const config = @import("config");

const xc = if (config.use_xc) @cImport(@cInclude("xc.h")) else struct {};

const real_vector = @import("real_vector.zig");

const RealVector = real_vector.RealVector;

/// Density funcitional input.
pub fn DensityFunctional(comptime T: type) type {
    return struct {
        exchange: ?[:0]const u8 = "LDA_X",
        correlation: ?[:0]const u8 = "LDA_C_VWN",
        params: ?[]const T = null
    };
}

/// Enumeration of available exchange functionals.
pub fn ExchangeFunctional(comptime T: type) type {
    return union(enum) {
        slater: SlaterExchange(T),
    };
}

/// Enumeration of available correlation functionals.
pub fn CorrelationFunctional(comptime T: type) type {
    return union(enum) {
        chachiyo: ChachiyoCorrelation(T),
        vwn5: VWN5Correlation(T),
    };
}

/// Get the functional kind.
pub fn getFunctionalKind(name: [:0]const u8) !enum{lda, gga, mgga} {
    if (comptime config.use_xc) {

        const func_id = xc.xc_functional_get_number(name.ptr); var func: xc.xc_func_type = undefined;

        if (xc.xc_func_init(&func, func_id, xc.XC_UNPOLARIZED) != 0) {
            std.log.err("FUNCTIONAL '{s}' NOT FOUND IN LIBXC\n", .{name}); return error.InputError;
        }

        switch (func.info.*.family) {
            xc.XC_FAMILY_LDA => return .lda,
            xc.XC_FAMILY_GGA, xc.XC_FAMILY_HYB_GGA => return .gga,
            xc.XC_FAMILY_MGGA, xc.XC_FAMILY_HYB_MGGA => return .mgga,
            else => {std.log.err("CANNOT IDENTIFY FUNCTIONAL TYPE", .{}); return error.InputError;}
        }
    }

    const kind_map = .{
        .{"LDA_X",          .lda},
        .{"LDA_C_CHACHIYO", .lda},
        .{"LDA_C_VWN",      .lda}
    };

    inline for (kind_map) |entry| {
        if (std.ascii.eqlIgnoreCase(name, entry[0])) {return entry[1];}
    } else {std.log.err("CANNOT DETERMINA '{s}' FUNCTIONAL KIND", .{name}); return error.InputError;}
}

/// Evaluate functional using libxc.
pub fn evaluateLibxc(comptime T: type, result: *RealVector(T), name: ?[:0]const u8, rho: RealVector(T), comptime derivative: u32) !void {
    if (name == null) return result.zero();

    const func_id = xc.xc_functional_get_number(name.?.ptr); var func: xc.xc_func_type = undefined;

    if (xc.xc_func_init(&func, func_id, xc.XC_UNPOLARIZED) != 0) {
        std.log.err("FUNCTIONAL '{s}' NOT FOUND IN LIBXC\n", .{name.?}); return error.InputError;
    }

    if (derivative == 0) switch (func.info.*.family) {
        xc.XC_FAMILY_LDA => xc.xc_lda_exc(&func, rho.len, rho.data.ptr, result.data.ptr),
        else => {std.log.err("THIS FUNCTIONAL FAMILY IS NOT YET IMPLEMENTED", .{}); return error.InputError;}
    };

    if (derivative == 1) switch (func.info.*.family) {
        xc.XC_FAMILY_LDA => xc.xc_lda_vxc(&func, rho.len, rho.data.ptr, result.data.ptr),
        else => {std.log.err("THIS FUNCTIONAL FAMILY IS NOT YET IMPLEMENTED", .{}); return error.InputError;}
    };

    if (derivative == 2) switch (func.info.*.family) {
        xc.XC_FAMILY_LDA => xc.xc_lda_fxc(&func, rho.len, rho.data.ptr, result.data.ptr),
        else => {std.log.err("THIS FUNCTIONAL FAMILY IS NOT YET IMPLEMENTED", .{}); return error.InputError;}
    };

    xc.xc_func_end(&func);
}

/// Compute the exchange-correlation energy per particle for a given electron density `rho` using the specified DFT functional.
pub fn computeExchangeCorrelationArray(comptime T: type, exc: *RealVector(T), cor: *RealVector(T), functional: DensityFunctional(T), rho: RealVector(T), comptime derivative: u32) !void {
    if (comptime config.use_xc) {

        try evaluateLibxc(T, exc, functional.exchange,    rho, derivative);
        try evaluateLibxc(T, cor, functional.correlation, rho, derivative);

        return;
    }

    var exchange: ?ExchangeFunctional(T) = null; var correlation: ?CorrelationFunctional(T) = null;

    const exchange_map = [_]struct {name: []const u8, val: ExchangeFunctional(T)}{
        .{.name = "LDA_X", .val = .{.slater = .{}}}
    };

    const correlation_map = [_]struct {name: []const u8, val: CorrelationFunctional(T)}{
        .{.name = "LDA_C_CHACHIYO", .val = .{.chachiyo = .{}}},
        .{.name = "LDA_C_VWN",      .val = .{.vwn5     = .{}}},
    };

    if (functional.exchange) |x| {
        for (exchange_map) |entry| {
            if (std.ascii.eqlIgnoreCase(x, entry.name)) {exchange = entry.val; break;}
        } else {std.log.err("EXCHANGE FUNCTIONAL '{s}' NOT FOUND", .{x}); return error.InputError;}
    }

    if (functional.correlation) |c| {
        for (correlation_map) |entry| {
            if (std.ascii.eqlIgnoreCase(c, entry.name)) {correlation = entry.val; break;}
        } else {std.log.err("CORRELATION FUNCTIONAL '{s}' NOT FOUND", .{c}); return error.InputError;}
    }

    for (0..rho.len) |i| exc.ptr(i).*, cor.ptr(i).* = computeExchangeCorrelation(T, exchange, correlation, rho.at(i), derivative);
}

/// Compute the exchange-correlation energy per particle for a given electron density `rho` using the specified DFT functional.
pub fn computeExchangeCorrelation(comptime T: type, exchange: ?ExchangeFunctional(T), correlation: ?CorrelationFunctional(T), rho: T, comptime derivative: u32) struct{T, T} {
    if (rho <= 1e-12) return .{0, 0};

    const x = if (exchange) |ex| switch (ex) {
        inline else => |x| switch (derivative) {
            0 => x.evaluateD0(rho),
            1 => x.evaluateD1(rho),
            2 => x.evaluateD2(rho),
            else => @compileError("INVALID FUNCTIONAL DERIVATIVE ORDER")
        }
    } else 0;

    const c = if (correlation) |corr| switch (corr) {
        inline else => |c| switch (derivative) {
            0 => c.evaluateD0(rho),
            1 => c.evaluateD1(rho),
            2 => c.evaluateD2(rho),
            else => @compileError("INVALID FUNCTIONAL DERIVATIVE ORDER")
        }
    } else 0;

    return .{x, c};
}

/// The Local Density Approximation (LDA) functional for exchange-correlation energy.
pub fn SlaterExchange(comptime T: type) type {
    return struct {

        /// Evaluate function.
        pub fn evaluateD0(_: @This(), rho: T) T {
            const Cx = 0.75 * std.math.pow(T, 3.0 / std.math.pi, 1.0 / 3.0);

            const eps = -Cx * std.math.pow(T, rho, 1.0 / 3.0);

            return eps;
        }

        /// First derivative of the exchange energy density with respect to the electron density `rho`.
        pub fn evaluateD1(self: @This(), rho: T) T {
            const Cx = 0.75 * std.math.pow(T, 3.0 / std.math.pi, 1.0 / 3.0);

            const v = self.evaluateD0(rho) - Cx * std.math.pow(T, rho, 1.0 / 3.0) / 3;

            return v;
        }

        /// Second derivative of the exchange energy density with respect to the electron density `rho`.
        pub fn evaluateD2(_: @This(), rho: T) T {
            const Cx = 0.75 * std.math.pow(T, 3.0 / std.math.pi, 1.0 / 3.0);

            const f = -4 * Cx * std.math.pow(T, rho, -2.0 / 3.0) / 9;

            return f;
        }
    };
}

/// The Local Density Approximation (LDA) functional for correlation energy by Chachiyo. https://doi.org/10.1063/1.4958669
pub fn ChachiyoCorrelation(comptime T: type) type {
    return struct {

        /// Evaluate function.
        pub fn evaluateD0(_: @This(), rho: T) T {
            const a: T = (std.math.ln2 - 1) / (2 * std.math.pi * std.math.pi); const b: T = 20.4562557;

            const rs = std.math.pow(T, 3 / (4 * rho * std.math.pi), 1.0 / 3.0);

            const eps = a * std.math.log(T, std.math.e, 1 + b / rs + b / (rs * rs));

            return eps;
        }

        /// First derivative of the correlation energy density with respect to the electron density `rho`.
        pub fn evaluateD1(self: @This(), rho: T) T {
            const a: T = (std.math.ln2 - 1) / (2 * std.math.pi * std.math.pi); const b: T = 20.4562557;

            const rs = std.math.pow(T, 3 / (4 * rho * std.math.pi), 1.0 / 3.0);

            const v = self.evaluateD0(rho) + a * b * (2 + rs) / (3 * (b + b * rs + rs * rs));

            return v;
        }

        /// Second derivative of the correlation energy density with respect to the electron density `rho`.
        pub fn evaluateD2(_: @This(), rho: T) T {
            const a: T = (std.math.ln2 - 1) / (2 * std.math.pi * std.math.pi); const b: T = 20.4562557;

            const rs = std.math.pow(T, 3 / (4 * rho * std.math.pi), 1.0 / 3.0);

            const f_term1 = 4 * a * b * std.math.pi * std.math.pow(T, rs, 3);
            const f_term2 = 6 * b + 10 * b * rs + 10 * rs * rs + 3 * b * rs * rs + 4 * rs * rs * rs;
            const f_term3 = 27 * (b + b * rs + rs * rs) * (b + b * rs + rs * rs);

            const f = f_term1 * f_term2 / f_term3;

            return f;
        }
    };
}

/// The VWN5 Local Density Approximation (LDA) functional for correlation energy. https://doi.org/10.1139/p80-159
pub fn VWN5Correlation(comptime T: type) type {
    return struct {

        /// Evaluate function.
        pub fn evaluateD0(_: @This(), rho: T) T {
            const A: T = 0.0310907; const x0: T = -0.10498; const b: T = 3.72744; const c: T = 12.9352;

            const rs = std.math.pow(T, 3 / (4 * std.math.pi * rho), 1.0 / 3.0);

            const x = std.math.sqrt(rs); const Q = std.math.sqrt(4 * c - b * b);

            const X = struct {pub fn eval(t: T) T {return t * t + b * t + c;}}.eval;

            const eps_term1 = std.math.log(T, std.math.e, x * x / X(x));
            const eps_term2 = 2 * b / Q * std.math.atan(Q / (2 * x + b));
            const eps_term3 = std.math.log(T, std.math.e, (x - x0) * (x - x0) / X(x));
            const eps_term4 = 2 * (b + 2 * x0) / Q * std.math.atan(Q / (2 * x + b));
            
            const eps = A * (eps_term1 + eps_term2 - b * x0 / X(x0) * (eps_term3 + eps_term4));

            return eps;
        }

        /// Evaluate function.
        pub fn evaluateD1(self: @This(), rho: T) T {
            const A: T = 0.0310907; const x0: T = -0.10498; const b: T = 3.72744; const c: T = 12.9352;

            const rs = std.math.pow(T, 3 / (4 * std.math.pi * rho), 1.0 / 3.0);

            const x = std.math.sqrt(rs);

            const X = struct {pub fn eval(t: T) T {return t * t + b * t + c;}}.eval;
            const dX = struct {pub fn eval(t: T) T {return 2 * t + b;}}.eval;

            const v_term1 = 2 * b * x0 / (x - x0);
            const v_term2 = (X(x0) * (b + dX(x)) - b * x0 * (dX(x) + dX(x0))) / X(x);

            const v = self.evaluateD0(rho) + A * (x * (v_term1 + v_term2) / X(x0) - 2) / 6;

            return v;
        }

        /// Evaluate function.
        pub fn evaluateD2(_: @This(), rho: T) T {
            const A: T = 0.0310907; const x0: T = -0.10498; const b: T = 3.72744; const c: T = 12.9352;

            const rs = std.math.pow(T, 3 / (4 * std.math.pi * rho), 1.0 / 3.0);

            const x = std.math.sqrt(rs);

            const X = struct {pub fn eval(t: T) T {return t * t + b * t + c;}}.eval;
            const dX = struct {pub fn eval(t: T) T {return 2 * t + b;}}.eval;

            const f_term1 = 2 * b * x0 * (5 * x0 - 6 * x) / ((x - x0) * (x - x0));
            const f_term2 = b * (2 * c + b * x) * X(x0) / (X(x) * X(x));
            const f_term3 = x * (b * x0 - X(x0)) * dX(x) * dX(x) / (X(x) * X(x));
            const f_term4 = b * (5 * c + x * (6 * b + 7 * x)) * x0 * dX(x0) / (X(x) * X(x));
            const f_term5 = X(x0) * (2 * x - 7 * b - 5 * dX(x)) / X(x);
            const f_term6 = b * x0 * (5 * dX(x) - 2 * x) / X(x);

            const f = -A * std.math.pi * std.math.pow(T, x, 7) * (12 / x + (f_term1 + f_term2 + f_term3 + f_term4 + f_term5 + f_term6) / X(x0)) / 27;

            return f;
        }
    };
}
