//! Wrapper for the libxc library, which provides exchange-correlation functionals for density functional theory (DFT) calculations.

const std = @import("std");

const config = @import("config");

const xc = if (config.use_xc) @cImport(@cInclude("xc.h")) else struct {};

const real_vector = @import("real_vector.zig");

const RealVector = real_vector.RealVector;

/// Evaluate functional using libxc.
pub fn computeExchangeCorrelationArrayLibxc(comptime T: type, result: *RealVector(T), name: ?[:0]const u8, rho: RealVector(T), comptime derivative: u32) !void {
    if (name == null) return result.zero();

    const func_id = xc.xc_functional_get_number(name.?.ptr); var func: xc.xc_func_type = undefined;

    if (xc.xc_func_init(&func, func_id, xc.XC_UNPOLARIZED) != 0) {
        std.log.err("FUNCTIONAL '{s}' NOT FOUND IN LIBXC\n", .{name.?}); return error.InputError;
    }

    defer xc.xc_func_end(&func);

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
}

/// Get the functional family from the libxc functional name.
pub fn getFunctionalFamilyLibxc(name: [:0]const u8) !enum{lda, gga, mgga} {
    if (comptime !config.use_xc) @compileError("LIBXC SUPPORT IS DISABLED IN CONFIG");

    const func_id = xc.xc_functional_get_number(name.ptr); var func: xc.xc_func_type = undefined;

    if (xc.xc_func_init(&func, func_id, xc.XC_UNPOLARIZED) != 0) {
        std.log.err("FUNCTIONAL '{s}' NOT FOUND IN LIBXC", .{name}); return error.InputError;
    }

    defer xc.xc_func_end(&func);

    switch (func.info.*.family) {
        xc.XC_FAMILY_LDA, xc.XC_FAMILY_HYB_LDA => return .lda,
        xc.XC_FAMILY_GGA, xc.XC_FAMILY_HYB_GGA => return .gga,
        xc.XC_FAMILY_MGGA, xc.XC_FAMILY_HYB_MGGA => return .mgga,
        else => {std.log.err("CANNOT IDENTIFY FUNCTIONAL FAMILY", .{}); return error.InputError;}
    }
}

/// Get the functional kind from the libxc functional name.
pub fn getFunctionalKindLibxc(name: [:0]const u8) !enum{exchange, correlation, exchange_correlation, kinetic} {
    if (comptime !config.use_xc) @compileError("LIBXC SUPPORT IS DISABLED IN CONFIG");

    const func_id = xc.xc_functional_get_number(name.ptr); var func: xc.xc_func_type = undefined;

    if (xc.xc_func_init(&func, func_id, xc.XC_UNPOLARIZED) != 0) {
        std.log.err("FUNCTIONAL '{s}' NOT FOUND IN LIBXC", .{name}); return error.InputError;
    }

    defer xc.xc_func_end(&func);

    switch (func.info.*.kind) {
        xc.XC_EXCHANGE => return .exchange,
        xc.XC_CORRELATION => return .correlation,
        xc.XC_EXCHANGE_CORRELATION => return .exchange_correlation,
        xc.XC_KINETIC => return .kinetic,
        else => {std.log.err("CANNOT IDENTIFY FUNCTIONAL KIND", .{}); return error.InputError;}
    }
}
