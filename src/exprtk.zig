const std = @import("std");

const exprtk = @import("cimport.zig").exprtk;

const primType = @import("value.zig").primType;

pub fn Expression(comptime T: type) type {
    return struct {
        ptr: *exprtk.Expression,

        pub fn init(expr: []const u8, nvar: usize) !@This() {
            if (comptime primType(T) != f64) @compileError("EXPRTK NOW ONLY SUPPORTS F64 NUMBERS");

            const ptr = exprtk.exprtk_init(expr.ptr, expr.len, nvar) orelse return error.InitializationFailed;
            errdefer exprtk.exprtk_deinit(ptr);

            return .{ .ptr = ptr };
        }

        pub fn deinit(self: @This()) void {
            exprtk.exprtk_deinit(self.ptr);
        }

        pub fn evaluate_d0(self: @This(), vars: []const T, time: T) T {
            std.debug.assert(vars.len == self.ndim());

            return exprtk.exprtk_evaluate_d0(self.ptr, vars.ptr, time);
        }

        pub fn evaluate_d1(self: @This(), vars: []const T, time: T, idx: usize) T {
            std.debug.assert(vars.len == self.ndim());

            return exprtk.exprtk_evaluate_d1(self.ptr, vars.ptr, time, idx);
        }

        pub fn ndim(self: @This()) usize {
            return exprtk.exprtk_nvar(self.ptr);
        }
    };
}
