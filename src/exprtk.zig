//! Binding to ExprTk, a mathematical expression parser and evaluation engine for high-performance computations.

const std = @import("std");

const exprtk = @import("cimport.zig").exprtk;

const primType = @import("value.zig").primType;

/// Returns an expression type parser that parses and evaluates mathematical expression strings.
pub fn Expression(comptime T: type) type {
    return struct {
        ptr: *exprtk.Expression,

        /// Parses and compiles a mathematical expression string given the number of independent variables.
        pub fn init(expr: []const u8, nvar: usize) !@This() {
            if (comptime primType(T) != f64) @compileError("EXPRTK NOW ONLY SUPPORTS F64 NUMBERS");

            const ptr = exprtk.exprtk_init(expr.ptr, expr.len, nvar) orelse return error.InitializationFailed;
            errdefer exprtk.exprtk_deinit(ptr);

            return .{ .ptr = ptr };
        }

        /// Frees the allocated resources associated with the compiled mathematical expression.
        pub fn deinit(self: @This()) void {
            exprtk.exprtk_deinit(self.ptr);
        }

        /// Evaluates the zeroth derivative (value) of the function for the given variables and time.
        pub fn evaluate_d0(self: @This(), vars: []const T, time: T) T {
            std.debug.assert(vars.len == self.ndim());

            return exprtk.exprtk_evaluate_d0(self.ptr, vars.ptr, time);
        }

        /// Evaluates the first derivative of the function with respect to the variable at index idx.
        pub fn evaluate_d1(self: @This(), vars: []const T, time: T, idx: usize) T {
            std.debug.assert(vars.len == self.ndim());

            return exprtk.exprtk_evaluate_d1(self.ptr, vars.ptr, time, idx);
        }

        /// Returns the number of independent variables defined in the mathematical expression.
        pub fn ndim(self: @This()) usize {
            return exprtk.exprtk_nvar(self.ptr);
        }
    };
}
