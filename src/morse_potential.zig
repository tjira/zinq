//! File with Morse potential struct and functions.

const std = @import("std");

const error_handling = @import("error_handling.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const throw = error_handling.throw;

/// Struct holding parameters for the 1D Morse potential.
pub fn MorsePotential(comptime T: type) type {
    return struct {
        De: T = 0.5,
        a: T = 1,
        r: T = 0,

        /// Diabatic potential evaluator.
        pub fn evaluateDiabatic(self: @This(), U: *RealMatrix(T), position: RealVector(T), time: T) void {
            U.ptr(0, 0).* = self.evaluateDiabaticElementComptime(0, 0, position, time);
        }

        /// Diabatic potential matrix element evaluator.
        pub fn evaluateDiabaticElement(self: @This(), i: usize, j: usize, position: RealVector(T), time: T) !T {
            if (i >= 1 or j >= 1) return throw(T, "INVALID INDEX WHEN EVALUATING DIABATIC MATRIX ELEMENT", .{});

            return self.evaluateDiabaticElementComptime(0, 0, position, time);
        }

        /// Comptime potential matrix element evaluator.
        pub fn evaluateDiabaticElementComptime(self: @This(), comptime i: usize, comptime j: usize, position: RealVector(T), _: T) T {
            if (i >= 1 or j >= 1) @compileError("INVALID INDEX WHEN EVALUATING DIABATIC MATRIX ELEMENT");

            return self.De * std.math.pow(T, 1 - std.math.exp(-self.a * (position.at(0) - self.r)), 2);
        }
    };
}
