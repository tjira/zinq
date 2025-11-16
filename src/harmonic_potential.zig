//! File with harmonic potential struct and functions.

const std = @import("std");

const error_handling = @import("error_handling.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const throw = error_handling.throw;

/// Struct holding parameters for the multidimensional harmonic potential.
pub fn HarmonicPotential(comptime T: type) type {
    return struct {
        k: []const T = &[_]T{1},

        /// Diabatic potential evaluator.
        pub fn evaluateDiabatic(self: @This(), U: *RealMatrix(T), position: RealVector(T), _: T) void {
            U.ptr(0, 0).* = 0;

            for (0..self.k.len) |i| {
                U.ptr(0, 0).* += 0.5 * self.k[i] * position.at(i) * position.at(i);
            }
        }

        /// Diabatic potential matrix element evaluator.
        pub fn evaluateDiabaticElement(self: @This(), i: usize, j: usize, position: RealVector(T), _: T) !T {
            if (i >= 1 or j >= 1) return throw(T, "INVALID INDEX WHEN EVALUATING DIABATIC MATRIX ELEMENT", .{});

            var value: T = 0;

            for (0..self.k.len) |k| {
                value += 0.5 * self.k[k] * position.at(k) * position.at(k);
            }

            return value;
        }
    };
}
