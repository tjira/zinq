//! File with Morse potential struct and functions.

const std = @import("std");

const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

/// Struct holding parameters for the 1D Morse potential.
pub fn MorsePotential(comptime T: type) type {
    return struct {
        De: T = 0.5,
        a: T = 1,
        r: T = 0,

        /// Diabatic potential evaluator.
        pub fn evaluateDiabatic(self: @This(), U: *RealMatrix(T), position: RealVector(T), _: T) void {
            U.ptr(0, 0).* = self.De * std.math.pow(T, 1 - std.math.exp(-self.a * (position.at(0) - self.r)), 2);
        }
    };
}
