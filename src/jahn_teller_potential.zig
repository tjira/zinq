//! File with the Jahn-Teller potential implementation.

const std = @import("std");

const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

/// Struct holding parameters for the Jahn-Teller potential.
pub fn JahnTellerPotential(comptime T: type) type {
    return struct {
        k: T = 1, g: T = 1,

        /// Diabatic potential matrix evaluator.
        pub fn evaluateDiabatic(self: @This(), U: *RealMatrix(T), position: RealVector(T), _: T) void {
            U.ptr(0, 0).* = self.g * position.at(0) + 0.5 * self.k * (position.at(0) * position.at(0) + position.at(1) * position.at(1));
            U.ptr(0, 1).* = self.g * position.at(1);
            U.ptr(1, 0).* = self.g * position.at(1);
            U.ptr(1, 1).* = 0.5 * self.k * (position.at(0) * position.at(0) + position.at(1) * position.at(1)) - self.g * position.at(0);
        }
    };
}
