//! File with Tully's potential structs and functions.

const std = @import("std");

const math_functions = @import("math_functions.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const sgn = math_functions.sgn;

/// Struct holding parameters for the first Tully's potential.
pub fn TullyPotential1(comptime T: type) type {
    return struct {
        A: T = 0.01, B: T = 1.6, C: T = 0.005, D: T = 1,

        /// Diabatic potential matrix evaluator.
        pub fn evaluateDiabatic(self: @This(), U: *RealMatrix(T), position: RealVector(T), time: T) void {
            _ = time;

            U.ptr(0, 0).* = sgn(position.at(0)) * self.A * (1 - std.math.exp(-sgn(position.at(0)) * self.B * position.at(0)));
            U.ptr(0, 1).* = self.C * std.math.exp(-self.D * position.at(0) * position.at(0));
            U.ptr(1, 0).* = U.at(0, 1);
            U.ptr(1, 1).* = -U.at(0, 0);
        }
    };
}
