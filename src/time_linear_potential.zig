//! File with time-linear potential model used in Landau-Zener problem.

const std = @import("std");

const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

/// Struct holding parameters for the time-linear potential.
pub fn TimeLinearPotential(comptime T: type) type {
    return struct {
        coupling: T = 2,
        slope: T = 10,

        /// Diabatic potential evaluator.
        pub fn evaluateDiabatic(self: @This(), U: *RealMatrix(T), _: RealVector(T), time: T) void {
            U.ptr(0, 0).* = self.slope * (time - self.slope);
            U.ptr(0, 1).* = self.coupling;
            U.ptr(1, 0).* = self.coupling;
            U.ptr(1, 1).* = -U.at(0, 0);
        }
    };
}
