//! File with time-linear potential model used in Landau-Zener problem.

const std = @import("std");

const error_handling = @import("error_handling.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const throw = error_handling.throw;

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

        /// Diabatic potential matrix element evaluator.
        pub fn evaluateDiabaticElement(self: @This(), i: usize, j: usize, _: RealVector(T), time: T) !T {
            if (i == 0 and j == 0) return self.slope * (time - self.slope);
            if (i + j == 1) return self.coupling;
            if (i == 1 and j == 1) return self.slope * (self.slope - time);

            return throw(T, "INVALID INDEX WHEN EVALUATING DIABATIC MATRIX ELEMENT", .{});
        }
    };
}
