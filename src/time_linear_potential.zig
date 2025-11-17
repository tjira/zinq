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
        pub fn evaluateDiabatic(self: @This(), U: *RealMatrix(T), position: RealVector(T), time: T) !void {
            U.ptr(0, 0).* = self.evaluateDiabaticElementComptime(0, 0, position, time);
            U.ptr(0, 1).* = self.evaluateDiabaticElementComptime(0, 1, position, time);
            U.ptr(1, 0).* = self.evaluateDiabaticElementComptime(1, 0, position, time);
            U.ptr(1, 1).* = self.evaluateDiabaticElementComptime(1, 1, position, time);
        }

        /// Diabatic potential matrix element evaluator.
        pub fn evaluateDiabaticElement(self: @This(), i: usize, j: usize, position: RealVector(T), time: T) !T {
            return switch (i + j) {
                0 => self.evaluateDiabaticElementComptime(0, 0, position, time),
                1 => self.evaluateDiabaticElementComptime(0, 1, position, time),
                2 => self.evaluateDiabaticElementComptime(1, 1, position, time),
                else => throw(T, "INVALID INDEX WHEN EVALUATING DIABATIC MATRIX ELEMENT", .{})
            };
        }

        /// Comptime diabatic matrix element evaluator.
        pub fn evaluateDiabaticElementComptime(self: @This(), comptime i: usize, comptime j: usize, _: RealVector(T), time: T) T {
            return switch (i + j) {
                0 => return self.slope * (time - self.slope),
                1 => return self.coupling,
                2 => return self.slope * (self.slope - time),
                else => @compileError("INVALID INDEX WHEN EVALUATING DIABATIC MATRIX ELEMENT"),
            };
        }
    };
}
