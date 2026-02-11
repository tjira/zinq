//! File with Richardson's avoided crossing potential.

const std = @import("std");

const error_handling = @import("error_handling.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const throw = error_handling.throw;

/// Struct holding parameters for the Richardson's avoided crossing potential, along with methods to evaluate the potential matrix and its elements.
pub fn AvoidedCrossingPotential(comptime T: type) type {
    return struct {
        m: T = 1,
        omega: T = 1,
        kappa: T = 3,
        epsilon: T = 4,
        delta: T = 2,

        /// Diabatic potential matrix evaluator.
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

        /// Comptime potential matrix element evaluator.
        pub fn evaluateDiabaticElementComptime(self: @This(), comptime i: usize, comptime j: usize, position: RealVector(T), _: T) T {
            return switch (i + j) {
                0 => 0.5 * self.m * self.omega * self.omega * position.at(0) * position.at(0) + self.kappa * position.at(0) + self.epsilon,
                1 => self.delta,
                2 => 0.5 * self.m * self.omega * self.omega * position.at(0) * position.at(0) - self.kappa * position.at(0) - self.epsilon,
                else => @compileError("INVALID INDEX WHEN EVALUATING DIABATIC MATRIX ELEMENT"),
            };
        }
    };
}
