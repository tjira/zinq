//! File with Tully's potential structs and functions.

const std = @import("std");

const error_handling = @import("error_handling.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const throw = error_handling.throw;

/// Struct holding parameters for the first Tully's potential.
pub fn TullyPotential1(comptime T: type) type {
    return struct {
        A: T = 0.01,
        B: T = 1.6,
        C: T = 0.005,
        D: T = 1,

        /// Diabatic potential matrix evaluator.
        pub fn evaluateDiabatic(self: @This(), U: *RealMatrix(T), position: RealVector(T), time: T) void {
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
                0 => std.math.sign(position.at(0)) * self.A * (1 - std.math.exp(-self.B * @abs(position.at(0)))),
                1 => self.C * std.math.exp(-self.D * position.at(0) * position.at(0)),
                2 => std.math.sign(position.at(0)) * self.A * (std.math.exp(-self.B * @abs(position.at(0))) - 1),
                else => @compileError("INVALID INDEX WHEN EVALUATING DIABATIC MATRIX ELEMENT"),
            };
        }
    };
}
