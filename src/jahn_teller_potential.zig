//! File with the Jahn-Teller potential implementation.

const std = @import("std");

const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

/// Struct holding parameters for the Jahn-Teller potential.
pub fn JahnTellerPotential(comptime T: type) type {
    return struct {
        k: T = 1, g: T = 1, d: usize = 2, abs_coupling: bool = false, coupling_sum: usize = 1,

        /// Diabatic potential matrix evaluator.
        pub fn evaluateDiabatic(self: @This(), U: *RealMatrix(T), position: RealVector(T), time: T) !void {
            U.ptr(0, 0).* = try self.evaluateDiabaticElementComptime(0, 0, position, time);
            U.ptr(0, 1).* = try self.evaluateDiabaticElementComptime(0, 1, position, time);
            U.ptr(1, 0).* = try self.evaluateDiabaticElementComptime(1, 0, position, time);
            U.ptr(1, 1).* = try self.evaluateDiabaticElementComptime(1, 1, position, time);
        }

        /// Diabatic potential matrix element evaluator.
        pub fn evaluateDiabaticElement(self: @This(), i: usize, j: usize, position: RealVector(T), time: T) !T {
            return switch (i + j) {
                0 => try self.evaluateDiabaticElementComptime(0, 0, position, time),
                1 => try self.evaluateDiabaticElementComptime(0, 1, position, time),
                2 => try self.evaluateDiabaticElementComptime(1, 1, position, time),
                else => {

                    std.log.err("INVALID INDEX ({d}, {d}) WHEN EVALUATING DIABATIC MATRIX ELEMENT, THE POTENTIAL MATRIX IS {d}X{d}", .{i, j, 2, 2});

                    return error.ProgrammingError;
                }
            };
        }

        /// Comptime potential matrix element evaluator.
        pub fn evaluateDiabaticElementComptime(self: @This(), comptime i: usize, comptime j: usize, position: RealVector(T), _: T) !T {
            return switch (i + j) {
                0 => {

                    var sumsq: T = 0; for (position.data) |q| sumsq += q * q;

                    return self.g * position.at(0) + 0.5 * self.k * sumsq;
                },
                1 => {

                    if (self.coupling_sum > position.len - 1) {

                        std.log.err("COUPLING SUM HAS TO BE LESS THAN OR EQUAL TO THE NUMBER OF NUCLEAR COORDINATES MINUS ONE", .{});

                        return error.InputError;
                    }

                    var sum: T = 0;

                    for (0..self.coupling_sum) |k| sum += if (self.abs_coupling) @abs(position.at(1 + k)) else position.at(1 + k);

                    return self.g * sum;
                },
                2 => {
                    
                    var sumsq: T = 0; for (position.data) |q| sumsq += q * q;

                    return 0.5 * self.k * sumsq - self.g * position.at(0);
                },
                else => @compileError("INVALID INDEX WHEN EVALUATING DIABATIC MATRIX ELEMENT")
            };
        }
    };
}
