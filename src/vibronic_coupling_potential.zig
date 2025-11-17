//! File with vibronic potential structs and functions.

const std = @import("std");

const error_handling = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const throw = error_handling.throw;

const AU2EV = global_variables.AU2EV;
const EV2RCM = global_variables.EV2RCM;

/// Struct holding parameters for the vibronic coupling potential.
pub fn VibronicCouplingPotential(comptime T: type) type {
    return struct {
        const MorsePotential = struct {
            d: []const []const T,
            a: []const []const T,
            q: []const []const T,
            e: []const []const T
        };

        energy: []const T, omega: []const T,

        morse_potential: MorsePotential,
        diagonal_linear: []const []const T,
        diagonal_quadratic: []const []const T,
        diagonal_quartic: []const []const T,
        nondiagonal_linear: []const []const T,

        /// Diabatic potential matrix evaluator.
        pub fn evaluateDiabatic(self: @This(), U: *RealMatrix(T), position: RealVector(T), time: T) !void {
            for (0..self.nstate()) |i| for (i..self.nstate()) |j| {
                U.ptr(i, j).* = try self.evaluateDiabaticElement(i, j, position, time); U.ptr(j, i).* = U.at(i, j);
            };
        }

        /// Diabatic potential matrix element evaluator.
        pub fn evaluateDiabaticElement(self: @This(), i: usize, j: usize, position: RealVector(T), _: T) !T {
            if (i >= self.nstate() or j >= self.nstate()) return throw(T, "INVALID INDEX WHEN EVALUATING DIABATIC MATRIX ELEMENT", .{});

            var value = if (i == j) self.energy[i] else 0;

            for (0..position.len) |k| {

                const qi = position.at(k) * std.math.sqrt(self.omega[k] / EV2RCM / AU2EV);

                if (i == j) {

                    const d = self.morse_potential.d[i][k];
                    const a = self.morse_potential.a[i][k];
                    const q = self.morse_potential.q[i][k];
                    const e = self.morse_potential.e[i][k];

                    const harmonic = d == 0 and a == 0 and q == 0 and e == 0;

                    value += if (harmonic) self.omega[k] * qi * qi / (2 * EV2RCM) else d * std.math.pow(T, std.math.exp(-a * (qi - q)) - 1, 2) + e;

                    value += self.diagonal_linear[i][k] * qi;
                    value += self.diagonal_quadratic[i][k] * qi * qi / 2;
                    value += self.diagonal_quartic[i][k] * qi * qi * qi * qi / 24;
                }

                if (i != j) {
                    value += self.nondiagonal_linear[i * (2 * self.nstate() - i - 1) / 2 + (j - i - 1)][k] * qi;
                }
            }

            return value / AU2EV;
        }

        /// Dimension getter.
        pub fn ndim(self: @This()) usize {
            return self.omega.len;
        }

        /// State getter.
        pub fn nstate(self: @This()) usize {
            return self.energy.len;
        }
    };
}
