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
        pub fn evaluateDiabatic(self: @This(), U: *RealMatrix(T), position: RealVector(T), _: T) void {
            U.fill(0); for (0..self.nstate()) |i| U.ptr(i, i).* = self.energy[i];

            for (0..position.len) |i| {

                const qi = position.at(i) * std.math.sqrt(self.omega[i] / EV2RCM / AU2EV);

                for (0..self.nstate()) |j| {

                    const d = self.morse_potential.d[j][i];
                    const a = self.morse_potential.a[j][i];
                    const q = self.morse_potential.q[j][i];
                    const e = self.morse_potential.e[j][i];

                    const harmonic = d == 0 and a == 0 and q == 0 and e == 0;

                    U.ptr(j, j).* += if (harmonic) self.omega[i] * qi * qi / (2 * EV2RCM) else d * std.math.pow(T, std.math.exp(-a * (qi - q)) - 1, 2) + e;
                }

                for (0..self.nstate()) |j| U.ptr(j, j).* += self.diagonal_linear[j][i] * qi;
                for (0..self.nstate()) |j| U.ptr(j, j).* += self.diagonal_quadratic[j][i] * qi * qi / 2;
                for (0..self.nstate()) |j| U.ptr(j, j).* += self.diagonal_quartic[j][i] * qi * qi * qi * qi / 24;

                for (0..U.rows) |j| for (j + 1..U.cols) |k| {
                    U.ptr(j, k).* += self.nondiagonal_linear[j * (2 * self.nstate() - j - 1) / 2 + (k - j - 1)][i] * qi;
                };
            }

            for (0..self.nstate()) |i| for (i..self.nstate()) |j| {
                U.ptr(i, j).* /= AU2EV; U.ptr(j, i).* = U.at(i, j);
            };
        }

        /// Diabatic potential matrix element evaluator.
        pub fn evaluateDiabaticElement(self: @This(), i: usize, j: usize, position: RealVector(T), _: T) !T {
            if (i >= self.nstate() or j >= self.nstate()) return throw(T, "ELEMENT EVALUATION FOR VIBRONIC COUPLING POTENTIAL NOT IMPLEMENTED", .{});

            var value = if (i == j) self.energy[i] else 0;

            for (0..position.len) |k| {

                const qi = position.at(i) * std.math.sqrt(self.omega[k] / EV2RCM / AU2EV);

                if (i == j) {

                    for (0..self.nstate()) |l| {

                        const d = self.morse_potential.d[l][k];
                        const a = self.morse_potential.a[l][k];
                        const q = self.morse_potential.q[l][k];
                        const e = self.morse_potential.e[l][k];

                        const harmonic = d == 0 and a == 0 and q == 0 and e == 0;

                        value += if (harmonic) self.omega[k] * qi * qi / (2 * EV2RCM) else d * std.math.pow(T, std.math.exp(-a * (qi - q)) - 1, 2) + e;
                    }

                    for (0..self.nstate()) |l| value += self.diagonal_linear[l][k] * qi;
                    for (0..self.nstate()) |l| value += self.diagonal_quadratic[l][k] * qi * qi / 2;
                    for (0..self.nstate()) |l| value += self.diagonal_quartic[l][k] * qi * qi * qi * qi / 24;
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
