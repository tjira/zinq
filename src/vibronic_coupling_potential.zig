//! File with vibronic potential structs and functions.

const std = @import("std");

const global_variables = @import("global_variables.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

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
