//! File with harmonic potential struct and functions.

const std = @import("std");

const classical_particle = @import("classical_particle.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const ClassicalParticle = classical_particle.ClassicalParticle;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

/// Struct holding parameters for the multidimensional harmonic potential.
pub fn HarmonicPotential(comptime T: type) type {
    return struct {
        k: []const T = &[_]T{1},

        /// Diabatic potential evaluator.
        pub fn evaluateDiabatic(self: HarmonicPotential(T), U: *RealMatrix(T), position: RealVector(T), time: T) void {
            _ = time; U.ptr(0, 0).* = 0;

            for (0..self.k.len) |i| {
                U.ptr(0, 0).* += 0.5 * self.k[i] * position.at(i) * position.at(i);
            }
        }
    };
}
