//! General union type for different electronic potential models.

const std = @import("std");

const classical_particle = @import("classical_particle.zig");
const harmonic_potential = @import("harmonic_potential.zig");
const real_matrix = @import("real_matrix.zig");
const tully_potential = @import("tully_potential.zig");

const ClassicalParticle = classical_particle.ClassicalParticle;
const HarmonicPotential = harmonic_potential.HarmonicPotential;
const RealMatrix = real_matrix.RealMatrix;
const TullyPotential1 = tully_potential.TullyPotential1;

/// Electronic potential modeunion.
pub fn ElectronicPotential(comptime T: type) type {
    return union(enum) {
        harmonic: HarmonicPotential(T),
        tully_1: TullyPotential1(T),

        /// Evaluate the potential energy matrix at given system state and time.
        pub fn eval(self: @This(), U: *RealMatrix(T), system: ClassicalParticle(T), time: T) !void {
            switch (self) {
                inline else => |field| return try field.eval(U, system, time),
            }
        }

        /// Getter for number of dimensions.
        pub fn ndim(self: @This()) usize {
            return switch (self) {
                .harmonic => |field| field.k.len,
                inline else => |_| 1,
            };
        }

        /// Getter for number of electronic states.
        pub fn nstate(self: @This()) usize {
            return switch (self) {
                .tully_1 => |_| 2,
                inline else => |_| 1,
            };
        }
    };
}
