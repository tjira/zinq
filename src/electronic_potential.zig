//! General union type for different electronic potential models.

const std = @import("std");

const classical_particle = @import("classical_particle.zig");
const harmonic_potential = @import("harmonic_potential.zig");
const linear_algebra = @import("linear_algebra.zig");
const real_matrix = @import("real_matrix.zig");
const tully_potential = @import("tully_potential.zig");

const ClassicalParticle = classical_particle.ClassicalParticle;
const HarmonicPotential = harmonic_potential.HarmonicPotential;
const RealMatrix = real_matrix.RealMatrix;
const TullyPotential1 = tully_potential.TullyPotential1;

const eighReal = linear_algebra.eighReal;

/// Electronic potential mode union.
pub fn ElectronicPotential(comptime T: type) type {
    return union(enum) {
        harmonic: HarmonicPotential(T),
        tully_1: TullyPotential1(T),

        /// Evaluate the dabatic potential energy matrix at given system state and time.
        pub fn evaluateDiabatic(self: @This(), diabatic_potential: *RealMatrix(T), system: ClassicalParticle(T), time: T) void {
            switch (self) {
                inline else => |field| field.evaluateDiabatic(diabatic_potential, system, time)
            }
        }

        /// Evaluate the adabatic potential energy matrix at given system state and time.
        pub fn evaluateAdiabatic(self: @This(), adiabatic_potential: *RealMatrix(T), system: ClassicalParticle(T), time: T) void {
            self.evaluateDiabatic(adiabatic_potential, system, time);
        }

        /// Getter for number of dimensions.
        pub fn ndim(self: @This()) usize {
            return switch (self) {
                .harmonic => |field| field.k.len,
                .tully_1 => 1
            };
        }

        /// Getter for number of electronic states.
        pub fn nstate(self: @This()) usize {
            return switch (self) {
                .harmonic => |_| 1,
                .tully_1 => |_| 2
            };
        }
    };
}
