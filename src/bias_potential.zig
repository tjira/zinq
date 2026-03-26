//! Bias potential module.

const std = @import("std");

const real_matrix = @import("real_matrix.zig");

const RealMatrix = real_matrix.RealMatrix;

/// One dimensional bias potential struct.
pub fn BiasPotential(comptime T: type) type {
    return struct {
        pub const Function = struct {
            pub const Harmonic = struct {
                k: T
            };
        };
        pub const Variable = struct {
            pub const PotentialEnergy = struct {
                value: T
            };
            pub const PotentialEnergyDifference = struct {
                value: T,
                states: [2]u32
            };
        };

        variable: union(enum) {
            potential_energy: Variable.PotentialEnergy,
            potential_energy_difference: Variable.PotentialEnergyDifference,
        },

        function: union(enum) {
            harmonic: Function.Harmonic,
        },

        /// Evaluate the bias potential at given position.
        pub fn evaluate(self: @This(), variable: T) T {
            const value = switch (self.variable) {
                .potential_energy => self.variable.potential_energy.value,
                .potential_energy_difference => self.variable.potential_energy_difference.value,
            };

            return switch (self.function) {
                .harmonic => |field| 0.5 * field.k * (variable - value) * (variable - value),
            };
        }

        /// Evaluate the force from the bias potential for the specified coordinate index.
        pub fn evaluateForce(self: @This(), variable: T) T {
            const value = switch (self.variable) {
                .potential_energy => self.variable.potential_energy.value,
                .potential_energy_difference => self.variable.potential_energy_difference.value,
            };

            return switch (self.function) {
                .harmonic => |field| -field.k * (variable - value),
            };
        }

        /// Evaluate the force from the bias potential for the specified system state and coordinate index.
        pub fn force(self: @This(), adiabatic_potential: RealMatrix(T), state: usize, _: usize) !T {
            const variable = switch (self.variable) {
                .potential_energy => adiabatic_potential.at(state, state),
                .potential_energy_difference => |field| blk: {

                    const state1 = field.states[0];
                    const state2 = field.states[1];

                    if (state1 >= adiabatic_potential.rows or state2 >= adiabatic_potential.rows) {

                        std.log.err("STATE OUT OF BOUNDS IN BIAS POTENTIAL", .{});

                        return error.InvalidInput;
                    }

                    break :blk adiabatic_potential.at(state2, state2) - adiabatic_potential.at(state1, state1);
                }
            };

            return self.evaluateForce(variable);
         }
    };
}
