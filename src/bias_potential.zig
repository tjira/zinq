//! Bias potential module.

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
        };

        variable: union(enum) {
            potential_energy: Variable.PotentialEnergy,
        },

        function: union(enum) {
            harmonic: Function.Harmonic,
        },

        /// Evaluate the bias potential at given position.
        pub fn evaluate(self: @This(), variable: T) T {
            const value = switch (self.variable) {
                .potential_energy => self.variable.potential_energy.value,
            };

            return switch (self.function) {
                .harmonic => |field| 0.5 * field.k * (variable - value) * (variable - value),
            };
        }

        /// Evaluate the force from the bias potential for the specified coordinate index.
        pub fn evaluateForce(self: @This(), variable: T) T {
            const value = switch (self.variable) {
                .potential_energy => self.variable.potential_energy.value,
            };

            return switch (self.function) {
                .harmonic => |field| -field.k * (variable - value),
            };
        }

        /// Apply the bias potential to the adiabatic potential energy matrix.
        pub fn apply(self: @This(), adiabatic_potential: *RealMatrix(T)) void {
            for (0..adiabatic_potential.rows) |i| {

                const variable = switch (self.variable) {
                    .potential_energy => adiabatic_potential.at(i, i)
                };

                adiabatic_potential.ptr(i, i).* += self.evaluate(variable);
            }
        }

        /// Evaluate the force from the bias potential for the specified system state and coordinate index.
        pub fn force(self: @This(), adiabatic_potential: RealMatrix(T), state: usize, _: usize) T {
            const variable = switch (self.variable) {
                .potential_energy => adiabatic_potential.at(state, state)
            };

            return self.evaluateForce(variable);
         }
    };
}
