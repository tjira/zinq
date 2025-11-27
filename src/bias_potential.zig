//! Bias potential module.

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
    };
}
