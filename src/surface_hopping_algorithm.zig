//! General union type for surface hopping algorithms.

const std = @import("std");

const classical_particle = @import("classical_particle.zig");
const fewest_switches = @import("fewest_switches.zig");
const landau_zener = @import("landau_zener.zig");
const object_array = @import("object_array.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const ClassicalParticle = classical_particle.ClassicalParticle;
const FewestSwitches = fewest_switches.FewestSwitches;
const LandauZener = landau_zener.LandauZener;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;
const RingBufferArray = object_array.RingBufferArray;

/// Parameters struct for all the surface hopping algorithms.
pub fn Parameters(comptime T: type) type {
    return struct {
        fs_parameters: fewest_switches.Parameters(T),
        lz_parameters: landau_zener.Parameters(T)
    };
}

/// Electronic potential mode union.
pub fn SurfaceHoppingAlgorithm(comptime T: type) type {
    return union(enum) {
        fewest_switches: FewestSwitches(T),
        landau_zener: LandauZener(T),

        /// Get the jump probabilities for the current state.
        pub fn getJumpProbabilities(self: @This(), jump_probabilities: *RealVector(T), parameters: Parameters(T), current_state: usize) void {
            switch (self) {
                .fewest_switches => |field| field.getJumpProbabilities(jump_probabilities, parameters.fs_parameters, current_state),
                .landau_zener => |field| field.getJumpProbabilities(jump_probabilities, parameters.lz_parameters, current_state)
            }
        }

        /// Perform the jump based on the probabilities and adjust the momentum accordingly.
        pub fn jump(self: @This(), system: *ClassicalParticle(T), jump_probabilities: *RealVector(T), parameters: Parameters(T), adiabatic_potential: RealMatrix(T), state: usize, random: *std.Random) usize {
            const substeps = if (self == .fewest_switches) self.fewest_switches.substeps else 1;

            var new_state: usize = state;

            for (0..substeps) |_| {

                getJumpProbabilities(self, jump_probabilities, parameters, state);

                const random_number = random.float(T);

                if (jump_probabilities.sum() == 0) continue;

                var cumulative_probability: T = 0;

                for (0..jump_probabilities.len) |i| {

                    if (i == state) continue;

                    cumulative_probability += jump_probabilities.at(i);

                    if (random_number < cumulative_probability) new_state = @intCast(i);
                }
            }

            if (new_state != state) {

                const kinetic_energy = system.kineticEnergy();
                const current_energy = adiabatic_potential.at(state, state);
                const new_energy = adiabatic_potential.at(new_state, new_state);
                const energy_gap = new_energy - current_energy;

                if (kinetic_energy < energy_gap) return state;

                const factor = std.math.sqrt((kinetic_energy - energy_gap) / kinetic_energy);

                system.velocity.muls(factor);
            }

            return new_state;
        }
    };
}
