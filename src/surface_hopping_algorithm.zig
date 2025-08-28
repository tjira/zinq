//! General union type for surface hopping algorithms.

const std = @import("std");

const landau_zener = @import("landau_zener.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");
const object_array = @import("object_array.zig");
const classical_particle = @import("classical_particle.zig");

const LandauZener = landau_zener.LandauZener;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;
const ClassicalParticle = classical_particle.ClassicalParticle;
const RingBufferArray = object_array.RingBufferArray;

/// Parameters struct for all the surface hopping algorithms.
pub fn Parameters(comptime T: type) type {
    return struct {
        adiabatic_potential: RealMatrix(T),
        lz_parameters: landau_zener.Parameters(T)
    };
}

/// Electronic potential mode union.
pub fn SurfaceHoppingAlgorithm(comptime T: type) type {
    return union(enum) {
        landau_zener: LandauZener(T),

        /// Get the jump probabilities for the current state.
        pub fn getJumpProbabilities(self: @This(), jump_probabilities: *RealVector(T), parameters: Parameters(T), current_state: u32) void {
            switch (self) {
                .landau_zener => |field| field.getJumpProbabilities(jump_probabilities, parameters.lz_parameters, current_state)
            }
        }

        /// Perform the jump based on the probabilities and adjust the momentum accordingly.
        pub fn jump(self: @This(), system: *ClassicalParticle(T), jump_probabilities: *RealVector(T), parameters: Parameters(T), current_state: u32, random: *std.Random) u32 {
            getJumpProbabilities(self, jump_probabilities, parameters, current_state);

            if (jump_probabilities.sum() == 0) return current_state;

            var cumulative_probability: T = 0;
            var new_state: u32 = current_state;
            const random_number = random.float(T);

            for (0..jump_probabilities.len) |i| {

                if (i == current_state) continue;

                cumulative_probability += jump_probabilities.at(i);

                if (random_number < cumulative_probability) new_state = @intCast(i);
            }

            if (new_state != current_state) {

                const kinetic_energy = system.kineticEnergy();
                const current_energy = parameters.adiabatic_potential.at(current_state, current_state);
                const new_energy = parameters.adiabatic_potential.at(new_state, new_state);
                const energy_gap = new_energy - current_energy;

                if (kinetic_energy < energy_gap) return current_state;

                const factor = std.math.sqrt((kinetic_energy - energy_gap) / kinetic_energy);

                system.velocity.muls(factor);
            }

            return new_state;
        }
    };
}
