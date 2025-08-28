//! General union type for surface hopping algorithms.

const std = @import("std");

const landau_zener = @import("landau_zener.zig");
const real_vector = @import("real_vector.zig");
const object_array = @import("object_array.zig");
const classical_particle = @import("classical_particle.zig");

const LandauZener = landau_zener.LandauZener;
const RealVector = real_vector.RealVector;
const ClassicalParticle = classical_particle.ClassicalParticle;
const RingBufferArray = object_array.RingBufferArray;

/// Electronic potential mode union.
pub fn SurfaceHoppingAlgorithm(comptime T: type) type {
    return union(enum) {
        landau_zener: LandauZener(T),

        /// Get the jump probabilities for the current state.
        pub fn getJumpProbabilities(self: @This(), jump_probabilities: *RealVector(T), energy_gaps: RingBufferArray(T), current_state: u32, time_step: T) void {
            switch (self) {
                inline else => |field| field.getJumpProbabilities(jump_probabilities, energy_gaps, current_state, time_step)
            }
        }

        /// Perform the jump based on the probabilities and adjust the momentum accordingly.
        pub fn jump(self: @This(), system: *ClassicalParticle(T), jump_probabilities: *RealVector(T), energy_gaps: RingBufferArray(T), current_state: u32, time_step: T, random: *std.Random) u32 {
            getJumpProbabilities(self, jump_probabilities, energy_gaps, current_state, time_step);

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

                const coupling_index = current_state + new_state - 1;

                const kinetic_energy = system.kineticEnergy();
                const energy_gap = energy_gaps.at(coupling_index).last(0);

                if (kinetic_energy < energy_gap) return current_state;

                const sign: T = if (new_state > current_state) 1 else -1;

                const factor = std.math.sqrt((kinetic_energy - sign * energy_gap) / kinetic_energy);

                system.velocity.muls(factor);
            }

            return new_state;
        }
    };
}
