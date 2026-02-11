//! Struct holding the Landau-Zener surface hopping algorithm implementation.

const std = @import("std");

const real_vector = @import("real_vector.zig");
const object_array = @import("object_array.zig");

const RealVector = real_vector.RealVector;
const RingBufferArray = object_array.RingBufferArray;

/// Parameters for the Landau-Zener method.
pub fn Parameters(comptime T: type) type {
    return struct {
        energy_gaps: RingBufferArray(T),
        time_step: T
    };
}

/// Struct holding the Landau-Zener surface hopping algorithm implementation.
pub fn LandauZener(comptime T: type) type {
    return struct {
        pub const ThreeStateVariant = enum {
            maximum_curvature,
            maximum_probability,
            nearest,
            normalied_probability
        };

        three_state_variant: ThreeStateVariant = .maximum_curvature,

        /// Get the jump probabilities for the current state.
        pub fn getJumpProbabilities(self: @This(), jump_probabilities: *RealVector(T), parameters: Parameters(T), current_state: usize) void {
            jump_probabilities.zero(); var maxddZ1: T = -std.math.inf(T); var minZ0: T = std.math.inf(T);

            const nstate = (1 + std.math.sqrt(1 + 8 * parameters.energy_gaps.len)) / 2;

            for (0..jump_probabilities.len) |i| if (i != current_state) {

                const coupling_index = @min(current_state, i) * (2 * nstate - @min(current_state, i) - 1) / 2 + (@max(current_state, i) - @min(current_state, i) - 1);

                const Z0 = parameters.energy_gaps.at(coupling_index).last(0);
                const Z1 = parameters.energy_gaps.at(coupling_index).last(1);
                const Z2 = parameters.energy_gaps.at(coupling_index).last(2);

                const dZ0 = (Z0 - Z1) / parameters.time_step;
                const dZ1 = (Z1 - Z2) / parameters.time_step;

                const ddZ1 = (Z0 - 2 * Z1 + Z2) / parameters.time_step / parameters.time_step;

                if (dZ0 * dZ1 > 0 or (dZ0 * dZ1 < 0 and ddZ1 < 0)) continue;

                const A = ddZ1 / 2; const B = (Z2 - Z0) / (2 * parameters.time_step); const C = Z1;

                const t0 = -B / (2 * A); const g = A * t0 * t0 + B * t0 + C;

                const veff = std.math.sqrt(g * ddZ1); const delta: T = 0.25 * std.math.pow(T, g, 2) / veff;

                var p = std.math.exp(-2 * std.math.pi * delta); if (std.math.isNan(p)) p = 0;

                if (jump_probabilities.len > 2 and self.three_state_variant != .normalied_probability) {

                    if (self.three_state_variant == .maximum_curvature and ddZ1 > maxddZ1) {
                        maxddZ1 = ddZ1; jump_probabilities.zero(); jump_probabilities.ptr(i).* = p;
                    }

                    else if (self.three_state_variant == .maximum_probability and jump_probabilities.sum() < p) {
                        jump_probabilities.zero(); jump_probabilities.ptr(i).* = p;
                    }

                    else if (self.three_state_variant == .nearest and Z0 < minZ0) {
                        minZ0 = Z0; jump_probabilities.zero(); jump_probabilities.ptr(i).* = p;
                    }
                }

                else jump_probabilities.ptr(i).* = p;
            };

            const probability_sum = jump_probabilities.sum();

            if (probability_sum > 1) for (0..jump_probabilities.len) |i| {
                jump_probabilities.ptr(i).* /= probability_sum;
            };
        }
    };
}
