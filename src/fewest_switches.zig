//! File that contains the Fewest Switches Surface Hopping algorithm.

const std = @import("std");

const complex_vector = @import("complex_vector.zig");
const global_variables = @import("global_variables.zig");
const real_vector = @import("real_vector.zig");
const real_matrix = @import("real_matrix.zig");

const ComplexVector = complex_vector.ComplexVector;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const FSSH_DENOMINATOR_OFFSET = global_variables.FSSH_DENOMINATOR_OFFSET;

/// Parameters for the Fewest Switches method.
pub fn Parameters(comptime T: type) type {
    return struct {
        derivative_coupling: RealMatrix(T),
        time_step: T,
        amplitudes: ComplexVector(T)
    };
}

/// FSSH struct.
pub fn FewestSwitches(comptime T: type) type {
    return struct {
        substeps: u32 = 100,

        /// Get the jump probabilities for the current state.
        pub fn getJumpProbabilities(self: @This(), jump_probabilities: *RealVector(T), parameters: Parameters(T), current_state: usize) void {
            _ = self; jump_probabilities.zero();

            for (0..parameters.amplitudes.len) |i| if (i != current_state) {

                const re = parameters.amplitudes.at(i).mul(parameters.amplitudes.at(current_state).conjugate()).re;

                const denominator = std.math.pow(T, parameters.amplitudes.at(current_state).magnitude(), 2) + FSSH_DENOMINATOR_OFFSET;

                jump_probabilities.ptr(i).* = 2 * parameters.derivative_coupling.at(i, current_state) * re / denominator * parameters.time_step;
            };
        }
    };
}
