//! File that contains the Fewest Switches Surface Hopping algorithm.

const std = @import("std");

const complex_runge_kutta = @import("complex_runge_kutta.zig");
const complex_vector = @import("complex_vector.zig");
const global_variables = @import("global_variables.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const Complex = std.math.complex.Complex;
const ComplexRungeKutta = complex_runge_kutta.ComplexRungeKutta;
const ComplexVector = complex_vector.ComplexVector;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const FSSH_DENOMINATOR_OFFSET = global_variables.FSSH_DENOMINATOR_OFFSET;

/// Parameters for the Fewest Switches method.
pub fn Parameters(comptime T: type) type {
    return struct {
        adiabatic_potential: RealMatrix(T),
        amplitudes: *ComplexVector(T),
        derivative_coupling: RealMatrix(T),
        runge_kutta: ComplexRungeKutta(T),
        time_step: T
    };
}

/// FSSH struct.
pub fn FewestSwitches(comptime T: type) type {
    return struct {
        substeps: u32 = 10,

        /// Get the jump probabilities for the current state.
        pub fn getJumpProbabilities(self: @This(), jump_probabilities: *RealVector(T), parameters: Parameters(T), current_state: usize) void {
            jump_probabilities.zero();

            const quantum_step = parameters.time_step / @as(T, @floatFromInt(self.substeps));

            const coefficient_derivative = struct {
                pub fn get(k: *ComplexVector(T), amplitudes: ComplexVector(T), derivative_parameters: anytype) void {
                    const adiabatic_potential = derivative_parameters.adiabatic_potential;
                    const derivative_coupling = derivative_parameters.derivative_coupling;

                    for (0..amplitudes.len) |i| {
                        k.ptr(i).* = amplitudes.at(i).mul(Complex(T).init(adiabatic_potential.at(i, i), 0)).mulbyi().neg();
                    }

                    for (0..amplitudes.len) |i| for (0..amplitudes.len) |j| {
                        k.ptr(i).* = k.at(i).sub(amplitudes.at(j).mul(Complex(T).init(derivative_coupling.at(i, j), 0)));
                    };
                }
            };

            const coefficient_derivative_parameters = .{
                .adiabatic_potential = parameters.adiabatic_potential,
                .derivative_coupling = parameters.derivative_coupling
            };

            @constCast(&parameters.runge_kutta).propagate(parameters.amplitudes, coefficient_derivative.get, coefficient_derivative_parameters, quantum_step);

            for (0..parameters.amplitudes.len) |j| if (j != current_state) {

                const re = parameters.amplitudes.at(j).mul(parameters.amplitudes.at(current_state).conjugate()).re;

                const denominator = std.math.pow(T, parameters.amplitudes.at(current_state).magnitude(), 2) + FSSH_DENOMINATOR_OFFSET;

                const p = 2 * parameters.derivative_coupling.at(current_state, j) * re / denominator * quantum_step;

                jump_probabilities.ptr(j).* = @max(p, 0);
            };
        }
    };
}
