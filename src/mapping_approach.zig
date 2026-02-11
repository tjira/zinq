//! File that contains the Mapping Approach to Surface Hopping algorithm.

const std = @import("std");

const complex_runge_kutta = @import("complex_runge_kutta.zig");
const complex_vector = @import("complex_vector.zig");
const error_handling = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const Complex = std.math.complex.Complex;
const ComplexRungeKutta = complex_runge_kutta.ComplexRungeKutta;
const ComplexVector = complex_vector.ComplexVector;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const mm = matrix_multiplication.mm;

const FSSH_DENOMINATOR_OFFSET = global_variables.FSSH_DENOMINATOR_OFFSET;

/// Parameters for the Mapping Approach method.
pub fn Parameters(comptime T: type) type {
    return struct {
        adiabatic_potential: RealMatrix(T),
        bloch_vector: *RealVector(T),
        derivative_coupling: RealMatrix(T),
        time_step: T
    };
}

/// MASH struct.
pub fn MappingApproach(comptime T: type) type {
    return struct {

        /// Get the jump probabilities for the current state.
        pub fn getJumpProbabilities(_: @This(), jump_probabilities: *RealVector(T), parameters: Parameters(T), current_state: usize) !void {
            if (jump_probabilities.len != 2) return error_handling.throw(void, "MAPPING APPROACH ONLY IMPLEMENTED FOR 2 STATES", .{});

            jump_probabilities.zero(); var omega_data: [9]T = undefined; var new_bloch_data: [3]T = undefined;

            const V = (parameters.adiabatic_potential.at(1, 1) - parameters.adiabatic_potential.at(0, 0)) * parameters.time_step;
            const C = parameters.derivative_coupling.at(1, 0) * parameters.time_step;

            const a = V * V + 4 * C * C; const b = Complex(T).init(0, std.math.sqrt(a));
            const sinhb = std.math.complex.sinh(b); const coshb = std.math.complex.cosh(b);

            omega_data[0] = coshb.re;
            omega_data[1] = -V * sinhb.div(b).re;
            omega_data[2] = 2 * C * sinhb.div(b).re;
            omega_data[3] = -omega_data[1];
            omega_data[4] = (4 * C * C + V * V * coshb.re) / a;
            omega_data[5] = -2 * C * V * (coshb.re - 1) / a;
            omega_data[6] = -omega_data[2];
            omega_data[7] = omega_data[5];
            omega_data[8] = (V * V + 4 * C * C * coshb.re) / a;

            const omega = RealMatrix(T){.rows = 3, .cols = 3, .data = &omega_data};
            var new_bloch = RealMatrix(T){.rows = 3, .cols = 1, .data = &new_bloch_data};

            try mm(T, &new_bloch, omega, false, parameters.bloch_vector.*.asMatrix(), false); try new_bloch.asVector().copyTo(parameters.bloch_vector);

            if (current_state == 1 and parameters.bloch_vector.at(2) < 0) jump_probabilities.ptr(0).* = 1;
            if (current_state == 0 and parameters.bloch_vector.at(2) > 0) jump_probabilities.ptr(1).* = 1;
        }
    };
}
