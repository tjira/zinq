//! File to calculate time derivative coupling using the NACV method.

const std = @import("std");

const eigenproblem_solver = @import("eigenproblem_solver.zig");
const electronic_potential = @import("electronic_potential.zig");
const error_handling = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const fixGauge = eigenproblem_solver.fixGauge;
const throw = error_handling.throw;

const ElectronicPotential = electronic_potential.ElectronicPotential;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const MAX_NACV_STATES = global_variables.MAX_NACV_STATES;

/// Parameters for the Nonadiabatic Coupling Vector method.
pub fn Parameters(comptime T: type) type {
    return struct {
        adiabatic_potential: RealMatrix(T),
        diabatic_potential: RealMatrix(T),
        adiabatic_eigenvectors: RealMatrix(T),
        electronic_potential: ElectronicPotential(T),
        position: RealVector(T),
        velocity: RealVector(T),
        time: T
    };
}

/// Nonadiabatic Coupling Vector (NACV) method for calculating time derivative couplings.
pub fn NonadiabaticCouplingVector(comptime T: type) type {
    return struct {
        finite_differences_step: T = 1e-8,

        /// Evaluate the time derivative coupling.
        pub fn evaluate(self: @This(), derivative_coupling: *RealMatrix(T), parameters: Parameters(T)) !void {
            if (derivative_coupling.rows > MAX_NACV_STATES or derivative_coupling.cols > MAX_NACV_STATES) {
                return throw(void, "MAXIMUM NUMBER OF STATES FOR NACV METHOD IS {d}", .{MAX_NACV_STATES});
            }

            derivative_coupling.zero();

            const adiabatic_potential = parameters.adiabatic_potential;
            const diabatic_potential = parameters.diabatic_potential;
            const adiabatic_eigenvectors = parameters.adiabatic_eigenvectors;
            const position = parameters.position;
            const velocity = parameters.velocity;
            const time = parameters.time;

            var data_plus: [MAX_NACV_STATES * MAX_NACV_STATES]T = undefined; var data_minus: [MAX_NACV_STATES * MAX_NACV_STATES]T = undefined;

            var eigenvectors_plus = RealMatrix(T){.data = &data_plus, .rows = adiabatic_eigenvectors.rows, .cols = adiabatic_eigenvectors.cols};
            var eigenvectors_minus = RealMatrix(T){.data = &data_minus, .rows = adiabatic_eigenvectors.rows, .cols = adiabatic_eigenvectors.cols};

            for (0..position.len) |i| {

                const original_position = position.at(i);

                @constCast(&position).ptr(i).* = original_position + self.finite_differences_step; 

                try parameters.electronic_potential.evaluateEigensystem(@constCast(&diabatic_potential), @constCast(&adiabatic_potential), &eigenvectors_plus, position, time);

                @constCast(&position).ptr(i).* = original_position - self.finite_differences_step;

                try parameters.electronic_potential.evaluateEigensystem(@constCast(&diabatic_potential), @constCast(&adiabatic_potential), &eigenvectors_minus, position, time);

                @constCast(&position).ptr(i).* = original_position;

                try fixGauge(T, &eigenvectors_plus, adiabatic_eigenvectors);
                try fixGauge(T, &eigenvectors_minus, adiabatic_eigenvectors);

                for (0..derivative_coupling.rows) |j| for (j + 1..derivative_coupling.cols) |k| for (0..derivative_coupling.cols) |l| {

                    const bra = adiabatic_eigenvectors.at(l, j); const ket = (eigenvectors_plus.at(l, k) - eigenvectors_minus.at(l, k)) / (2 * self.finite_differences_step);

                    derivative_coupling.ptr(j, k).* += bra * ket * velocity.at(i);
                    derivative_coupling.ptr(k, j).* -= bra * ket * velocity.at(i);
                };
            }
        }
    };
}
