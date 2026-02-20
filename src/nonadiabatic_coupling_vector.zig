//! File to calculate time derivative coupling using the NACV method.

const std = @import("std");

const device_read = @import("device_read.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const electronic_potential = @import("electronic_potential.zig");
const error_handling = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const fixGauge = eigenproblem_solver.fixGauge;
const readRealMatrix = device_read.readRealMatrix;
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
        previous_nacv: RealVector(T),
        position: RealVector(T),
        velocity: RealVector(T),
        time: T,
        dir: std.fs.Dir,
        allocator: std.mem.Allocator
    };
}

/// Nonadiabatic Coupling Vector (NACV) method for calculating time derivative couplings.
pub fn NonadiabaticCouplingVector(comptime T: type) type {
    return struct {
        finite_differences_step: T = 1e-8,
        minimum_energy_gap: T = 0,

        /// Evaluate the time derivative coupling.
        pub fn evaluate(self: @This(), derivative_coupling: *RealMatrix(T), parameters: Parameters(T)) !void {
            if (parameters.electronic_potential == .ab_initio) return try self.evaluateAbInitio(derivative_coupling, parameters);

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

                for (0..derivative_coupling.rows) |j| for (j + 1..derivative_coupling.cols) |k| {

                    if (@abs(adiabatic_potential.at(j, j) - adiabatic_potential.at(k, k)) > self.minimum_energy_gap) continue;

                    for (0..derivative_coupling.cols) |l| {

                        const bra = adiabatic_eigenvectors.at(l, j); const ket = (eigenvectors_plus.at(l, k) - eigenvectors_minus.at(l, k)) / (2 * self.finite_differences_step);

                        derivative_coupling.ptr(j, k).* += bra * ket * velocity.at(i);
                        derivative_coupling.ptr(k, j).* -= bra * ket * velocity.at(i);
                    }
                };
            }
        }

        /// Evaluate nonadiabatic coupling vector for an ab initio potential.
        pub fn evaluateAbInitio(self: @This(), derivative_coupling: *RealMatrix(T), params: Parameters(T)) !void {
            const dirname = try params.dir.realpathAlloc(params.allocator, "."); defer params.allocator.free(dirname);
            const path = try std.mem.concat(params.allocator, u8, &.{dirname, "/NACV.mat"}); defer params.allocator.free(path);

            var NACV = try readRealMatrix(T, path, params.allocator); defer NACV.deinit(params.allocator);

            var l: usize = 0; derivative_coupling.zero();

            if (params.time > 0) for (0..NACV.rows / params.velocity.len) |i| {

                var overlap: T = 0;

                for (0..params.velocity.len) |j| overlap += params.previous_nacv.at(i * params.velocity.len + j) * NACV.at(i * params.velocity.len + j, 0);

                if (overlap < 0) for (0..params.velocity.len) |j| {NACV.ptr(i * params.velocity.len + j, 0).* *= -1;};
            };

            l = 0;

            for (0..derivative_coupling.rows) |i| for (i + 1..derivative_coupling.cols) |j| {

                if (@abs(params.adiabatic_potential.at(i, i) - params.adiabatic_potential.at(j, j)) > self.minimum_energy_gap) continue;

                for (0..params.velocity.len) |k| {
                    derivative_coupling.ptr(i, j).* += params.velocity.at(k) * NACV.at(l, 0); l += 1;
                }

                derivative_coupling.ptr(j, i).* = -derivative_coupling.at(i, j);
            };

            for (0..NACV.rows) |i| @constCast(&params.previous_nacv).ptr(i).* = NACV.at(i, 0);
        }
    };
}
