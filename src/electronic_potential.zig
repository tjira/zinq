//! General union type for different electronic potential models.

const std = @import("std");

const device_write = @import("device_write.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const global_variables = @import("global_variables.zig");
const harmonic_potential = @import("harmonic_potential.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");
const time_linear_potential = @import("time_linear_potential.zig");
const tully_potential = @import("tully_potential.zig");

const HarmonicPotential = harmonic_potential.HarmonicPotential;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;
const TimeLinearPotential = time_linear_potential.TimeLinearPotential;
const TullyPotential1 = tully_potential.TullyPotential1;

const diagonalizeSymmetric = eigenproblem_solver.diagonalizeSymmetric;
const eigensystemSymmetric = eigenproblem_solver.eigensystemSymmetric;
const printRealMatrix = device_write.printRealMatrix;

const FINITE_DIFFERENCES_STEP = global_variables.FINITE_DIFFERENCES_STEP;

/// Electronic potential mode union.
pub fn ElectronicPotential(comptime T: type) type {
    return union(enum) {
        harmonic: HarmonicPotential(T),
        time_linear: TimeLinearPotential(T),
        tully_1: TullyPotential1(T),

        /// Evaluate the adabatic potential energy matrix at given system state and time.
        pub fn evaluateAdiabatic(self: @This(), adiabatic_potential: *RealMatrix(T), position: RealVector(T), time: T) !void {
            self.evaluateDiabatic(adiabatic_potential, position, time);

            try diagonalizeSymmetric(T, adiabatic_potential);
        }

        /// Evaluate the dabatic potential energy matrix at given system state and time.
        pub fn evaluateDiabatic(self: @This(), diabatic_potential: *RealMatrix(T), position: RealVector(T), time: T) void {
            switch (self) {
                inline else => |field| field.evaluateDiabatic(diabatic_potential, position, time)
            }
        }

        /// Evaluate adiabatic eigensystem at given system state and time.
        pub fn evaluateEigensystem(self: @This(), potential_eigensystem: anytype, position: RealVector(T), time: T) !void {
            self.evaluateDiabatic(&potential_eigensystem.diabatic_potential, position, time);

            try eigensystemSymmetric(T, &potential_eigensystem.adiabatic_potential, &potential_eigensystem.adiabatic_eigenvectors, potential_eigensystem.diabatic_potential);
        }

        /// Evaluate the adabatic potential energy gradient for a specific coordinate index. The adiabatic potential matrix will hold the potential at the r + dr point.
        pub fn forceAdiabatic(self: @This(), adiabatic_potential: *RealMatrix(T), position: RealVector(T), time: T, state: usize, index: usize) !T {
            const original_position = position.at(index);

            @constCast(&position).ptr(index).* = original_position - FINITE_DIFFERENCES_STEP;

            try self.evaluateAdiabatic(adiabatic_potential, position, time);

            const energy_minus = adiabatic_potential.at(state, state);

            @constCast(&position).ptr(index).* = original_position + FINITE_DIFFERENCES_STEP;

            try self.evaluateAdiabatic(adiabatic_potential, position, time);

            const energy_plus = adiabatic_potential.at(state, state);

            @constCast(&position).ptr(index).* = original_position;

            return 0.5 * (energy_minus - energy_plus) / FINITE_DIFFERENCES_STEP;
        }

        /// Getter for number of dimensions.
        pub fn ndim(self: @This()) usize {
            return switch (self) {
                .harmonic => |field| field.k.len,
                .tully_1 => 1,
                .time_linear => 1
            };
        }

        /// Getter for number of electronic states.
        pub fn nstate(self: @This()) usize {
            return switch (self) {
                .harmonic => |_| 1,
                .tully_1 => |_| 2,
                .time_linear => |_| 2
            };
        }
    };
}
