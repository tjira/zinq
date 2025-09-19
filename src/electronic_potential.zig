//! General union type for different electronic potential models.

const std = @import("std");

const eigenproblem_solver = @import("eigenproblem_solver.zig");
const file_potential = @import("file_potential.zig");
const global_variables = @import("global_variables.zig");
const harmonic_potential = @import("harmonic_potential.zig");
const morse_potential = @import("morse_potential.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");
const time_linear_potential = @import("time_linear_potential.zig");
const tully_potential = @import("tully_potential.zig");

const FilePotential = file_potential.FilePotential;
const HarmonicPotential = harmonic_potential.HarmonicPotential;
const MorsePotential = morse_potential.MorsePotential;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;
const TimeLinearPotential = time_linear_potential.TimeLinearPotential;
const TullyPotential1 = tully_potential.TullyPotential1;

const diagonalizeSymmetric = eigenproblem_solver.diagonalizeSymmetric;
const eigensystemSymmetric = eigenproblem_solver.eigensystemSymmetric;

const FINITE_DIFFERENCES_STEP = global_variables.FINITE_DIFFERENCES_STEP;

/// Electronic potential mode union.
pub fn ElectronicPotential(comptime T: type) type {
    return union(enum) {
        file: FilePotential(T),
        harmonic: HarmonicPotential(T),
        morse: MorsePotential(T),
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
                .file => |field| field.evaluateDiabatic(diabatic_potential, position, time) catch unreachable,
                inline else => |field| field.evaluateDiabatic(diabatic_potential, position, time)
            }
        }

        /// Evaluate adiabatic eigensystem at given system state and time.
        pub fn evaluateEigensystem(self: @This(), diabatic: *RealMatrix(T), adiabatic: *RealMatrix(T), eigenvectors: *RealMatrix(T), position: RealVector(T), time: T) !void {
            self.evaluateDiabatic(diabatic, position, time);

            try eigensystemSymmetric(T, adiabatic, eigenvectors, diabatic.*);
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
                .file => |field| field.ndim,
                .harmonic => |field| field.k.len,
                .morse => 1,
                .tully_1 => 1,
                .time_linear => 1
            };
        }

        /// Getter for number of electronic states.
        pub fn nstate(self: @This()) usize {
            return switch (self) {
                .file => |field| field.nstate,
                .harmonic => 1,
                .morse => 1,
                .tully_1 => 2,
                .time_linear => 2
            };
        }
    };
}
