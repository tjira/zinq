//! General union type for different electronic potential models.

const std = @import("std");

const bias_potential = @import("bias_potential.zig");
const custom_potential = @import("custom_potential.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const file_potential = @import("file_potential.zig");
const harmonic_potential = @import("harmonic_potential.zig");
const jahn_teller_potential = @import("jahn_teller_potential.zig");
const morse_potential = @import("morse_potential.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");
const time_linear_potential = @import("time_linear_potential.zig");
const tully_potential = @import("tully_potential.zig");
const vibronic_coupling_potential = @import("vibronic_coupling_potential.zig");

const BiasPotential = bias_potential.BiasPotential;
const CustomPotential = custom_potential.CustomPotential;
const FilePotential = file_potential.FilePotential;
const HarmonicPotential = harmonic_potential.HarmonicPotential;
const JahnTellerPotential = jahn_teller_potential.JahnTellerPotential;
const MorsePotential = morse_potential.MorsePotential;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;
const TimeLinearPotential = time_linear_potential.TimeLinearPotential;
const TullyPotential1 = tully_potential.TullyPotential1;
const VibronicCouplingPotential = vibronic_coupling_potential.VibronicCouplingPotential;

const diagonalizeHermitian = eigenproblem_solver.diagonalizeHermitian;
const eigensystemHermitian = eigenproblem_solver.eigensystemHermitian;

/// Electronic potential mode union.
pub fn ElectronicPotential(comptime T: type) type {
    return union(enum) {
        custom: CustomPotential(T),
        file: FilePotential(T),
        harmonic: HarmonicPotential(T),
        jahn_teller: JahnTellerPotential(T),
        morse: MorsePotential(T),
        time_linear: TimeLinearPotential(T),
        tully_1: TullyPotential1(T),
        vibronic_coupling: VibronicCouplingPotential(T),

        /// Evaluate the adabatic potential energy matrix at given system state and time.
        pub fn evaluateAdiabatic(self: @This(), adiabatic_potential: *RealMatrix(T), position: RealVector(T), time: T, bias: ?BiasPotential(T)) !void {
            try self.evaluateDiabatic(adiabatic_potential, position, time);

            if (bias) |bias_struct| for (0..adiabatic_potential.rows) |i| {

                const variable = switch (bias_struct.variable) {
                    .potential_energy => adiabatic_potential.at(i, i)
                };

                adiabatic_potential.ptr(i, i).* += bias_struct.evaluate(variable);
            };

            try diagonalizeHermitian(T, adiabatic_potential);
        }

        /// Evaluate the dabatic potential energy matrix at given system state and time.
        pub fn evaluateDiabatic(self: @This(), diabatic_potential: *RealMatrix(T), position: RealVector(T), time: T) !void {
            switch (self) {
                inline else => |field| try field.evaluateDiabatic(diabatic_potential, position, time)
            }
        }

        /// Evaluate the dabatic potential energy matrix element.
        pub fn evaluateDiabaticElement(self: @This(), i: usize, j: usize, position: RealVector(T), time: T) !T {
            return switch (self) {
                inline else => |field| field.evaluateDiabaticElement(i, j, position, time)
            };
        }

        /// Evaluate the matrix element of potential derivative.
        pub fn evaluateDiabaticElementDerivative1(self: @This(), i: usize, j: usize, position: RealVector(T), time: T, index: usize, fdiff_step: T) !T {
            const original_position = position.at(index);

            @constCast(&position).ptr(index).* = original_position - fdiff_step;

            const energy_minus = try self.evaluateDiabaticElement(i, j, position, time);

            @constCast(&position).ptr(index).* = original_position + fdiff_step;

            const energy_plus = try self.evaluateDiabaticElement(i, j, position, time);

            @constCast(&position).ptr(index).* = original_position;

            return 0.5 * (energy_plus - energy_minus) / fdiff_step;
        }

        /// Evaluate the matrix element of potential second derivative.
        pub fn evaluateDiabaticElementDerivative2(self: @This(), i: usize, j: usize, position: RealVector(T), time: T, index: usize, fdiff_step: T) !T {
            const original_position = position.at(index);

            @constCast(&position).ptr(index).* = original_position - fdiff_step;

            const energy_minus = try self.evaluateDiabaticElement(i, j, position, time);

            @constCast(&position).ptr(index).* = original_position;

            const energy = try self.evaluateDiabaticElement(i, j, position, time);

            @constCast(&position).ptr(index).* = original_position + fdiff_step;

            const energy_plus = try self.evaluateDiabaticElement(i, j, position, time);

            @constCast(&position).ptr(index).* = original_position;

            return (energy_plus - 2 * energy + energy_minus) / (fdiff_step * fdiff_step);
        }

        /// Evaluate adiabatic eigensystem at given system state and time.
        pub fn evaluateEigensystem(self: @This(), diabatic: *RealMatrix(T), adiabatic: *RealMatrix(T), eigenvectors: *RealMatrix(T), position: RealVector(T), time: T) !void {
            try self.evaluateDiabatic(diabatic, position, time);

            try eigensystemHermitian(T, adiabatic, eigenvectors, diabatic.*);
        }

        /// Evaluate the adabatic potential energy gradient for a specific coordinate index. The adiabatic potential matrix will hold the potential at the r + dr point.
        pub fn forceAdiabatic(self: @This(), adiabatic_potential: *RealMatrix(T), position: RealVector(T), time: T, state: usize, index: usize, fdiff_step: T, bias: ?BiasPotential(T)) !T {
            const original_position = position.at(index);

            @constCast(&position).ptr(index).* = original_position - fdiff_step;

            try self.evaluateAdiabatic(adiabatic_potential, position, time, bias);

            const energy_minus = adiabatic_potential.at(state, state);

            @constCast(&position).ptr(index).* = original_position + fdiff_step;

            try self.evaluateAdiabatic(adiabatic_potential, position, time, bias);

            const energy_plus = adiabatic_potential.at(state, state);

            @constCast(&position).ptr(index).* = original_position;

            return 0.5 * (energy_minus - energy_plus) / fdiff_step;
        }

        /// Getter for number of dimensions.
        pub fn ndim(self: @This()) usize {
            return switch (self) {
                .custom => |field| field.variables.len,
                .file => |field| field.ndim,
                .harmonic => |field| field.k.len,
                .jahn_teller => 2,
                .morse => 1,
                .time_linear => 1,
                .tully_1 => 1,
                .vibronic_coupling => |field| field.ndim()
            };
        }

        /// Getter for number of electronic states.
        pub fn nstate(self: @This()) usize {
            return switch (self) {
                .custom => |field| field.matrix.len,
                .file => |field| field.nstate,
                .harmonic => 1,
                .jahn_teller => 2,
                .morse => 1,
                .time_linear => 2,
                .tully_1 => 2,
                .vibronic_coupling => |field| field.nstate()
            };
        }
    };
}
