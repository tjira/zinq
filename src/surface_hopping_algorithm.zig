//! General union type for surface hopping algorithms.

const std = @import("std");

const classical_particle = @import("classical_particle.zig");
const complex_vector = @import("complex_vector.zig");
const fewest_switches = @import("fewest_switches.zig");
const landau_zener = @import("landau_zener.zig");
const mapping_approach = @import("mapping_approach.zig");
const object_array = @import("object_array.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const ClassicalParticle = classical_particle.ClassicalParticle;
const Complex = std.math.complex.Complex;
const ComplexVector = complex_vector.ComplexVector;
const FewestSwitches = fewest_switches.FewestSwitches;
const LandauZener = landau_zener.LandauZener;
const MappingApproach = mapping_approach.MappingApproach;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;
const RingBufferArray = object_array.RingBufferArray;

/// Parameters struct for all the surface hopping algorithms.
pub fn Parameters(comptime T: type) type {
    return struct {
        fs_parameters: fewest_switches.Parameters(T),
        lz_parameters: landau_zener.Parameters(T),
        ma_parameters: mapping_approach.Parameters(T)
    };
}

/// Electronic potential mode union.
pub fn SurfaceHoppingAlgorithm(comptime T: type) type {
    return union(enum) {
        fewest_switches: FewestSwitches(T),
        landau_zener: LandauZener(T),
        mapping_approach: MappingApproach(T),

        /// Get the jump probabilities for the current state.
        pub fn getJumpProbabilities(self: @This(), jump_probabilities: *RealVector(T), parameters: Parameters(T), current_state: usize) !void {
            switch (self) {
                .fewest_switches => |field| field.getJumpProbabilities(jump_probabilities, parameters.fs_parameters, current_state),
                .landau_zener => |field| field.getJumpProbabilities(jump_probabilities, parameters.lz_parameters, current_state),
                .mapping_approach => |field| try field.getJumpProbabilities(jump_probabilities, parameters.ma_parameters, current_state)
            }
        }

        /// Perform the jump based on the probabilities and adjust the momentum accordingly.
        pub fn jump(self: @This(), system: *ClassicalParticle(T), jump_probabilities: *RealVector(T), parameters: Parameters(T), adiabatic_potential: RealMatrix(T), state: usize, random: *std.Random) !usize {
            const substeps = if (self == .fewest_switches) self.fewest_switches.substeps else 1;

            var new_state: usize = state;

            for (0..substeps) |_| {

                try getJumpProbabilities(self, jump_probabilities, parameters, new_state);

                const random_number = random.float(T);

                if (jump_probabilities.sum() == 0) continue;

                var cumulative_probability: T = 0;

                if (new_state == state) for (0..jump_probabilities.len) |i| {

                    if (i == state) continue;

                    cumulative_probability += jump_probabilities.at(i);

                    if (random_number < cumulative_probability) {
                        new_state = @intCast(i); break;
                    }
                };
            }

            if (new_state != state) {

                const kinetic_energy = system.kineticEnergy();
                const current_energy = adiabatic_potential.at(state, state);
                const new_energy = adiabatic_potential.at(new_state, new_state);
                const energy_gap = new_energy - current_energy;

                if (kinetic_energy < energy_gap) return state;

                const factor = std.math.sqrt((kinetic_energy - energy_gap) / kinetic_energy);

                system.velocity.muls(factor);
            }

            return new_state;
        }
    };
}

/// Function to apply decoherence correction to the coefficients after a jump.
pub fn applyDecoherenceCorrection(comptime T: type, coefficients: *ComplexVector(T), adiabatic_potential: RealMatrix(T), kinetic_energy: T, current_state: usize, time_step: T, alpha: T) void {
    for (0..coefficients.len) |j| if (j != current_state) {

        const tau = (1 + alpha / kinetic_energy) / @abs(adiabatic_potential.at(j, j) - adiabatic_potential.at(current_state, current_state));

        coefficients.ptr(j).* = coefficients.at(j).mul(Complex(T).init(std.math.exp(-0.5 * time_step / tau), 0));
    };

    var sumc: T = 0; for (0..coefficients.len) |j| if (j != current_state) {
        sumc += coefficients.at(j).magnitude() * coefficients.at(j).magnitude();
    };

    if (coefficients.at(current_state).magnitude() > 0) {

        const num = 1 - sumc; const denom = coefficients.at(current_state).magnitude() * coefficients.at(current_state).magnitude();

        coefficients.ptr(current_state).* = coefficients.at(current_state).mul(Complex(T).init(std.math.sqrt(num / denom), 0));
    }
}
