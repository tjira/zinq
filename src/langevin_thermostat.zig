//! File that contains the Langevin thermostat implementation.

const std = @import("std");

const global_variables = @import("global_variables.zig");
const real_vector = @import("real_vector.zig");

const RealVector = real_vector.RealVector;

const AU2K = global_variables.AU2K;

/// Parameters for the Langevin thermostat.
pub fn Parameters(comptime T: type) type {
    return struct {
        time_step: T,
        masses: RealVector(T),
        random: *std.Random,
    };
}

/// Langevin thermostat implementation.
pub fn Langevin(comptime T: type) type {
    return struct {
        tau: T = 0.1,
        temperature: T = 298.15,

        /// Apply the Langevin thermostat to the velocities. Should be applied twice per time step, once before and once after the velocity Verlet integration.
        pub fn apply(self: @This(), velocities: *RealVector(T), parameters: Parameters(T)) !void {
            const c1 = std.math.exp(-parameters.time_step / self.tau / 2); const c2 = self.temperature * (1 - c1 * c1) / AU2K;

            for (0..velocities.len) |i| {
                velocities.ptr(i).* = velocities.at(i) * c1 + std.math.sqrt(c2 / parameters.masses.at(i)) * parameters.random.floatNorm(T);
            }
        }
    };
}
