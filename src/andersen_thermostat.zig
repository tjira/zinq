//! File that contains the Andersen thermostat implementation.

const std = @import("std");

const global_variables = @import("global_variables.zig");
const real_vector = @import("real_vector.zig");

const RealVector = real_vector.RealVector;

const AU2K = global_variables.AU2K;

/// Parameters for the Andersen thermostat.
pub fn Parameters(comptime T: type) type {
    return struct {
        time_step: T,
        masses: RealVector(T),
        random: *std.Random
    };
}

/// Andersen thermostat implementation.
pub fn Andersen(comptime T: type) type {
    return struct {
        tau: T = 0.1,
        temperature: T = 298.15,

        /// Apply the Andersen thermostat.
        pub fn apply(self: @This(), velocities: *RealVector(T), parameters: Parameters(T)) !void {
            const prob = 1 - std.math.exp(-parameters.time_step / self.tau);

            for (0..velocities.len) |i| {

                const rnorm = parameters.random.floatNorm(T);

                if (parameters.random.float(T) < prob) {
                    velocities.ptr(i).* = rnorm * std.math.sqrt(self.temperature / parameters.masses.at(i) / AU2K);
                }
            }
        }
    };
}
