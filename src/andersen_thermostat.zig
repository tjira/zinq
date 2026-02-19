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
        molecule: bool,
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

            for (0..if (parameters.molecule) velocities.len / 3 else velocities.len) |i| {

                const rnorm1 = parameters.random.floatNorm(T);
                const rnorm2 = parameters.random.floatNorm(T);
                const rnorm3 = parameters.random.floatNorm(T);

                if (parameters.random.float(T) < prob) {

                    if (parameters.molecule) {
                        velocities.ptr(i * 3 + 0).* = rnorm1 * std.math.sqrt(self.temperature / parameters.masses.at(3 * i + 0) / AU2K);
                        velocities.ptr(i * 3 + 1).* = rnorm2 * std.math.sqrt(self.temperature / parameters.masses.at(3 * i + 1) / AU2K);
                        velocities.ptr(i * 3 + 2).* = rnorm3 * std.math.sqrt(self.temperature / parameters.masses.at(3 * i + 2) / AU2K);
                    }

                    else velocities.ptr(i).* = rnorm1 * std.math.sqrt(self.temperature / parameters.masses.at(i) / AU2K);
                }
            }
        }
    };
}
