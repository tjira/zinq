//! File that contains the Berendsens thermostat implementation.

const std = @import("std");

const global_variables = @import("global_variables.zig");
const real_vector = @import("real_vector.zig");

const RealVector = real_vector.RealVector;

const BERENDSEN_TOLERANCE = 1e-6;

/// Parameters for the Berendsen thermostat.
pub fn Parameters(comptime T: type) type {
    return struct {
        time_step: T,
        temperature: T
    };
}

/// Berendsen thermostat implementation.
pub fn Berendsen(comptime T: type) type {
    return struct {
        tau: T = 0.1,
        temperature: T = 298.15,

        /// Apply the Berendsen thermostat.
        pub fn apply(self: @This(), velocities: *RealVector(T), parameters: Parameters(T)) !void {
            if (parameters.temperature < 1e-6) return;

            const lambda = std.math.sqrt(1 + parameters.time_step / self.tau * (self.temperature / parameters.temperature - 1));

            for (0..velocities.len) |i| velocities.ptr(i).* *= lambda;
        }
    };
}
