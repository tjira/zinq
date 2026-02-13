//! File that contains the Berendsens thermostat implementation.

const std = @import("std");

const real_vector = @import("real_vector.zig");

const RealVector = real_vector.RealVector;

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
            const lambda = std.math.sqrt(1 + parameters.time_step / self.tau * (self.temperature / parameters.temperature - 1));

            for (0..velocities.len) |i| velocities.ptr(i).* *= lambda;
        }
    };
}
