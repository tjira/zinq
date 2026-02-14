//! File that contains the Nose-Hoover thermostat implementation.

const std = @import("std");

const global_variables = @import("global_variables.zig");
const real_vector = @import("real_vector.zig");

const RealVector = real_vector.RealVector;

const AU2K = global_variables.AU2K;

/// Parameters for the Nose-Hoover thermostat.
pub fn Parameters(comptime T: type) type {
    return struct {
        time_step: T,
        ndof: T,
        masses: RealVector(T),
        xi: *T
    };
}

/// Nose-Hoover thermostat implementation.
pub fn NoseHoover(comptime T: type) type {
    return struct {
        tau: T = 0.1,
        temperature: T = 298.15,

        /// Apply the first stage of the Nose-Hoover thermostat (before velocity Verlet).
        pub fn applyBefore(self: @This(), velocities: *RealVector(T), parameters: Parameters(T)) !void {
            const dt = parameters.time_step; const ndof = parameters.ndof; var ekin: T = 0;

            for (0..velocities.len) |i| ekin += 0.5 * parameters.masses.at(i) * velocities.at(i) * velocities.at(i);

            const temp = self.temperature / AU2K; const Q = parameters.ndof * temp * self.tau * self.tau;

            parameters.xi.* += (dt / 2) * (2 * ekin - ndof * temp) / Q;

            for (0..velocities.len) |d| velocities.ptr(d).* *= std.math.exp(-parameters.xi.* * dt / 2);
        }

        /// Apply the second stage of the Nose-Hoover thermostat (after velocity Verlet).
        pub fn applyAfter(self: @This(), velocities: *RealVector(T), parameters: Parameters(T)) !void {
            const dt = parameters.time_step; const ndof = parameters.ndof; var ekin: T = 0;

            for (0..velocities.len) |d| velocities.ptr(d).* *= std.math.exp(-parameters.xi.* * dt / 2);

            for (0..velocities.len) |i| ekin += 0.5 * parameters.masses.at(i) * velocities.at(i) * velocities.at(i);

            const temp = self.temperature / AU2K; const Q = parameters.ndof * temp * self.tau * self.tau;

            parameters.xi.* += (dt / 2) * (2 * ekin - ndof * temp) / Q;
        }
    };
}
