//! Struct holding parameters for a complex absorbing potential (CAP).

const std = @import("std");

const real_vector = @import("real_vector.zig");

const RealVector = real_vector.RealVector;

/// Struct holding parameters for a complex absorbing potential (CAP).
pub fn ComplexAbsorbingPotential(comptime T: type) type {
    return struct {
        limits: []const []const T,
        strength: T = 0.001,

        /// Get the value of the CAP at the specified point.
        pub fn apply(self: @This(), position: RealVector(T)) T {
            var cap: T = 0;

            for (0..position.len) |i| {
                if (position.at(i) < self.limits[i][0]) cap -= self.strength * (std.math.exp(self.limits[i][0] - position.at(i)) - 1);
                if (position.at(i) > self.limits[i][1]) cap -= self.strength * (std.math.exp(position.at(i) - self.limits[i][1]) - 1);
            }

            return cap;
        }
    };
}
