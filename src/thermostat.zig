//! General union type for thermostats.

const std = @import("std");

const andersen_thermostat = @import("andersen_thermostat.zig");
const berendsen_thermostat = @import("berendsen_thermostat.zig");
const real_vector = @import("real_vector.zig");

const Andersen = andersen_thermostat.Andersen;
const Berendsen = berendsen_thermostat.Berendsen;
const RealVector = real_vector.RealVector;

/// Parameters struct for all the thermostats.
pub fn Parameters(comptime T: type) type {
    return struct {
        andersen_params: andersen_thermostat.Parameters(T),
        berendsen_params: berendsen_thermostat.Parameters(T)
    };
}

/// Electronic potential mode union.
pub fn Thermostat(comptime T: type) type {
    return union(enum) {
        andersen: Andersen(T),
        berendsen: Berendsen(T),

        /// Apply the thermostat.
        pub fn apply(self: @This(), velocities: *RealVector(T), parameters: Parameters(T)) !void {
            switch (self) {
                .andersen => |field| try field.apply(velocities, parameters.andersen_params),
                .berendsen => |field| try field.apply(velocities, parameters.berendsen_params)
            }
        }
    };
}
