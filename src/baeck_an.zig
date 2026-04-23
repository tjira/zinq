//! Baeck and An's method for computing time derivative couplings in nonadiabatic molecular dynamics simulations.

const std = @import("std");

const object_array = @import("object_array.zig");
const real_matrix = @import("real_matrix.zig");

const Complex = std.math.Complex;
const RealMatrix = real_matrix.RealMatrix;
const RingBufferArray = object_array.RingBufferArray;

/// Parameters for the Baeck and An method.
pub fn Parameters(comptime T: type) type {
    return struct { energy_gaps: RingBufferArray(T), time_step: T };
}

/// Baeck and An's method for computing time derivative couplings.
pub fn BaeckAn(comptime T: type) type {
    return struct {
        /// Evaluate the time derivative coupling.
        pub fn evaluate(_: @This(), derivative_coupling: *RealMatrix(T), parameters: Parameters(T)) !void {
            derivative_coupling.zero();

            const time_step = parameters.time_step;

            for (0..derivative_coupling.rows) |i| for (i + 1..derivative_coupling.cols) |j| {
                const coupling_index = i * (2 * derivative_coupling.cols - i - 1) / 2 + (j - i - 1);

                const Z0 = parameters.energy_gaps.at(coupling_index).last(0);
                const Z1 = parameters.energy_gaps.at(coupling_index).last(1);
                const Z2 = parameters.energy_gaps.at(coupling_index).last(2);

                const ddZ0 = (Z0 - 2 * Z1 + Z2) / time_step / time_step;

                const sigma = if (ddZ0 > 1e-14) 0.5 * std.math.sqrt(ddZ0 / Z0) else 0;

                derivative_coupling.ptr(i, j).* = sigma;
                derivative_coupling.ptr(j, i).* = -derivative_coupling.at(i, j);
            };
        }
    };
}
