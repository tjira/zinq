//! Norm Preserving interpolation coupling implementation.

const std = @import("std");

const real_matrix = @import("real_matrix.zig");

const Complex = std.math.Complex;
const RealMatrix = real_matrix.RealMatrix;

/// Norm Preserving interpolation coupling implementation.
pub fn HammesSchifferTully(comptime T: type) type {
    return struct {

        /// Evaluate the time derivative coupling.
        pub fn evaluate(self: @This(), derivative_coupling: *RealMatrix(T), eigenvector_overlap: RealMatrix(T), time_step: T) !void {
            _ = self; derivative_coupling.zero();

            for (0..derivative_coupling.rows) |i| for (i + 1..derivative_coupling.cols) |j| {

                const overlap_difference = eigenvector_overlap.at(i, j) - eigenvector_overlap.at(j, i);

                derivative_coupling.ptr(i, j).* = overlap_difference / (2 * time_step);

                derivative_coupling.ptr(j, i).* = -derivative_coupling.at(i, j);
            };
        }
    };
}
