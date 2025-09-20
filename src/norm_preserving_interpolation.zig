//! Norm Preserving interpolation coupling implementation.

const std = @import("std");

const error_handling = @import("error_handling.zig");
const real_matrix = @import("real_matrix.zig");

const Complex = std.math.Complex;
const RealMatrix = real_matrix.RealMatrix;

const throw = error_handling.throw;

/// Norm Preserving interpolation coupling implementation.
pub fn NormPreservingInterpolation(comptime T: type) type {
    return struct {

        /// Evaluate the time derivative coupling.
        pub fn evaluate(_: @This(), derivative_coupling: *RealMatrix(T), eigenvector_overlap: RealMatrix(T), time_step: T) !void {
            if (derivative_coupling.rows != 2) return throw(void, "NORM PRESERVING INTERPOLATION ONLY IMPLEMENTED FOR 2 STATES", .{});

            derivative_coupling.zero();

            if (eigenvector_overlap.at(0, 1) == 0) return;

            const a = Complex(T).init(eigenvector_overlap.at(0, 0), 0);
            const b = Complex(T).init(eigenvector_overlap.at(0, 1), 0);
            const c = Complex(T).init(eigenvector_overlap.at(1, 1), 0);

            const tau = std.math.complex.sqrt(a.sub(c).mul(a.sub(c)).sub(b.mul(b).mul(Complex(T).init(4, 0))));

            const log1 = std.math.complex.log(a.add(c).add(tau));
            const log2 = std.math.complex.log(a.add(c).sub(tau));

            derivative_coupling.ptr(0, 1).* = b.mul(log1.sub(log2)).div(tau).re / time_step;

            derivative_coupling.ptr(1, 0).* = -derivative_coupling.at(0, 1);
        }
    };
}
