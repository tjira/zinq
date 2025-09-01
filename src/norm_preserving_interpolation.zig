//! Norm Preserving interpolation coupling implementation.

const std = @import("std");

const real_matrix = @import("real_matrix.zig");

const Complex = std.math.Complex;
const RealMatrix = real_matrix.RealMatrix;

/// Norm Preserving interpolation coupling implementation.
pub fn NormPreservingInterpolation(comptime T: type) type {
    return struct {

        /// Evaluate the time derivative coupling.
        pub fn evaluate(self: @This(), TDC: *RealMatrix(T), S: RealMatrix(T), time_step: T) !void {
            if (TDC.rows != 2) return error.NotImplemented;

            const a = Complex(T).init(S.at(0, 0), 0);
            const b = Complex(T).init(S.at(0, 1), 0);
            const c = Complex(T).init(S.at(1, 1), 0);

            _ = self; TDC.zero();

            const tau = std.math.complex.sqrt(a.sub(c).mul(a.sub(c)).sub(b.mul(b).mul(Complex(T).init(4, 0))));

            const log1 = std.math.complex.log(a.add(c).add(tau));
            const log2 = std.math.complex.log(a.add(c).sub(tau));

            TDC.ptr(0, 1).* = b.mul(log1.sub(log2)).div(tau).re / time_step;

            TDC.ptr(1, 0).* = -TDC.at(0, 1);
        }
    };
}
