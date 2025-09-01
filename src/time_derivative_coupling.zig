//! General union type for different time derivative coupling implementations.

const std = @import("std");

const norm_preserving_interpolation = @import("norm_preserving_interpolation.zig");
const real_matrix = @import("real_matrix.zig");

const NormPreservingInterpolation = norm_preserving_interpolation.NormPreservingInterpolation;
const RealMatrix = real_matrix.RealMatrix;

/// Electronic potential mode union.
pub fn TimeDerivativeCoupling(comptime T: type) type {
    return union(enum) {
        npi: NormPreservingInterpolation(T),

        /// Evaluate the time derivative coupling.
        pub fn evaluate(self: @This(), TDC: *RealMatrix(T), S: RealMatrix(T), time_step: T) !void {
            switch (self) {
                inline else => |field| try field.evaluate(TDC, S, time_step)
            }
        }
    };
}
