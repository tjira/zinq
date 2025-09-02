//! General union type for different time derivative coupling implementations.

const std = @import("std");

const hammes_schiffer_tully = @import("hammes_schiffer_tully.zig");
const norm_preserving_interpolation = @import("norm_preserving_interpolation.zig");
const real_matrix = @import("real_matrix.zig");

const HammesSchifferTully = hammes_schiffer_tully.HammesSchifferTully;
const NormPreservingInterpolation = norm_preserving_interpolation.NormPreservingInterpolation;
const RealMatrix = real_matrix.RealMatrix;

/// Time derivative coupling union.
pub fn DerivativeCoupling(comptime T: type) type {
    return union(enum) {
        hst: HammesSchifferTully(T),
        npi: NormPreservingInterpolation(T),

        /// Evaluate the time derivative coupling.
        pub fn evaluate(self: @This(), derivative_coupling: *RealMatrix(T), S: RealMatrix(T), time_step: T) !void {
            switch (self) {
                inline else => |field| try field.evaluate(derivative_coupling, S, time_step)
            }
        }
    };
}
