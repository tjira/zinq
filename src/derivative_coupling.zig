//! General union type for different time derivative coupling implementations.

const std = @import("std");

const hammes_schiffer_tully = @import("hammes_schiffer_tully.zig");
const nonadiabatic_coupling_vector = @import("nonadiabatic_coupling_vector.zig");
const norm_preserving_interpolation = @import("norm_preserving_interpolation.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const HammesSchifferTully = hammes_schiffer_tully.HammesSchifferTully;
const NonadiabaticCouplingVector = nonadiabatic_coupling_vector.NonadiabaticCouplingVector;
const NormPreservingInterpolation = norm_preserving_interpolation.NormPreservingInterpolation;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

/// Parameters struct for all the time derivative coupling methods.
pub fn Parameters(comptime T: type) type {
    return struct {
        hst_parameters: hammes_schiffer_tully.Parameters(T),
        nacv_parameters: nonadiabatic_coupling_vector.Parameters(T),
        npi_parameters: norm_preserving_interpolation.Parameters(T)
    };
}

/// Time derivative coupling union.
pub fn DerivativeCoupling(comptime T: type) type {
    return union(enum) {
        hst: HammesSchifferTully(T),
        nacv: NonadiabaticCouplingVector(T),
        npi: NormPreservingInterpolation(T),

        /// Evaluate the time derivative coupling.
        pub fn evaluate(self: @This(), derivative_coupling: *RealMatrix(T), parameters: Parameters(T)) !void {
            switch (self) {
                .hst => |field| try field.evaluate(derivative_coupling, parameters.hst_parameters),
                .nacv => |field| try field.evaluate(derivative_coupling, parameters.nacv_parameters),
                .npi => |field| try field.evaluate(derivative_coupling, parameters.npi_parameters),
            }
        }
    };
}
