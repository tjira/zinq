//! File with functions to perform DFT integration, including the evaluation of the exchange-correlation energy.

const std = @import("std");

const basis_set = @import("basis_set.zig");
const dft_functional = @import("dft_functional.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const ExchangeFunctional = dft_functional.ExchangeFunctional;
const CorrelationFunctional = dft_functional.CorrelationFunctional;
const BasisSet = basis_set.BasisSet;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const computeExchangeCorrelation = dft_functional.computeExchangeCorrelation;

/// Evaluate the exchange-correlation energy for a given set of points and weights, using the provided density matrix and basis set. The specific functional to use is determined by the `dft` parameter, which can be used to specify different functionals (e.g., LDA, GGA, etc.). The function returns the computed exchange-correlation energy.
pub fn evaluateXC(comptime T: type, Vxc: *RealMatrix(T), P: RealMatrix(T), basis: BasisSet(T), points: RealMatrix(T), weights: RealVector(T), exchange: ExchangeFunctional(T), correlation: CorrelationFunctional(T), generalized: bool, allocator: std.mem.Allocator) !T {
    Vxc.zero(); var Exc: T = 0; const factor: T = if (generalized) 1.0 else 2.0;
    
    var phi = try RealVector(T).init(basis.nbf(), allocator); defer phi.deinit(allocator);

    for (0..points.rows) |i| {

        var rho: T = 0;

        for (0..basis.nbf()) |j| {
            phi.ptr(j).* = basis.contracted_gaussians[j].evaluate(points.at(i, 0), points.at(i, 1), points.at(i, 2));
        }

        for (0..P.rows) |mu| {
            for (0..P.cols) |nu| {
                rho += P.at(mu, nu) * phi.at(mu) * phi.at(nu);
            }
        }

        rho *= factor;

        if (rho <= 1e-12) continue;

        const dft_result = computeExchangeCorrelation(T, exchange, correlation, rho);

        Exc += rho * dft_result.eps_xc * weights.at(i);

        for (0..Vxc.rows) |mu| {
            for (0..Vxc.cols) |nu| {
                Vxc.ptr(mu, nu).* += dft_result.v_xc * phi.at(mu) * phi.at(nu) * weights.at(i);
            }
        }
    }

    return Exc;
}

