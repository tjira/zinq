//! File with functions to perform DFT integration, including the evaluation of the exchange-correlation energy.

const std = @import("std");

const basis_set = @import("basis_set.zig");
const dft_functional = @import("dft_functional.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const BasisSet = basis_set.BasisSet;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const computeExchangeCorrelationArray = dft_functional.computeExchangeCorrelationArray;

/// Density functional theory context with necessary arguments.
pub fn DensityContext(comptime T: type) type {
    return struct {
        basis: BasisSet(T), points: RealMatrix(T), weights: RealVector(T), exchange: ?dft_functional.ExchangeFunctional(T), correlation: ?dft_functional.CorrelationFunctional(T), generalized: bool
    };
}

/// Evaluate the exchange-correlation energy for a given set of points and weights, using the provided density matrix and basis set. The specific functional to use is determined by the `dft` parameter, which can be used to specify different functionals (e.g., LDA, GGA, etc.). The function returns the computed exchange-correlation energy.
pub fn evaluateXC(comptime T: type, Vxc: *RealMatrix(T), P: RealMatrix(T), context: DensityContext(T), allocator: std.mem.Allocator) !T {
    Vxc.zero(); var Exc: T = 0; const factor: T = if (context.generalized) 1.0 else 2.0;

    const basis = context.basis;
    const points = context.points; const weights = context.weights;
    const exchange = context.exchange; const correlation = context.correlation;
    
    var phi = try RealMatrix(T).initZero(points.rows, basis.nbf(), allocator); defer phi.deinit(allocator);

    var rho = try RealVector(T).initZero(points.rows, allocator); defer rho.deinit(allocator);
    var eps = try RealVector(T).initZero(points.rows, allocator); defer eps.deinit(allocator);
    var vxc = try RealVector(T).initZero(points.rows, allocator); defer vxc.deinit(allocator);

    for (0..rho.len) |i| {

        for (0..basis.nbf()) |j| {
            phi.ptr(i, j).* = basis.contracted_gaussians[j].evaluate(points.at(i, 0), points.at(i, 1), points.at(i, 2));
        }

        for (0..P.rows) |mu| {
            for (0..P.cols) |nu| {
                rho.ptr(i).* += P.at(mu, nu) * phi.at(i, mu) * phi.at(i, nu);
            }
        }

        rho.ptr(i).* *= factor;
    }

    computeExchangeCorrelationArray(T, &eps, exchange, correlation, rho, 0);
    computeExchangeCorrelationArray(T, &vxc, exchange, correlation, rho, 1);

    for (0..points.rows) |i| {

        Exc += rho.at(i) * eps.at(i) * weights.at(i);

        for (0..Vxc.rows) |mu| for (0..Vxc.cols) |nu| {
            Vxc.ptr(mu, nu).* += vxc.at(i) * phi.at(i, mu) * phi.at(i, nu) * weights.at(i);
        };
    }

    return Exc;
}

/// Evaluate the spatial exchange-correlation kernel matrix (ia|fxc|jb) for TDDFT.
pub fn evaluateXCKernel(comptime T: type, K: *RealMatrix(T), P: RealMatrix(T), C: RealMatrix(T), context: DensityContext(T), nocc: usize, allocator: std.mem.Allocator) !void {
    K.zero(); const factor: T = if (context.generalized) 1.0 else 2.0; const nvir = context.basis.nbf() - nocc;

    const basis = context.basis;
    const points = context.points; const weights = context.weights;
    const exchange = context.exchange; const correlation = context.correlation;

    var chi = try RealMatrix(T).initZero(points.rows, basis.nbf(), allocator); defer chi.deinit(allocator);
    var phi = try RealMatrix(T).initZero(points.rows, basis.nbf(), allocator); defer phi.deinit(allocator);
    
    var rho = try RealVector(T).initZero(points.rows, allocator); defer rho.deinit(allocator);
    var fxc = try RealVector(T).initZero(points.rows, allocator); defer fxc.deinit(allocator);

    for (0..rho.len) |i| {

        for (0..basis.nbf()) |mu| {
            chi.ptr(i, mu).* = basis.contracted_gaussians[mu].evaluate(points.at(i, 0), points.at(i, 1), points.at(i, 2));
        }

        for (0..P.rows) |mu| {
            for (0..P.cols) |nu| {
                rho.ptr(i).* += P.at(mu, nu) * chi.at(i, mu) * chi.at(i, nu);
            }
        }

        rho.ptr(i).* *= factor;

        for (0..basis.nbf()) |p| for (0..basis.nbf()) |mu| {
            phi.ptr(i, p).* += C.at(mu, p) * chi.at(i, mu);
        };
    }

    computeExchangeCorrelationArray(T, &fxc, exchange, correlation, rho, 2);

    for (0..points.rows) |g| {

        if (fxc.at(g) == 0) continue;

        for (0..nocc) |i| for (0..nvir) |a| for (0..nocc) |j| for (0..nvir) |b| {

            const ia = i * nvir + a;
            const jb = j * nvir + b;
            
            if (jb < ia) continue;

            const phi_ia = phi.at(g, i) * phi.at(g, nocc + a);
            const phi_jb = phi.at(g, j) * phi.at(g, nocc + b);

            const val = weights.at(g) * fxc.at(g) * phi_ia * phi_jb;
            
            K.ptr(ia, jb).* += val; if (ia != jb) K.ptr(jb, ia).* += val;
        };
    }
}
