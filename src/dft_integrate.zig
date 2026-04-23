//! File with functions to perform DFT integration, including the evaluation of the exchange-correlation energy.

const std = @import("std");

const basis_set = @import("basis_set.zig");
const dft_functional = @import("dft_functional.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const BasisSet = basis_set.BasisSet;
const DensityFunctional = dft_functional.DensityFunctional;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const computeExchangeCorrelationArray = dft_functional.computeExchangeCorrelationArray;
const getFunctionalFamily = dft_functional.getFunctionalFamily;

/// Density functional theory context with necessary arguments.
pub fn DensityIntegrateContext(comptime T: type) type {
    return struct { basis: BasisSet(T), points: RealMatrix(T), weights: RealVector(T), functional: DensityFunctional(T), generalized: bool };
}

/// Evaluate the exchange-correlation energy for a given set of points and weights, using the provided density matrix and basis set. The specific functional to use is determined by the `dft` parameter, which can be used to specify different functionals (e.g., LDA, GGA, etc.). The function returns the computed exchange-correlation energy.
pub fn evaluateXC(comptime T: type, Vxc: *RealMatrix(T), P: RealMatrix(T), context: DensityIntegrateContext(T), allocator: std.mem.Allocator) !T {
    Vxc.zero();
    var Exc: T = 0;
    const factor: T = if (context.generalized) 1.0 else 2.0;

    const basis = context.basis;
    const points = context.points;
    const weights = context.weights;
    const functional = context.functional;

    if (functional.exchange) |x| if (try getFunctionalFamily(x) != .lda) {
        std.log.err("CURRENTLY ONLY LDA FUNCTIONALS ARE IMPLEMENTED", .{});
        return error.InputError;
    };

    if (functional.correlation) |c| if (try getFunctionalFamily(c) != .lda) {
        std.log.err("CURRENTLY ONLY LDA FUNCTIONALS ARE IMPLEMENTED", .{});
        return error.InputError;
    };

    var phi = try RealMatrix(T).initZero(points.rows, basis.nbf(), allocator);
    defer phi.deinit(allocator);

    var rho = try RealVector(T).initZero(points.rows, allocator);
    defer rho.deinit(allocator);

    var eps_x = try RealVector(T).initZero(points.rows, allocator);
    defer eps_x.deinit(allocator);

    var eps_c = try RealVector(T).initZero(points.rows, allocator);
    defer eps_c.deinit(allocator);

    var vxc_x = try RealVector(T).initZero(points.rows, allocator);
    defer vxc_x.deinit(allocator);

    var vxc_c = try RealVector(T).initZero(points.rows, allocator);
    defer vxc_c.deinit(allocator);

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

    try computeExchangeCorrelationArray(T, &eps_x, &eps_c, functional, rho, 0);
    try computeExchangeCorrelationArray(T, &vxc_x, &vxc_c, functional, rho, 1);

    for (0..points.rows) |i| {
        const eps = eps_x.at(i) + eps_c.at(i);
        const vxc = vxc_x.at(i) + vxc_c.at(i);

        Exc += rho.at(i) * eps * weights.at(i);

        for (0..Vxc.rows) |mu| for (0..Vxc.cols) |nu| {
            Vxc.ptr(mu, nu).* += vxc * phi.at(i, mu) * phi.at(i, nu) * weights.at(i);
        };
    }

    return Exc;
}

/// Evaluate the spatial exchange-correlation kernel matrix (ia|fxc|jb) for TDDFT.
pub fn evaluateXCKernel(comptime T: type, K: *RealMatrix(T), P: RealMatrix(T), C: RealMatrix(T), context: DensityIntegrateContext(T), nocc: usize, allocator: std.mem.Allocator) !void {
    K.zero();
    const factor: T = if (context.generalized) 1.0 else 2.0;
    const nvir = context.basis.nbf() - nocc;

    const basis = context.basis;
    const points = context.points;
    const weights = context.weights;
    const functional = context.functional;

    if (functional.exchange) |x| if (try getFunctionalFamily(x) != .lda) {
        std.log.err("CURRENTLY ONLY LDA FUNCTIONALS ARE IMPLEMENTED", .{});
        return error.InputError;
    };

    if (functional.correlation) |c| if (try getFunctionalFamily(c) != .lda) {
        std.log.err("CURRENTLY ONLY LDA FUNCTIONALS ARE IMPLEMENTED", .{});
        return error.InputError;
    };

    var chi = try RealMatrix(T).initZero(points.rows, basis.nbf(), allocator);
    defer chi.deinit(allocator);

    var phi = try RealMatrix(T).initZero(points.rows, basis.nbf(), allocator);
    defer phi.deinit(allocator);

    var rho = try RealVector(T).initZero(points.rows, allocator);
    defer rho.deinit(allocator);

    var fxc_x = try RealVector(T).initZero(points.rows, allocator);
    defer fxc_x.deinit(allocator);

    var fxc_c = try RealVector(T).initZero(points.rows, allocator);
    defer fxc_c.deinit(allocator);

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

    try computeExchangeCorrelationArray(T, &fxc_x, &fxc_c, functional, rho, 2);

    for (0..points.rows) |g| {
        const fxc = fxc_x.at(g) + fxc_c.at(g);

        if (fxc == 0) continue;

        for (0..nocc) |i| for (0..nvir) |a| for (0..nocc) |j| for (0..nvir) |b| {
            const ia = i * nvir + a;
            const jb = j * nvir + b;

            if (jb < ia) continue;

            const phi_ia = phi.at(g, i) * phi.at(g, nocc + a);
            const phi_jb = phi.at(g, j) * phi.at(g, nocc + b);

            const val = weights.at(g) * fxc * phi_ia * phi_jb;

            K.ptr(ia, jb).* += val;
            if (ia != jb) K.ptr(jb, ia).* += val;
        };
    }
}
