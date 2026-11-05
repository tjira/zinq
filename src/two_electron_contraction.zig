//! Contracts two-electron repulsion integrals and derivatives with density matrices.

const std = @import("std");

const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;

/// Contracts two-electron repulsion integrals with the one-particle density matrix to form the Coulomb matrix.
pub fn contractCoulomb(comptime T: type, J: *Matrix(T), g: Tensor(T, 4), P: Matrix(T)) void {
    const nbf = J.nrow();

    std.debug.assert(g.shape[0] == nbf);
    std.debug.assert(g.shape[1] == nbf);
    std.debug.assert(g.shape[2] == nbf);
    std.debug.assert(g.shape[3] == nbf);
    std.debug.assert(P.shape[0] == nbf);
    std.debug.assert(P.shape[1] == nbf);

    for (0..nbf) |mu| for (0..nbf) |nu| for (0..nbf) |lambda| for (0..nbf) |sigma| {
        J.ptr(mu, nu).* += P.at(lambda, sigma) * g.at(.{ mu, lambda, nu, sigma });
    };
}

/// Contracts two-electron repulsion integrals with the one-particle density matrix to form the exchange matrix.
pub fn contractExchange(comptime T: type, K: *Matrix(T), g: Tensor(T, 4), P: Matrix(T)) void {
    const nbf = K.nrow();

    std.debug.assert(g.shape[0] == nbf);
    std.debug.assert(g.shape[1] == nbf);
    std.debug.assert(g.shape[2] == nbf);
    std.debug.assert(g.shape[3] == nbf);
    std.debug.assert(P.shape[0] == nbf);
    std.debug.assert(P.shape[1] == nbf);

    for (0..nbf) |mu| for (0..nbf) |nu| for (0..nbf) |lambda| for (0..nbf) |sigma| {
        K.ptr(mu, nu).* += P.at(lambda, sigma) * g.at(.{ mu, lambda, sigma, nu });
    };
}

/// Evaluates Coulomb and exchange contractions directly to accumulate two-electron contributions into Fock matrix.
pub fn contractTwoElectronFock(comptime T: type, F: *Matrix(T), g: Tensor(T, 4), P: Matrix(T), exch_factor: T) void {
    const nbf = F.nrow();

    std.debug.assert(g.shape[0] == nbf);
    std.debug.assert(g.shape[1] == nbf);
    std.debug.assert(g.shape[2] == nbf);
    std.debug.assert(g.shape[3] == nbf);
    std.debug.assert(P.shape[0] == nbf);
    std.debug.assert(P.shape[1] == nbf);

    for (0..nbf) |mu| for (0..nbf) |nu| for (0..nbf) |lambda| for (0..nbf) |sigma| {
        const j_val = g.at(.{ mu, lambda, nu, sigma });
        const k_val = g.at(.{ mu, lambda, sigma, nu });

        F.ptr(mu, nu).* += P.at(lambda, sigma) * (j_val - exch_factor * k_val);
    };
}

/// Contracts two-electron integral derivatives with density matrices to evaluate nuclear gradient contributions.
pub fn contractTwoElectronGradient(comptime T: type, dg: Tensor(T, 5), coord: usize, P: Matrix(T), exch_factor: T) T {
    const nbf = P.nrow();

    var g_val: T = 0;

    for (0..nbf) |mu| for (0..nbf) |nu| for (0..nbf) |lambda| for (0..nbf) |sigma| {
        const dg1 = dg.at(.{ coord, mu, lambda, nu, sigma });
        const dg2 = dg.at(.{ coord, mu, lambda, sigma, nu });

        g_val += 0.5 * P.at(mu, nu) * P.at(lambda, sigma) * (dg1 - exch_factor * dg2);
    };

    return g_val;
}

/// Contracts two-electron integral derivatives with the density matrix to form first-order Fock response matrices.
pub fn contractTwoElectronResponse(comptime T: type, F_x: *Matrix(T), dg: Tensor(T, 5), coord: usize, P: Matrix(T), exch_factor: T) void {
    const nbf = F_x.nrow();

    for (0..nbf) |lambda| for (0..nbf) |sigma| for (0..nbf) |mu| for (0..nbf) |nu| {
        const coul = dg.at(.{ coord, lambda, mu, sigma, nu });
        const exch = dg.at(.{ coord, lambda, sigma, mu, nu });

        F_x.ptr(mu, nu).* += P.at(lambda, sigma) * (coul - exch_factor * exch);
    };
}