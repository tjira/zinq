//! Contracts two-electron repulsion integrals and derivatives with density matrices.

const std = @import("std");

const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;

/// Contracts two-electron repulsion integrals with density using 8-fold permutational symmetry and contiguous traversal.
pub fn contractFock(comptime T: type, F: *Matrix(T), g: Tensor(T, 4), P: Matrix(T), exch_factor: T) void {
    const nbf = F.nrow();

    std.debug.assert(g.shape[0] == nbf and g.shape[1] == nbf and g.shape[2] == nbf and g.shape[3] == nbf);

    std.debug.assert(P.shape[0] == nbf and P.shape[1] == nbf);

    const g_data = g.data;
    const P_data = P.data;
    const F_data = F.data;

    for (0..nbf) |i| {
        for (0..i + 1) |k| {
            for (0..i + 1) |j| {
                const max_l, const base_idx = .{ if (k == i) j else k, ((i * nbf + k) * nbf + j) * nbf };

                for (0..max_l + 1) |l| {
                    const v = g_data[base_idx + l];

                    const ex_v, var perms: [8][4]usize = .{ exch_factor * v, undefined };

                    perms[0] = .{ i, j, k, l };

                    var n_perms: usize = 1;

                    if (i != j) {
                        perms[1] = .{ j, i, k, l };

                        n_perms = 2;
                    }

                    if (k != l) {
                        for (0..n_perms) |idx| {
                            perms[n_perms + idx] = .{ perms[idx][0], perms[idx][1], l, k };
                        }

                        n_perms *= 2;
                    }

                    if (k != i or l != j) {
                        for (0..n_perms) |idx| {
                            perms[n_perms + idx] = .{ perms[idx][2], perms[idx][3], perms[idx][0], perms[idx][1] };
                        }

                        n_perms *= 2;
                    }

                    for (0..n_perms) |idx| {
                        const p0 = perms[idx][0];
                        const p1 = perms[idx][1];
                        const p2 = perms[idx][2];
                        const p3 = perms[idx][3];

                        F_data[p0 * nbf + p1] += P_data[p2 * nbf + p3] * v;

                        if (exch_factor != 0) {
                            F_data[p0 * nbf + p3] -= ex_v * P_data[p2 * nbf + p1];
                        }
                    }
                }
            }
        }
    }
}

/// Contracts two-electron integral derivatives with density matrices to evaluate nuclear gradient contributions.
pub fn contractGradient(comptime T: type, dg: Tensor(T, 5), coord: usize, P: Matrix(T), exch_factor: T) T {
    const nbf = P.nrow();

    var g_val: T = 0;

    for (0..nbf) |mu| for (0..nbf) |nu| for (0..nbf) |lambda| for (0..nbf) |sigma| {
        const G = dg.at(.{ coord, mu, lambda, nu, sigma }) - exch_factor * dg.at(.{ coord, mu, lambda, sigma, nu });

        g_val += 0.5 * P.at(mu, nu) * P.at(lambda, sigma) * G;
    };

    return g_val;
}

/// Contracts two-electron integral derivatives with the density matrix to form first-order Fock response matrices.
pub fn contractFockDerivative(comptime T: type, F_x: *Matrix(T), dg: Tensor(T, 5), coord: usize, P: Matrix(T), exch_factor: T) void {
    const nbf = F_x.nrow();

    for (0..nbf) |lambda| for (0..nbf) |sigma| for (0..nbf) |mu| for (0..nbf) |nu| {
        const G = dg.at(.{ coord, mu, lambda, nu, sigma }) - exch_factor * dg.at(.{ coord, mu, sigma, lambda, nu });

        F_x.ptr(mu, nu).* += P.at(lambda, sigma) * G;
    };
}

/// Contracts two-electron repulsion integrals with the response density to evaluate orbital Hessian response.
pub fn contractResponseHessian(comptime T: type, V_x: *Matrix(T), g: Tensor(T, 4), P_x: Matrix(T), exch_factor: T) void {
    const nbf = V_x.nrow();

    std.debug.assert(g.shape[0] == nbf);
    std.debug.assert(g.shape[1] == nbf);
    std.debug.assert(g.shape[2] == nbf);
    std.debug.assert(g.shape[3] == nbf);

    std.debug.assert(P_x.shape[0] == nbf);
    std.debug.assert(P_x.shape[1] == nbf);

    for (0..nbf) |mu| for (0..nbf) |nu| for (0..nbf) |lambda| for (0..nbf) |sigma| {
        V_x.ptr(lambda, sigma).* += P_x.at(mu, nu) * (g.at(.{ mu, lambda, nu, sigma }) - exch_factor * g.at(.{ mu, nu, lambda, sigma }));
    };
}
