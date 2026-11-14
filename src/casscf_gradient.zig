//! Evaluates analytic nuclear gradients for state-averaged and state-specific CASSCF and RI-CASSCF methods.

const std = @import("std");

const Allocator = std.mem.Allocator;

const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const MetricDecomposition = @import("density_fitting.zig").MetricDecomposition;
const OrbitalPair = @import("casscf_orbital.zig").OrbitalPair;
const ScalarDual = @import("dual.zig").ScalarDual;
const Tensor = @import("tensor.zig").Tensor;

const ao2mo_1e = @import("integral_transform.zig").ao2mo_1e;
const ao2mo_2e = @import("integral_transform.zig").ao2mo_2e;
const ao2mo_3c = @import("integral_transform.zig").ao2mo_3c;
const buildActRdm = @import("reduced_density_matrix.zig").buildActRdm;
const buildGeneralizedFock = @import("casscf_orbital.zig").buildGeneralizedFock;
const buildMetricDerivKernel = @import("density_fitting.zig").buildMetricDerivKernel;
const buildNonRedundantPairs = @import("casscf_orbital.zig").buildNonRedundantPairs;
const buildOrbitalHessian = @import("casscf_orbital.zig").buildOrbitalHessian;
const buildStateAveragedActRdm = @import("reduced_density_matrix.zig").buildStateAveragedActRdm;
const evaluateFittingTensorDeriv = @import("density_fitting.zig").evaluateFittingTensorDeriv;
const evaluateMetricDeriv = @import("density_fitting.zig").evaluateMetricDeriv;
const gemm = @import("linear_algebra.zig").gemm;
const getrf = @import("linear_algebra.zig").getrf;
const getrs = @import("linear_algebra.zig").getrs;
const slater = @import("slater_condon.zig").slater;
const transformInts = @import("integral_transform.zig").transformInts;

/// Computes analytic nuclear gradients of CASSCF state energies using dual number differentiation and CP-MCSCF response.
pub fn calculateCanonicalCasscfGradient(comptime T: type, hfres: HartreeFockResult(T), C: Matrix(T), C_ci: Matrix(T), dets: []const []const usize, weights: []const T, ncore: usize, nact: usize, states: []const u32, gpa: Allocator) ![]Matrix(T) {
    var grads = try gpa.alloc(Matrix(T), states.len);
    errdefer gpa.free(grads);

    var grads_done: usize = 0;
    errdefer for (0..grads_done) |k| grads[k].deinit(gpa);

    for (0..states.len) |idx| {
        grads[idx] = try hfres.ints.sys.nuclearRepulsionGradient(gpa);

        grads_done += 1;
    }

    const nbf, const nmo = .{ C.nrow(), C.ncol() };

    const natoms = hfres.ints.sys.atoms.len;

    const pairs = try buildNonRedundantPairs(ncore, nact, nmo, gpa);
    defer gpa.free(pairs);

    const nrot = pairs.len;

    var D_act_sa, var d_act_sa = try buildStateAveragedActRdm(T, dets, C_ci, weights, ncore, nact, gpa);

    defer {
        D_act_sa.deinit(gpa);
        d_act_sa.deinit(gpa);
    }

    var Z_mat = try gpa.alloc(Matrix(T), states.len);

    defer {
        for (0..states.len) |k| Z_mat[k].deinit(gpa);

        gpa.free(Z_mat);
    }

    if (weights.len > 1 and nrot > 0) {
        var h_mo_conv = try Matrix(T).init(nmo, nmo, gpa);
        defer h_mo_conv.deinit(gpa);

        ao2mo_1e(T, &h_mo_conv, hfres.ints.H.?, C);

        var g_mo_conv = try Tensor(T, 4).init(.{ nmo, nmo, nmo, nmo }, gpa);
        defer g_mo_conv.deinit(gpa);

        try ao2mo_2e(T, &g_mo_conv, hfres.ints.g.?, C, gpa);

        var A_sa = try Matrix(T).initZero(nrot, nrot, gpa);
        defer A_sa.deinit(gpa);

        try buildOrbitalHessian(T, &A_sa, h_mo_conv, g_mo_conv, D_act_sa, d_act_sa, pairs, ncore, gpa);

        var F_gen_s = try Matrix(T).initZero(nmo, nmo, gpa);
        defer F_gen_s.deinit(gpa);

        const ipiv = try gpa.alloc(i32, nrot);
        defer gpa.free(ipiv);

        var A_LU = try Matrix(T).init(nrot, nrot, gpa);
        defer A_LU.deinit(gpa);

        @memcpy(A_LU.data, A_sa.data);

        try getrf(T, &A_LU, ipiv);

        for (states, 0..) |s, s_idx| {
            Z_mat[s_idx] = try Matrix(T).initZero(nrot, 1, gpa);

            var D_s, var d_s = try buildActRdm(T, dets, C_ci, s, ncore, nact, gpa);

            defer {
                D_s.deinit(gpa);
                d_s.deinit(gpa);
            }

            try buildGeneralizedFock(T, &F_gen_s, h_mo_conv, g_mo_conv, D_s, d_s, ncore, gpa);

            var neg_g = try Matrix(T).init(nrot, 1, gpa);
            defer neg_g.deinit(gpa);

            for (pairs, 0..) |pair, i| {
                neg_g.ptr(i, 0).* = -2 * (F_gen_s.at(pair.p, pair.q) - F_gen_s.at(pair.q, pair.p));
            }

            try getrs(T, &Z_mat[s_idx], A_LU, ipiv, neg_g);
        }
    }

    if (weights.len == 1 or nrot == 0) {
        for (states, 0..) |_, s_idx| {
            Z_mat[s_idx] = try Matrix(T).initZero(nrot, 1, gpa);
        }
    }

    var H_d = try Matrix(ScalarDual(T)).init(nbf, nbf, gpa);
    defer H_d.deinit(gpa);

    var g_d = try Tensor(ScalarDual(T), 4).init(.{ nbf, nbf, nbf, nbf }, gpa);
    defer g_d.deinit(gpa);

    var C_d = try Matrix(ScalarDual(T)).init(nbf, nmo, gpa);
    defer C_d.deinit(gpa);

    var dS_ao = try Matrix(T).init(nbf, nbf, gpa);
    defer dS_ao.deinit(gpa);

    var tmp = try Matrix(T).init(nbf, nmo, gpa);
    defer tmp.deinit(gpa);

    var S_mo = try Matrix(T).init(nmo, nmo, gpa);
    defer S_mo.deinit(gpa);

    var dC = try Matrix(T).init(nbf, nmo, gpa);
    defer dC.deinit(gpa);

    for (0..3 * natoms) |coord| {
        for (0..nbf) |mu| for (0..nbf) |nu| {
            dS_ao.ptr(mu, nu).* = hfres.ints.dS.?.at(.{ coord, mu, nu });
        };

        gemm(T, &tmp, dS_ao, C, 1, 0, false, false);

        gemm(T, &S_mo, C, tmp, 1, 0, true, false);

        gemm(T, &dC, C, S_mo, -0.5, 0, false, false);

        for (0..nbf) |mu| for (0..nmo) |p| {
            C_d.ptr(mu, p).* = ScalarDual(T).init(C.at(mu, p), dC.at(mu, p));
        };

        for (0..nbf) |mu| for (0..nbf) |nu| {
            H_d.ptr(mu, nu).* = ScalarDual(T).init(hfres.ints.H.?.at(mu, nu), hfres.ints.dH.?.at(.{ coord, mu, nu }));
        };

        for (0..nbf) |mu| for (0..nbf) |lambda| for (0..nbf) |nu| for (0..nbf) |sigma| {
            const indices = .{ mu, lambda, nu, sigma };

            const val = hfres.ints.g.?.at(indices);

            const der = hfres.ints.dg.?.at(.{ coord, mu, lambda, nu, sigma });

            g_d.ptr(indices).* = ScalarDual(T).init(val, der);
        };

        var h_ms_d, var g_ms_d = try transformInts(ScalarDual(T), C_d, H_d, g_d, false, gpa);

        defer {
            h_ms_d.deinit(gpa);
            g_ms_d.deinit(gpa);
        }

        for (0..dets.len) |m| for (m..dets.len) |n| {
            const factor: T = if (m == n) 1 else 2;

            const val = slater(ScalarDual(T), dets[m], dets[n], h_ms_d, g_ms_d);

            for (states, 0..) |s, s_idx| {
                grads[s_idx].ptr(coord / 3, coord % 3).* += factor * C_ci.at(m, s) * C_ci.at(n, s) * val.der;
            }
        };

        if (weights.len > 1 and nrot > 0) {
            var dh_mo = try Matrix(T).init(nmo, nmo, gpa);
            defer dh_mo.deinit(gpa);

            var dg_mo = try Tensor(T, 4).init(.{ nmo, nmo, nmo, nmo }, gpa);
            defer dg_mo.deinit(gpa);

            for (0..nmo) |p| for (0..nmo) |q| {
                dh_mo.ptr(p, q).* = h_ms_d.at(2 * p, 2 * q).der;
            };

            for (0..nmo) |p| for (0..nmo) |q| for (0..nmo) |r| for (0..nmo) |s| {
                dg_mo.ptr(.{ p, q, r, s }).* = g_ms_d.at(.{ 2 * p, 2 * q, 2 * r, 2 * s }).der;
            };

            var dF_gen = try Matrix(T).initZero(nmo, nmo, gpa);
            defer dF_gen.deinit(gpa);

            try buildGeneralizedFock(T, &dF_gen, dh_mo, dg_mo, D_act_sa, d_act_sa, ncore, gpa);

            for (states, 0..) |_, s_idx| {
                var resp: T = 0;

                for (pairs, 0..) |pair, k| {
                    const B_k = 2 * (dF_gen.at(pair.p, pair.q) - dF_gen.at(pair.q, pair.p));

                    resp += Z_mat[s_idx].at(k, 0) * B_k;
                }

                grads[s_idx].ptr(coord / 3, coord % 3).* -= resp;
            }
        }
    }

    return grads;
}

/// Computes analytic nuclear gradients of CASSCF state energies using density fitting and dual numbers.
pub fn calculateRiCasscfGradient(comptime T: type, hfres: HartreeFockResult(T), C: Matrix(T), C_ci: Matrix(T), dets: []const []const usize, weights: []const T, ncore: usize, nact: usize, states: []const u32, gpa: Allocator) ![]Matrix(T) {
    var grads = try gpa.alloc(Matrix(T), states.len);
    errdefer gpa.free(grads);

    var grads_done: usize = 0;
    errdefer for (0..grads_done) |k| grads[k].deinit(gpa);

    for (0..states.len) |idx| {
        grads[idx] = try hfres.ints.sys.nuclearRepulsionGradient(gpa);

        grads_done += 1;
    }

    const J, const B = .{ hfres.ints.J orelse unreachable, hfres.ints.B.? };

    const dI = hfres.ints.dI orelse unreachable;
    const dJ = hfres.ints.dJ orelse unreachable;

    const I_3c = hfres.ints.I orelse unreachable;

    const nbf, const naux = .{ B.shape[0], B.shape[2] };

    const nmo = C.ncol();

    const natoms, const nocc = .{ hfres.ints.sys.atoms.len, ncore + nact };

    const nsp = 2 * nocc;

    const pairs = try buildNonRedundantPairs(ncore, nact, nmo, gpa);
    defer gpa.free(pairs);

    const nrot = pairs.len;

    var D_act_sa, var d_act_sa = try buildStateAveragedActRdm(T, dets, C_ci, weights, ncore, nact, gpa);

    defer {
        D_act_sa.deinit(gpa);
        d_act_sa.deinit(gpa);
    }

    var Z_mat = try gpa.alloc(Matrix(T), states.len);

    defer {
        for (0..states.len) |k| Z_mat[k].deinit(gpa);

        gpa.free(Z_mat);
    }

    if (weights.len > 1 and nrot > 0) {
        var B_mo = try Tensor(T, 3).initZero(.{ nmo, nmo, naux }, gpa);
        defer B_mo.deinit(gpa);

        try ao2mo_3c(T, &B_mo, B, C, gpa);

        var h_mo = try Matrix(T).init(nmo, nmo, gpa);
        defer h_mo.deinit(gpa);

        ao2mo_1e(T, &h_mo, hfres.ints.H.?, C);

        var g_mo = try Tensor(T, 4).init(.{ nmo, nmo, nmo, nmo }, gpa);
        defer g_mo.deinit(gpa);

        for (0..nmo) |p| for (0..nmo) |r| for (0..nmo) |q| for (0..nmo) |s| {
            var sum: T = 0;

            for (0..naux) |P| {
                sum += B_mo.at(.{ p, r, P }) * B_mo.at(.{ q, s, P });
            }

            g_mo.ptr(.{ p, q, r, s }).* = sum;
        };

        var A_sa = try Matrix(T).initZero(nrot, nrot, gpa);
        defer A_sa.deinit(gpa);

        try buildOrbitalHessian(T, &A_sa, h_mo, g_mo, D_act_sa, d_act_sa, pairs, ncore, gpa);

        var F_gen_s = try Matrix(T).initZero(nmo, nmo, gpa);
        defer F_gen_s.deinit(gpa);

        const ipiv = try gpa.alloc(i32, nrot);
        defer gpa.free(ipiv);

        var A_LU = try Matrix(T).init(nrot, nrot, gpa);
        defer A_LU.deinit(gpa);

        @memcpy(A_LU.data, A_sa.data);

        try getrf(T, &A_LU, ipiv);

        for (states, 0..) |s, s_idx| {
            Z_mat[s_idx] = try Matrix(T).initZero(nrot, 1, gpa);

            var D_s, var d_s = try buildActRdm(T, dets, C_ci, s, ncore, nact, gpa);

            defer {
                D_s.deinit(gpa);
                d_s.deinit(gpa);
            }

            try buildGeneralizedFock(T, &F_gen_s, h_mo, g_mo, D_s, d_s, ncore, gpa);

            var neg_g = try Matrix(T).init(nrot, 1, gpa);
            defer neg_g.deinit(gpa);

            for (pairs, 0..) |pair, i| {
                neg_g.ptr(i, 0).* = -2 * (F_gen_s.at(pair.p, pair.q) - F_gen_s.at(pair.q, pair.p));
            }

            try getrs(T, &Z_mat[s_idx], A_LU, ipiv, neg_g);
        }
    }

    if (weights.len == 1 or nrot == 0) {
        for (states, 0..) |_, s_idx| {
            Z_mat[s_idx] = try Matrix(T).initZero(nrot, 1, gpa);
        }
    }

    var decomp = try MetricDecomposition(T).init(J, gpa);
    defer decomp.deinit(gpa);

    var M_kernel = try buildMetricDerivKernel(T, decomp.W, gpa);
    defer M_kernel.deinit(gpa);

    var dJ_inv_half = try Matrix(T).init(naux, naux, gpa);
    defer dJ_inv_half.deinit(gpa);

    var dB_ao = try Tensor(T, 3).initZero(B.shape, gpa);
    defer dB_ao.deinit(gpa);

    var H_d = try Matrix(ScalarDual(T)).init(nbf, nbf, gpa);
    defer H_d.deinit(gpa);

    var B_d = try Tensor(ScalarDual(T), 3).init(B.shape, gpa);
    defer B_d.deinit(gpa);

    var C_d = try Matrix(ScalarDual(T)).init(nbf, nmo, gpa);
    defer C_d.deinit(gpa);

    var dS_ao = try Matrix(T).init(nbf, nbf, gpa);
    defer dS_ao.deinit(gpa);

    var tmp = try Matrix(T).init(nbf, nmo, gpa);
    defer tmp.deinit(gpa);

    var S_mo = try Matrix(T).init(nmo, nmo, gpa);
    defer S_mo.deinit(gpa);

    var dC = try Matrix(T).init(nbf, nmo, gpa);
    defer dC.deinit(gpa);

    var h_mo_d = try Matrix(ScalarDual(T)).init(nmo, nmo, gpa);
    defer h_mo_d.deinit(gpa);

    var B_mo_d = try Tensor(ScalarDual(T), 3).initZero(.{ nmo, nmo, naux }, gpa);
    defer B_mo_d.deinit(gpa);

    var h_ms_d = try Matrix(ScalarDual(T)).initZero(nsp, nsp, gpa);
    defer h_ms_d.deinit(gpa);

    var g_ms_d = try Tensor(ScalarDual(T), 4).initZero(.{ nsp, nsp, nsp, nsp }, gpa);
    defer g_ms_d.deinit(gpa);

    var dh_mo = if (weights.len > 1 and nrot > 0) try Matrix(T).init(nmo, nmo, gpa) else undefined;
    defer if (weights.len > 1 and nrot > 0) dh_mo.deinit(gpa);

    var dg_mo = if (weights.len > 1 and nrot > 0) try Tensor(T, 4).init(.{ nmo, nmo, nmo, nmo }, gpa) else undefined;
    defer if (weights.len > 1 and nrot > 0) dg_mo.deinit(gpa);

    var dF_gen = if (weights.len > 1 and nrot > 0) try Matrix(T).initZero(nmo, nmo, gpa) else undefined;
    defer if (weights.len > 1 and nrot > 0) dF_gen.deinit(gpa);

    for (0..3 * natoms) |coord| {
        for (0..nbf) |mu| for (0..nbf) |nu| {
            dS_ao.ptr(mu, nu).* = hfres.ints.dS.?.at(.{ coord, mu, nu });
        };

        gemm(T, &tmp, dS_ao, C, 1, 0, false, false);

        gemm(T, &S_mo, C, tmp, 1, 0, true, false);

        gemm(T, &dC, C, S_mo, -0.5, 0, false, false);

        for (0..nbf) |mu| for (0..nmo) |p| {
            C_d.ptr(mu, p).* = ScalarDual(T).init(C.at(mu, p), dC.at(mu, p));
        };

        for (0..nbf) |mu| for (0..nbf) |nu| {
            H_d.ptr(mu, nu).* = ScalarDual(T).init(hfres.ints.H.?.at(mu, nu), hfres.ints.dH.?.at(.{ coord, mu, nu }));
        };

        try evaluateMetricDeriv(T, &dJ_inv_half, dJ.slice(coord).asMatrix(), decomp.U, M_kernel, gpa);

        evaluateFittingTensorDeriv(T, &dB_ao, dI.slice(coord), I_3c, decomp.J_inv_half, dJ_inv_half);

        for (0..nbf) |mu| for (0..nbf) |nu| for (0..naux) |P| {
            B_d.ptr(.{ mu, nu, P }).* = ScalarDual(T).init(B.at(.{ mu, nu, P }), dB_ao.at(.{ mu, nu, P }));
        };

        ao2mo_1e(ScalarDual(T), &h_mo_d, H_d, C_d);

        try ao2mo_3c(ScalarDual(T), &B_mo_d, B_d, C_d, gpa);

        for (0..nocc) |p| for (0..nocc) |q| {
            h_ms_d.ptr(2 * p + 0, 2 * q + 0).* = h_mo_d.at(p, q);
            h_ms_d.ptr(2 * p + 1, 2 * q + 1).* = h_mo_d.at(p, q);
        };

        for (0..nocc) |p| for (0..nocc) |r| {
            for (0..nocc) |q| for (0..nocc) |s| {
                var sum = ScalarDual(T).init(0, 0);

                for (0..naux) |P| {
                    sum = sum.add(B_mo_d.at(.{ p, r, P }).mul(B_mo_d.at(.{ q, s, P })));
                }

                g_ms_d.ptr(.{ 2 * p + 0, 2 * q + 0, 2 * r + 0, 2 * s + 0 }).* = sum;
                g_ms_d.ptr(.{ 2 * p + 0, 2 * q + 1, 2 * r + 0, 2 * s + 1 }).* = sum;
                g_ms_d.ptr(.{ 2 * p + 1, 2 * q + 0, 2 * r + 1, 2 * s + 0 }).* = sum;
                g_ms_d.ptr(.{ 2 * p + 1, 2 * q + 1, 2 * r + 1, 2 * s + 1 }).* = sum;
            };
        };

        for (0..dets.len) |m| for (m..dets.len) |n| {
            const factor: T = if (m == n) 1 else 2;

            const val = slater(ScalarDual(T), dets[m], dets[n], h_ms_d, g_ms_d);

            for (states, 0..) |s, s_idx| {
                grads[s_idx].ptr(coord / 3, coord % 3).* += factor * C_ci.at(m, s) * C_ci.at(n, s) * val.der;
            }
        };

        if (weights.len > 1 and nrot > 0) {
            for (0..nmo) |p| for (0..nmo) |q| {
                dh_mo.ptr(p, q).* = h_mo_d.at(p, q).der;
            };

            for (0..nmo) |p| for (0..nmo) |r| for (0..nmo) |q| for (0..nmo) |s| {
                var sum: T = 0;

                for (0..naux) |P| {
                    const b_pr = B_mo_d.at(.{ p, r, P });
                    const b_qs = B_mo_d.at(.{ q, s, P });

                    sum += b_pr.der * b_qs.val + b_pr.val * b_qs.der;
                }

                dg_mo.ptr(.{ p, q, r, s }).* = sum;
            };

            @memset(dF_gen.data, 0);

            try buildGeneralizedFock(T, &dF_gen, dh_mo, dg_mo, D_act_sa, d_act_sa, ncore, gpa);

            for (states, 0..) |_, s_idx| {
                var resp: T = 0;

                for (pairs, 0..) |pair, k| {
                    const B_k = 2 * (dF_gen.at(pair.p, pair.q) - dF_gen.at(pair.q, pair.p));

                    resp += Z_mat[s_idx].at(k, 0) * B_k;
                }

                grads[s_idx].ptr(coord / 3, coord % 3).* -= resp;
            }
        }
    }

    return grads;
}
