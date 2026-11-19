//! Evaluates analytic nuclear gradients for state-averaged and state-specific CASSCF and RI-CASSCF methods.

const std = @import("std");

const Allocator = std.mem.Allocator;

const CasSpace = @import("casscf_orbital.zig").CasSpace;
const CasState = @import("casscf_orbital.zig").CasState;
const DualMatrix = @import("dual.zig").DualMatrix;
const DualTensor = @import("dual.zig").DualTensor;
const FittingDerivativeWorkspace = @import("density_fitting.zig").FittingDerivativeWorkspace;
const GeneralizedFockRiWorkspace = @import("casscf_orbital.zig").GeneralizedFockRiWorkspace;
const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const NonRelaxedOrbitalWorkspace = @import("orbital_response.zig").NonRelaxedOrbitalWorkspace;
const OrbitalPair = @import("casscf_orbital.zig").OrbitalPair;
const Tensor = @import("tensor.zig").Tensor;

const ao2mo_1e = @import("integral_transform.zig").ao2mo_1e;
const ao2mo_2e = @import("integral_transform.zig").ao2mo_2e;
const ao2mo_3c = @import("integral_transform.zig").ao2mo_3c;
const b2ms = @import("integral_transform.zig").b2ms;
const buildActRdm = @import("reduced_density_matrix.zig").buildActRdm;
const buildGeneralizedFock = @import("casscf_orbital.zig").buildGeneralizedFock;
const buildGeneralizedFockDerivRi = @import("casscf_orbital.zig").buildGeneralizedFockDerivRi;
const buildNonRedundantPairs = @import("casscf_orbital.zig").buildNonRedundantPairs;
const buildOrbitalHessian = @import("casscf_orbital.zig").buildOrbitalHessian;
const buildStateAveragedActRdm = @import("reduced_density_matrix.zig").buildStateAveragedActRdm;
const getrf = @import("linear_algebra.zig").getrf;
const getrs = @import("linear_algebra.zig").getrs;
const slater = @import("slater_condon.zig").slater;

/// Computes analytic nuclear gradients of CASSCF state energies using dual number differentiation and CP-MCSCF response.
pub fn calculateCasscfGradient(comptime T: type, hfres: HartreeFockResult(T), C: Matrix(T), C_ci: Matrix(T), cas: CasState(T), states: []const u32, gpa: Allocator) ![]Matrix(T) {
    const dets, const weights, const space = .{ cas.dets, cas.weights, cas.space };

    var grads = try gpa.alloc(Matrix(T), states.len);
    errdefer gpa.free(grads);

    var grads_done: usize = 0;
    errdefer for (0..grads_done) |k| grads[k].deinit(gpa);

    for (0..states.len) |idx| {
        grads[idx] = try hfres.ints.sys.nuclearRepulsionGradient(gpa);

        grads_done += 1;
    }

    const natoms, const nocc = .{ hfres.ints.sys.atoms.len, space.nocc() };

    const nbf, const nmo, const nsp = .{ C.nrow(), space.nmo, 2 * nocc };

    const pairs = try buildNonRedundantPairs(space, gpa);
    defer gpa.free(pairs);

    const nrot = pairs.len;

    var D_act_sa, var d_act_sa = try buildStateAveragedActRdm(T, dets, C_ci, weights, space.ncore, space.nact, gpa);

    defer {
        D_act_sa.deinit(gpa);
        d_act_sa.deinit(gpa);
    }

    var Z_mat = try gpa.alloc(Matrix(T), states.len);

    defer {
        for (0..states.len) |k| Z_mat[k].deinit(gpa);

        gpa.free(Z_mat);
    }

    const is_ri, var naux: usize = .{ hfres.ints.B != null, 0 };

    if (is_ri) {
        naux = hfres.ints.B.?.shape[2];
    }

    var h_mo_conv: Matrix(T) = undefined;

    var B_mo_conv: Tensor(T, 3) = undefined;
    var g_mo_conv: Tensor(T, 4) = undefined;

    if (weights.len > 1 and nrot > 0) {
        h_mo_conv = try Matrix(T).init(nmo, nmo, gpa);

        if (is_ri) {
            B_mo_conv = try Tensor(T, 3).initZero(.{ nmo, nmo, naux }, gpa);
        }

        if (!is_ri) {
            g_mo_conv = try Tensor(T, 4).init(.{ nmo, nmo, nmo, nmo }, gpa);
        }
    }

    defer if (weights.len > 1 and nrot > 0) {
        h_mo_conv.deinit(gpa);

        if (is_ri) B_mo_conv.deinit(gpa) else g_mo_conv.deinit(gpa);
    };

    if (weights.len > 1 and nrot > 0) {
        try ao2mo_1e(T, &h_mo_conv, hfres.ints.H.?, C, gpa);

        var A_sa = try Matrix(T).initZero(nrot, nrot, gpa);
        defer A_sa.deinit(gpa);

        if (is_ri) {
            try ao2mo_3c(T, &B_mo_conv, hfres.ints.B.?, C, gpa);

            try buildOrbitalHessian(T, &A_sa, h_mo_conv, B_mo_conv, .{ D_act_sa, d_act_sa }, pairs, space, gpa);
        }

        if (!is_ri) {
            try ao2mo_2e(T, &g_mo_conv, hfres.ints.g.?, C, gpa);

            try buildOrbitalHessian(T, &A_sa, h_mo_conv, g_mo_conv, .{ D_act_sa, d_act_sa }, pairs, space, gpa);
        }

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

            var D_s, var d_s = try buildActRdm(T, dets, C_ci, s, space.ncore, space.nact, gpa);

            defer {
                D_s.deinit(gpa);
                d_s.deinit(gpa);
            }

            if (is_ri) {
                try buildGeneralizedFock(T, &F_gen_s, h_mo_conv, B_mo_conv, D_s, d_s, space, gpa);
            }

            if (!is_ri) {
                try buildGeneralizedFock(T, &F_gen_s, h_mo_conv, g_mo_conv, D_s, d_s, space, gpa);
            }

            var neg_g = try Matrix(T).init(nrot, 1, gpa);
            defer neg_g.deinit(gpa);

            for (pairs, 0..) |pair, i| {
                neg_g.ptr(i, 0).* = -2 * (F_gen_s.at(pair.p, pair.q) - F_gen_s.at(pair.q, pair.p));
            }

            try getrs(T, &Z_mat[s_idx], A_LU, ipiv, neg_g);
        }
    }

    if (weights.len <= 1 or nrot == 0) {
        for (states, 0..) |_, s_idx| {
            Z_mat[s_idx] = try Matrix(T).initZero(nrot, 1, gpa);
        }
    }

    var deriv_ws: ?FittingDerivativeWorkspace(T) = null;

    if (is_ri) {
        deriv_ws = try FittingDerivativeWorkspace(T).init(hfres.ints.J.?, hfres.ints.B.?.shape, gpa);
    }

    defer if (deriv_ws) |*w| w.deinit(gpa);

    var fock_ws: ?GeneralizedFockRiWorkspace(T) = null;

    if (is_ri and weights.len > 1 and nrot > 0) {
        fock_ws = try GeneralizedFockRiWorkspace(T).init(B_mo_conv, D_act_sa, d_act_sa, space, gpa);
    }

    defer if (fock_ws) |*w| w.deinit(gpa);

    var unrel_ws = try NonRelaxedOrbitalWorkspace(T).init(nbf, nmo, gpa);
    defer unrel_ws.deinit(gpa);

    var dC = try Matrix(T).init(nbf, nmo, gpa);
    defer dC.deinit(gpa);

    var h_mo_d = try DualMatrix(T).init(nmo, nmo, gpa);
    defer h_mo_d.deinit(gpa);

    var B_mo_d: DualTensor(T, 3) = undefined;
    var g_mo_d: DualTensor(T, 4) = undefined;
    var g_ms_d: DualTensor(T, 4) = undefined;

    if (is_ri) {
        B_mo_d = try DualTensor(T, 3).initZero(.{ nmo, nmo, naux }, gpa);

        g_ms_d = try DualTensor(T, 4).initZero(.{ nsp, nsp, nsp, nsp }, gpa);
    }

    if (!is_ri) {
        g_mo_d = try DualTensor(T, 4).init(.{ nmo, nmo, nmo, nmo }, gpa);
    }

    defer {
        if (is_ri) {
            B_mo_d.deinit(gpa);
            g_ms_d.deinit(gpa);
        }

        if (!is_ri) {
            g_mo_d.deinit(gpa);
        }
    }

    var h_ms = try Matrix(T).initZero(nsp, nsp, gpa);
    defer h_ms.deinit(gpa);

    var g_ms = try Tensor(T, 4).initZero(.{ nsp, nsp, nsp, nsp }, gpa);
    defer g_ms.deinit(gpa);

    var dF_gen: Matrix(T) = undefined;

    if (weights.len > 1 and nrot > 0) {
        dF_gen = try Matrix(T).initZero(nmo, nmo, gpa);
    }

    defer if (weights.len > 1 and nrot > 0) dF_gen.deinit(gpa);

    for (0..3 * natoms) |coord| {
        unrel_ws.evaluate(&dC, C, hfres.ints.dS.?.slice(coord).asMatrix());

        const C_d = DualMatrix(T).fromMatrices(C, dC);

        const H_d = DualMatrix(T).fromMatrices(hfres.ints.H.?, hfres.ints.dH.?.slice(coord).asMatrix());

        try ao2mo_1e(T, &h_mo_d, H_d, C_d, gpa);

        h_ms.zero();

        for (0..nocc) |p| for (0..nocc) |q| {
            h_ms.ptr(2 * p + 0, 2 * q + 0).* = h_mo_d.der.at(p, q);
            h_ms.ptr(2 * p + 1, 2 * q + 1).* = h_mo_d.der.at(p, q);
        };

        if (is_ri) {
            const B_d = try deriv_ws.?.evaluateCoord(hfres.ints, coord, gpa);

            try ao2mo_3c(T, &B_mo_d, B_d, C_d, gpa);

            try b2ms(T, &g_ms_d, B_mo_d, gpa);
        }

        if (!is_ri) {
            const g_d = DualTensor(T, 4).fromTensors(hfres.ints.g.?, hfres.ints.dg.?.slice(coord));

            try ao2mo_2e(T, &g_mo_d, g_d, C_d, gpa);

            g_ms.zero();

            for (0..nocc) |p| for (0..nocc) |r| for (0..nocc) |q| for (0..nocc) |s| {
                const val = g_mo_d.der.at(.{ p, q, r, s });

                g_ms.ptr(.{ 2 * p + 0, 2 * q + 0, 2 * r + 0, 2 * s + 0 }).* = val;
                g_ms.ptr(.{ 2 * p + 0, 2 * q + 1, 2 * r + 0, 2 * s + 1 }).* = val;
                g_ms.ptr(.{ 2 * p + 1, 2 * q + 0, 2 * r + 1, 2 * s + 0 }).* = val;
                g_ms.ptr(.{ 2 * p + 1, 2 * q + 1, 2 * r + 1, 2 * s + 1 }).* = val;
            };
        }

        var g_ms_der: Tensor(T, 4) = undefined;

        if (is_ri) {
            g_ms_der = g_ms_d.der;
        }

        if (!is_ri) {
            g_ms_der = g_ms;
        }

        for (0..dets.len) |m| for (m..dets.len) |n| {
            const factor: T = if (m == n) 1 else 2;

            const der = slater(T, dets[m], dets[n], h_ms, g_ms_der);

            for (states, 0..) |s, s_idx| {
                grads[s_idx].ptr(coord / 3, coord % 3).* += factor * C_ci.at(m, s) * C_ci.at(n, s) * der;
            }
        };

        if (weights.len > 1 and nrot > 0) {
            if (is_ri) {
                try buildGeneralizedFockDerivRi(T, &dF_gen, h_mo_d.der, .{ .B = B_mo_conv, .dB = B_mo_d.der }, .{ D_act_sa, d_act_sa }, space, fock_ws.?, gpa);
            }

            if (!is_ri) {
                try buildGeneralizedFock(T, &dF_gen, h_mo_d.der, g_mo_d.der, D_act_sa, d_act_sa, space, gpa);
            }

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
