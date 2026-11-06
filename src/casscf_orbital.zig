//! Implements CASSCF orbital optimization: generalized Fock matrices, orbital gradient, Hessian, and Cayley rotation.

const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

const ao2mo_1e = @import("integral_transform.zig").ao2mo_1e;
const ao2mo_2e = @import("integral_transform.zig").ao2mo_2e;
const gemm = @import("linear_algebra.zig").gemm;
const getrf = @import("linear_algebra.zig").getrf;
const getrs = @import("linear_algebra.zig").getrs;
const syevd = @import("linear_algebra.zig").syevd;

/// Orbital index pair representing a non-redundant rotation parameter between subspaces.
pub const OrbitalPair = struct {
    p: usize,
    q: usize,
};

/// Applies a unitary Cayley transform to update molecular orbital coefficients from rotation generator kappa.
pub fn applyCayleyTransform(comptime T: type, C: *Matrix(T), kappa: Matrix(T), gpa: Allocator) !void {
    const nmo = C.ncol();
    const nbf = C.nrow();

    var A = try Matrix(T).initZero(nmo, nmo, gpa);
    defer A.deinit(gpa);

    var U = try Matrix(T).initZero(nmo, nmo, gpa);
    defer U.deinit(gpa);

    for (0..nmo) |p| for (0..nmo) |q| {
        const delta: T = if (p == q) 1 else 0;

        A.ptr(p, q).* = delta - 0.5 * kappa.at(p, q);
        U.ptr(p, q).* = delta + 0.5 * kappa.at(p, q);
    };

    const ipiv = try gpa.alloc(i32, nmo);
    defer gpa.free(ipiv);

    try getrf(T, &A, ipiv);

    try getrs(T, &U, A, ipiv, U);

    var C_new = try Matrix(T).initZero(nbf, nmo, gpa);
    defer C_new.deinit(gpa);

    gemm(T, &C_new, C.*, U, 1, 0, false, false);

    @memcpy(C.data, C_new.data);
}

/// Constructs the core Fock matrix by contracting one-electron integrals with closed-shell core orbitals.
pub fn buildCoreFock(comptime T: type, F_core: *Matrix(T), h_mo: Matrix(T), g_mo: Tensor(T, 4), ncore: usize) void {
    const nmo = h_mo.nrow();

    for (0..nmo) |p| for (0..nmo) |q| {
        var sum = h_mo.at(p, q);

        for (0..ncore) |i| {
            sum += 2 * g_mo.at(.{ p, i, q, i }) - g_mo.at(.{ p, i, i, q });
        }

        F_core.ptr(p, q).* = sum;
    };
}

/// Constructs diagonal Fock matrix elements for orbital rotation preconditioner.
pub fn buildDiagonalFock(comptime T: type, F_diag: *Vector(T), h_mo: Matrix(T), g_mo: Tensor(T, 4), D_act: Matrix(T), ncore: usize, gpa: Allocator) !void {
    const nmo, const nact = .{ h_mo.nrow(), D_act.nrow() };

    var F_core = try Matrix(T).init(nmo, nmo, gpa);
    defer F_core.deinit(gpa);

    buildCoreFock(T, &F_core, h_mo, g_mo, ncore);

    for (0..nmo) |p| {
        var sum_diag = F_core.at(p, p);

        for (0..nact) |u| for (0..nact) |v| {
            const U, const V = .{ ncore + u, ncore + v };

            sum_diag += D_act.at(u, v) * (g_mo.at(.{ p, U, p, V }) - 0.5 * g_mo.at(.{ p, U, V, p }));
        };

        F_diag.ptr(p).* = sum_diag;
    }
}

/// Constructs the generalized Fock matrix using active reduced density matrices.
pub fn buildGeneralizedFock(comptime T: type, F_gen: *Matrix(T), h_mo: Matrix(T), g_mo: Tensor(T, 4), D_act: Matrix(T), d_act: Tensor(T, 4), ncore: usize, gpa: Allocator) !void {
    const nmo = h_mo.nrow();

    var F_core = try Matrix(T).init(nmo, nmo, gpa);
    defer F_core.deinit(gpa);

    buildCoreFock(T, &F_core, h_mo, g_mo, ncore);

    @memset(F_gen.data, 0);

    buildCoreColumns(T, F_gen, F_core, g_mo, D_act, ncore);

    buildActColumns(T, F_gen, F_core, g_mo, D_act, d_act, ncore);
}

/// Enumerates non-redundant orbital pairs across core-virtual, active-virtual, and core-active subspaces.
pub fn buildNonRedundantPairs(ncore: usize, nact: usize, nmo: usize, gpa: Allocator) ![]OrbitalPair {
    const total = ncore * (nmo - ncore - nact) + nact * (nmo - ncore - nact) + ncore * nact;

    var pairs, var idx: usize = .{ try gpa.alloc(OrbitalPair, total), 0 };

    for (0..ncore) |j| for (ncore + nact..nmo) |b| {
        pairs[idx] = .{ .p = b, .q = j };

        idx += 1;
    };

    for (ncore..ncore + nact) |u| for (ncore + nact..nmo) |b| {
        pairs[idx] = .{ .p = b, .q = u };

        idx += 1;
    };

    for (0..ncore) |j| for (ncore..ncore + nact) |u| {
        pairs[idx] = .{ .p = u, .q = j };

        idx += 1;
    };

    return pairs;
}

/// Evaluates the exact analytic orbital Hessian of the state-averaged energy with respect to orbital rotations.
pub fn buildOrbitalHessian(comptime T: type, A: *Matrix(T), h_mo: Matrix(T), g_mo: Tensor(T, 4), D_act: Matrix(T), d_act: Tensor(T, 4), pairs: []const OrbitalPair, ncore: usize, gpa: Allocator) !void {
    const nmo, const nrot = .{ h_mo.ncol(), pairs.len };

    var dh = try Matrix(T).initZero(nmo, nmo, gpa);
    defer dh.deinit(gpa);

    var dg = try Tensor(T, 4).initZero(.{ nmo, nmo, nmo, nmo }, gpa);
    defer dg.deinit(gpa);

    var dF_gen = try Matrix(T).initZero(nmo, nmo, gpa);
    defer dF_gen.deinit(gpa);

    for (pairs, 0..) |pair_j, j| {
        const p, const q = .{ pair_j.p, pair_j.q };

        @memset(dh.data, 0);

        for (0..nmo) |idx| {
            dh.ptr(p, idx).* += h_mo.at(q, idx);
            dh.ptr(q, idx).* -= h_mo.at(p, idx);
            dh.ptr(idx, p).* += h_mo.at(idx, q);
            dh.ptr(idx, q).* -= h_mo.at(idx, p);
        }

        @memset(dg.data, 0);

        for (0..nmo) |idx_a| for (0..nmo) |idx_b| for (0..nmo) |idx_c| {
            dg.ptr(.{ p, idx_a, idx_b, idx_c }).* += g_mo.at(.{ q, idx_a, idx_b, idx_c });
            dg.ptr(.{ q, idx_a, idx_b, idx_c }).* -= g_mo.at(.{ p, idx_a, idx_b, idx_c });

            dg.ptr(.{ idx_a, p, idx_b, idx_c }).* += g_mo.at(.{ idx_a, q, idx_b, idx_c });
            dg.ptr(.{ idx_a, q, idx_b, idx_c }).* -= g_mo.at(.{ idx_a, p, idx_b, idx_c });

            dg.ptr(.{ idx_a, idx_b, p, idx_c }).* += g_mo.at(.{ idx_a, idx_b, q, idx_c });
            dg.ptr(.{ idx_a, idx_b, q, idx_c }).* -= g_mo.at(.{ idx_a, idx_b, p, idx_c });

            dg.ptr(.{ idx_a, idx_b, idx_c, p }).* += g_mo.at(.{ idx_a, idx_b, idx_c, q });
            dg.ptr(.{ idx_a, idx_b, idx_c, q }).* -= g_mo.at(.{ idx_a, idx_b, idx_c, p });
        };

        @memset(dF_gen.data, 0);

        try buildGeneralizedFock(T, &dF_gen, dh, dg, D_act, d_act, ncore, gpa);

        for (pairs, 0..) |pair_i, i| {
            A.ptr(i, j).* = 2 * (dF_gen.at(pair_i.p, pair_i.q) - dF_gen.at(pair_i.q, pair_i.p));
        }
    }

    for (0..nrot) |i| for (i + 1..nrot) |j| {
        const avg = 0.5 * (A.at(i, j) + A.at(j, i));

        A.ptr(i, j).* = avg;
        A.ptr(j, i).* = avg;
    };
}

/// Evaluates the anti-symmetric orbital rotation generator matrix kappa using a diagonal Hessian approximation.
pub fn buildOrbitalRotationGenerator(comptime T: type, kappa: *Matrix(T), g_orb: Matrix(T), F_diag: Vector(T), D_act: Matrix(T), ncore: usize) void {
    const nmo, const nact = .{ kappa.ncol(), D_act.nrow() };

    for (0..ncore) |j| for (ncore + nact..nmo) |b| {
        const denom = @max(0.1, @abs(4 * (F_diag.at(b) - F_diag.at(j))));

        const step = std.math.clamp(-g_orb.at(b, j) / denom, -0.2, 0.2);

        kappa.ptr(b, j).* = 0 + step;
        kappa.ptr(j, b).* = 0 - step;
    };

    for (0..nact) |u_idx| {
        const u = ncore + u_idx;

        const Duu = D_act.at(u_idx, u_idx);

        for (ncore + nact..nmo) |b| {
            const denom = @max(0.1, @abs(2 * Duu * (F_diag.at(b) - F_diag.at(u))));

            const step = std.math.clamp(-g_orb.at(b, u) / denom, -0.2, 0.2);

            kappa.ptr(b, u).* = 0 + step;
            kappa.ptr(u, b).* = 0 - step;
        }

        for (0..ncore) |j| {
            const denom = @max(0.1, @abs(2 * (2 - Duu) * (F_diag.at(u) - F_diag.at(j))));

            const step = std.math.clamp(-g_orb.at(u, j) / denom, -0.2, 0.2);

            kappa.ptr(u, j).* = 0 + step;
            kappa.ptr(j, u).* = 0 - step;
        }
    }
}

/// Canonicalizes core and virtual orbital blocks by diagonalizing their respective Fock matrices.
pub fn canonicalizeOrbitals(comptime T: type, C: *Matrix(T), ints: anytype, ncore: usize, nact: usize, gpa: Allocator) !void {
    const nmo = C.ncol();

    var h_mo = try Matrix(T).init(nmo, nmo, gpa);
    defer h_mo.deinit(gpa);

    ao2mo_1e(T, &h_mo, ints.H.?, C.*);

    var g_mo = try Tensor(T, 4).init(.{ nmo, nmo, nmo, nmo }, gpa);
    defer g_mo.deinit(gpa);

    try ao2mo_2e(T, &g_mo, ints.g.?, C.*, gpa);

    var F_core = try Matrix(T).init(nmo, nmo, gpa);
    defer F_core.deinit(gpa);

    buildCoreFock(T, &F_core, h_mo, g_mo, ncore);

    try canonicalizeCoreSubspace(T, C, F_core, ncore, gpa);

    try canonicalizeVirtualSubspace(T, C, F_core, ncore, nact, gpa);
}

/// Evaluates the orbital gradient matrix and returns the maximum gradient norm among non-redundant pairs.
pub fn computeOrbitalGradient(comptime T: type, g_orb: *Matrix(T), F_gen: Matrix(T), ncore: usize, nact: usize) T {
    const nmo, var max_norm: T = .{ F_gen.nrow(), 0 };

    for (0..nmo) |p| for (0..nmo) |q| {
        g_orb.ptr(p, q).* = 2 * (F_gen.at(p, q) - F_gen.at(q, p));
    };

    for (0..ncore) |j| for (ncore + nact..nmo) |b| {
        max_norm = @max(max_norm, @abs(g_orb.at(b, j)));
    };

    for (ncore..ncore + nact) |u| for (ncore + nact..nmo) |b| {
        max_norm = @max(max_norm, @abs(g_orb.at(b, u)));
    };

    for (0..ncore) |j| for (ncore..ncore + nact) |u| {
        max_norm = @max(max_norm, @abs(g_orb.at(u, j)));
    };

    return max_norm;
}

/// Constructs active subspace columns of the generalized Fock matrix from active 1- and 2-RDM contractions.
fn buildActColumns(comptime T: type, F_gen: *Matrix(T), F_core: Matrix(T), g_mo: Tensor(T, 4), D_act: Matrix(T), d_act: Tensor(T, 4), ncore: usize) void {
    const nmo, const nact = .{ F_core.nrow(), D_act.nrow() };

    for (0..nmo) |p| for (0..nact) |u| {
        const U = ncore + u;
        var sum: T = 0;

        for (0..nact) |v| {
            sum += F_core.at(p, ncore + v) * D_act.at(v, u);
        }

        for (0..nact) |v| for (0..nact) |w| for (0..nact) |x| {
            const V, const W, const X = .{ ncore + v, ncore + w, ncore + x };

            sum += g_mo.at(.{ p, W, V, X }) * d_act.at(.{ u, v, w, x });
        };

        F_gen.ptr(p, U).* = sum;
    };
}

/// Constructs core subspace columns of the generalized Fock matrix from active 1-RDM contractions.
fn buildCoreColumns(comptime T: type, F_gen: *Matrix(T), F_core: Matrix(T), g_mo: Tensor(T, 4), D_act: Matrix(T), ncore: usize) void {
    const nmo, const nact = .{ F_core.nrow(), D_act.nrow() };

    for (0..nmo) |p| for (0..ncore) |j| {
        var sum = F_core.at(p, j);

        for (0..nact) |u| for (0..nact) |v| {
            const U, const V = .{ ncore + u, ncore + v };

            sum += D_act.at(u, v) * (g_mo.at(.{ p, U, j, V }) - 0.5 * g_mo.at(.{ p, U, V, j }));
        };

        F_gen.ptr(p, j).* = 2 * sum;
    };
}

/// Canonicalizes the doubly occupied core orbital subspace by diagonalizing its core Fock sub-block.
fn canonicalizeCoreSubspace(comptime T: type, C: *Matrix(T), F_core: Matrix(T), ncore: usize, gpa: Allocator) !void {
    if (ncore <= 1) return;

    const nbf = C.nrow();

    var F_cc = try Matrix(T).init(ncore, ncore, gpa);
    defer F_cc.deinit(gpa);

    for (0..ncore) |i| for (0..ncore) |j| {
        F_cc.ptr(i, j).* = F_core.at(i, j);
    };

    var E_cc = try Vector(T).init(ncore, gpa);
    defer E_cc.deinit(gpa);

    var V_cc = try Matrix(T).initZero(ncore, ncore, gpa);
    defer V_cc.deinit(gpa);

    try syevd(T, &E_cc, &V_cc, F_cc);

    var C_core = try Matrix(T).init(nbf, ncore, gpa);
    defer C_core.deinit(gpa);

    for (0..nbf) |mu| for (0..ncore) |i| {
        C_core.ptr(mu, i).* = C.at(mu, i);
    };

    var C_new = try Matrix(T).initZero(nbf, ncore, gpa);
    defer C_new.deinit(gpa);

    gemm(T, &C_new, C_core, V_cc, 1.0, 0.0, false, false);

    for (0..nbf) |mu| for (0..ncore) |i| {
        C.ptr(mu, i).* = C_new.at(mu, i);
    };
}

/// Canonicalizes the unoccupied virtual orbital subspace by diagonalizing its core Fock sub-block.
fn canonicalizeVirtualSubspace(comptime T: type, C: *Matrix(T), F_core: Matrix(T), ncore: usize, nact: usize, gpa: Allocator) !void {
    const nmo = C.ncol();

    const nvirt = nmo - ncore - nact;

    if (nvirt <= 1) return;

    const nbf, const v_start = .{ C.nrow(), ncore + nact };

    var F_vv = try Matrix(T).init(nvirt, nvirt, gpa);
    defer F_vv.deinit(gpa);

    for (0..nvirt) |a| for (0..nvirt) |b| {
        F_vv.ptr(a, b).* = F_core.at(v_start + a, v_start + b);
    };

    var E_vv = try Vector(T).init(nvirt, gpa);
    defer E_vv.deinit(gpa);

    var V_vv = try Matrix(T).initZero(nvirt, nvirt, gpa);
    defer V_vv.deinit(gpa);

    try syevd(T, &E_vv, &V_vv, F_vv);

    var C_virt = try Matrix(T).init(nbf, nvirt, gpa);
    defer C_virt.deinit(gpa);

    for (0..nbf) |mu| for (0..nvirt) |a| {
        C_virt.ptr(mu, a).* = C.at(mu, v_start + a);
    };

    var C_new = try Matrix(T).initZero(nbf, nvirt, gpa);
    defer C_new.deinit(gpa);

    gemm(T, &C_new, C_virt, V_vv, 1, 0, false, false);

    for (0..nbf) |mu| for (0..nvirt) |a| {
        C.ptr(mu, v_start + a).* = C_new.at(mu, a);
    };
}
