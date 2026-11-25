//! Implements CASSCF orbital optimization: generalized Fock matrices, orbital gradient, Hessian, and Cayley rotation.

const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;
const Value = @import("value.zig").Value;
const Vector = @import("tensor.zig").Vector;

const ao2mo_1e = @import("integral_transform.zig").ao2mo_1e;
const ao2mo_2e = @import("integral_transform.zig").ao2mo_2e;
const ao2mo_3c = @import("integral_transform.zig").ao2mo_3c;
const gemm = @import("linear_algebra.zig").gemm;
const getrf = @import("linear_algebra.zig").getrf;
const getrs = @import("linear_algebra.zig").getrs;
const syevd = @import("linear_algebra.zig").syevd;

/// Defines orbital subspace dimensions partitioned into core, active, and molecular orbital counts.
pub const CasSpace = struct {
    ncore: usize,

    nact: usize,

    nmo: usize,

    /// Returns the total number of occupied spatial orbitals including core and active subspaces.
    pub fn nocc(self: @This()) usize {
        return self.ncore + self.nact;
    }

    /// Returns the total number of unoccupied virtual spatial orbitals.
    pub fn nvir(self: @This()) usize {
        return self.nmo - self.ncore - self.nact;
    }
};

/// Orbital index pair representing a non-redundant rotation parameter between subspaces.
pub const OrbitalPair = struct {
    p: usize,
    q: usize,
};

/// Holds multiconfiguration determinant basis, state weights, and active space partitioning.
pub fn CasState(comptime T: type) type {
    return struct {
        dets: []const []const usize,
        weights: []const T,
        space: CasSpace,
    };
}

/// Holds 3-center density fitting integrals and their nuclear derivative tensor.
pub fn FittingDeriv(comptime T: type) type {
    return struct { B: Tensor(T, 3), dB: Tensor(T, 3) };
}

/// Intermediate density fitting contractions with active density matrices for generalized Fock derivative evaluations.
pub fn GeneralizedFockRiWorkspace(comptime T: type) type {
    return struct {
        v_core: Vector(T),

        W: Tensor(T, 3),
        Q: Tensor(T, 3),

        v_act: Vector(T),

        /// Precomputes and caches 3-center contractions with active 1- and 2-RDMs for derivative evaluations.
        pub fn init(B_mo: Tensor(T, 3), D_act: Matrix(T), d_act: Tensor(T, 4), space: CasSpace, gpa: Allocator) !@This() {
            @setFloatMode(.optimized);

            const nact, const naux = .{ space.nact, B_mo.shape[2] };

            var v_core = try Vector(T).initZero(naux, gpa);
            errdefer v_core.deinit(gpa);

            for (0..space.ncore) |i| for (0..naux) |P| {
                v_core.ptr(P).* += 2 * B_mo.at(.{ i, i, P });
            };

            var v_act = try Vector(T).initZero(naux, gpa);
            errdefer v_act.deinit(gpa);

            for (0..nact) |u| for (0..nact) |v| {
                const val_D = D_act.at(u, v);

                if (val_D == 0) continue;

                const U, const V = .{ space.ncore + u, space.ncore + v };

                for (0..naux) |P| {
                    v_act.ptr(P).* += val_D * B_mo.at(.{ U, V, P });
                }
            };

            var W = try Tensor(T, 3).initZero(.{ nact, space.ncore, naux }, gpa);
            errdefer W.deinit(gpa);

            for (0..nact) |u| for (0..nact) |v| {
                const val_D = D_act.at(u, v);

                if (val_D == 0) continue;

                const U = space.ncore + u;

                for (0..space.ncore) |j| for (0..naux) |P| {
                    W.ptr(.{ v, j, P }).* += val_D * B_mo.at(.{ U, j, P });
                };
            };

            var Q = try Tensor(T, 3).initZero(.{ nact, nact, naux }, gpa);
            errdefer Q.deinit(gpa);

            for (0..nact) |u| for (0..nact) |v| for (0..nact) |w| for (0..nact) |x| {
                const val_d = d_act.at(.{ u, v, w, x });

                if (val_d == 0) continue;

                const W_idx, const X_idx = .{ space.ncore + w, space.ncore + x };

                for (0..naux) |P| {
                    Q.ptr(.{ u, v, P }).* += val_d * B_mo.at(.{ W_idx, X_idx, P });
                }
            };

            return .{ .v_core = v_core, .v_act = v_act, .W = W, .Q = Q };
        }

        /// Releases cached density fitting intermediate tensors.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.v_core.deinit(gpa);

            self.W.deinit(gpa);
            self.Q.deinit(gpa);

            self.v_act.deinit(gpa);
        }
    };
}

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

/// Constructs core Fock matrix in molecular orbital basis from 4-center or 3-center integrals.
pub fn buildCoreFock(comptime T: type, F_core: *Matrix(T), h_mo: Matrix(T), two_e: anytype, space: CasSpace, gpa: Allocator) !void {
    const nmo = space.nmo;

    if (two_e.shape.len == 4) {
        for (0..nmo) |p| for (0..nmo) |q| {
            var sum = h_mo.at(p, q);

            for (0..space.ncore) |i| {
                sum += 2 * two_e.at(.{ p, i, q, i }) - two_e.at(.{ p, i, i, q });
            }

            F_core.ptr(p, q).* = sum;
        };

        return;
    }

    const naux = two_e.shape[2];

    var v_core = try Vector(T).initZero(naux, gpa);
    defer v_core.deinit(gpa);

    for (0..space.ncore) |i| for (0..naux) |P| {
        v_core.ptr(P).* += 2 * two_e.at(.{ i, i, P });
    };

    for (0..nmo) |p| for (0..nmo) |q| {
        var sum = h_mo.at(p, q);

        for (0..naux) |P| {
            sum += two_e.at(.{ p, q, P }) * v_core.at(P);
        }

        for (0..space.ncore) |i| for (0..naux) |P| {
            sum -= two_e.at(.{ p, i, P }) * two_e.at(.{ q, i, P });
        };

        F_core.ptr(p, q).* = sum;
    };
}

/// Constructs diagonal Fock matrix elements for orbital rotation preconditioner from 4-center or 3-center integrals.
pub fn buildDiagonalFock(comptime T: type, F_diag: *Vector(T), F_core: Matrix(T), two_e: anytype, D_act: Matrix(T), space: CasSpace, gpa: Allocator) !void {
    const nmo, const nact = .{ space.nmo, space.nact };

    if (two_e.shape.len == 4) {
        for (0..nmo) |p| {
            var sum_diag = F_core.at(p, p);

            for (0..nact) |u| for (0..nact) |v| {
                const U, const V = .{ space.ncore + u, space.ncore + v };

                sum_diag += D_act.at(u, v) * (two_e.at(.{ p, U, p, V }) - 0.5 * two_e.at(.{ p, U, V, p }));
            };

            F_diag.ptr(p).* = sum_diag;
        }

        return;
    }

    const naux = two_e.shape[2];

    var v_act = try Vector(T).initZero(naux, gpa);
    defer v_act.deinit(gpa);

    for (0..nact) |u| for (0..nact) |v| {
        const val_D = D_act.at(u, v);

        if (val_D == 0) continue;

        const U, const V = .{ space.ncore + u, space.ncore + v };

        for (0..naux) |P| {
            v_act.ptr(P).* += val_D * two_e.at(.{ U, V, P });
        }
    };

    for (0..nmo) |p| {
        var sum_diag = F_core.at(p, p);

        for (0..naux) |P| {
            sum_diag += two_e.at(.{ p, p, P }) * v_act.at(P);
        }

        for (0..nact) |u| for (0..nact) |v| {
            const val_D = D_act.at(u, v);

            if (val_D == 0) continue;

            const U, const V = .{ space.ncore + u, space.ncore + v };

            for (0..naux) |P| {
                sum_diag -= 0.5 * val_D * two_e.at(.{ p, U, P }) * two_e.at(.{ p, V, P });
            }
        };

        F_diag.ptr(p).* = sum_diag;
    }
}

/// Constructs generalized Fock matrix from active RDMs and 4-center or 3-center integrals.
pub fn buildGeneralizedFock(comptime T: type, F_gen: *Matrix(T), h_mo: Matrix(T), two_e: anytype, D_act: Matrix(T), d_act: Tensor(T, 4), space: CasSpace, gpa: Allocator) !void {
    const nmo = space.nmo;

    var F_core = try Matrix(T).init(nmo, nmo, gpa);
    defer F_core.deinit(gpa);

    try buildCoreFock(T, &F_core, h_mo, two_e, space, gpa);

    @memset(F_gen.data, 0);

    if (two_e.shape.len == 4) {
        buildCoreColumns(T, F_gen, F_core, two_e, D_act, space);

        buildActColumns(T, F_gen, F_core, two_e, D_act, d_act, space);
    }

    if (two_e.shape.len == 3) {
        try buildCoreColumnsRi(T, F_gen, F_core, two_e, D_act, space, gpa);

        try buildActColumnsRi(T, F_gen, F_core, two_e, D_act, d_act, space, gpa);
    }
}

/// Evaluates the generalized Fock matrix derivative using density fitting and precomputed intermediate contractions.
pub fn buildGeneralizedFockDerivRi(comptime T: type, dF: *Matrix(T), dh: Matrix(T), b: FittingDeriv(T), rdm: anytype, sp: CasSpace, ws: anytype, gpa: Allocator) !void {
    const nmo, const D_act = .{ sp.nmo, rdm[0] };

    var dF_core = try Matrix(T).init(nmo, nmo, gpa);
    defer dF_core.deinit(gpa);

    try buildCoreFockDerivRi(T, &dF_core, dh, b, ws.v_core, sp, gpa);

    @memset(dF.data, 0);

    try buildCoreColumnsDerivRi(T, dF, dF_core, b, D_act, ws, sp, gpa);

    try buildActColumnsDerivRi(T, dF, dF_core, b, rdm, ws.Q, sp, gpa);
}

/// Enumerates non-redundant orbital pairs across core-virtual, active-virtual, and core-active subspaces.
pub fn buildNonRedundantPairs(space: CasSpace, gpa: Allocator) ![]OrbitalPair {
    const total = space.ncore * space.nvir() + space.nact * space.nvir() + space.ncore * space.nact;

    var pairs, var idx: usize = .{ try gpa.alloc(OrbitalPair, total), 0 };

    for (0..space.ncore) |j| for (space.nocc()..space.nmo) |b| {
        pairs[idx] = .{ .p = b, .q = j };

        idx += 1;
    };

    for (space.ncore..space.nocc()) |u| for (space.nocc()..space.nmo) |b| {
        pairs[idx] = .{ .p = b, .q = u };

        idx += 1;
    };

    for (0..space.ncore) |j| for (space.ncore..space.nocc()) |u| {
        pairs[idx] = .{ .p = u, .q = j };

        idx += 1;
    };

    return pairs;
}

/// Evaluates the analytic orbital Hessian of the state-averaged energy with respect to orbital rotations.
pub fn buildOrbitalHessian(comptime T: type, A: *Matrix(T), h_mo: Matrix(T), two_e: anytype, rdm: anytype, pairs: []const OrbitalPair, space: CasSpace, gpa: Allocator) !void {
    const nmo = space.nmo;

    const D_act, const d_act = rdm;

    var dh = try Matrix(T).initZero(nmo, nmo, gpa);
    defer dh.deinit(gpa);

    var dF_gen = try Matrix(T).initZero(nmo, nmo, gpa);
    defer dF_gen.deinit(gpa);

    if (two_e.shape.len == 4) {
        var dg = try Tensor(T, 4).initZero(.{ nmo, nmo, nmo, nmo }, gpa);
        defer dg.deinit(gpa);

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
                dg.ptr(.{ p, idx_a, idx_b, idx_c }).* += two_e.at(.{ q, idx_a, idx_b, idx_c });
                dg.ptr(.{ q, idx_a, idx_b, idx_c }).* -= two_e.at(.{ p, idx_a, idx_b, idx_c });

                dg.ptr(.{ idx_a, p, idx_b, idx_c }).* += two_e.at(.{ idx_a, q, idx_b, idx_c });
                dg.ptr(.{ idx_a, q, idx_b, idx_c }).* -= two_e.at(.{ idx_a, p, idx_b, idx_c });

                dg.ptr(.{ idx_a, idx_b, p, idx_c }).* += two_e.at(.{ idx_a, idx_b, q, idx_c });
                dg.ptr(.{ idx_a, idx_b, q, idx_c }).* -= two_e.at(.{ idx_a, idx_b, p, idx_c });

                dg.ptr(.{ idx_a, idx_b, idx_c, p }).* += two_e.at(.{ idx_a, idx_b, idx_c, q });
                dg.ptr(.{ idx_a, idx_b, idx_c, q }).* -= two_e.at(.{ idx_a, idx_b, idx_c, p });
            };

            @memset(dF_gen.data, 0);

            try buildGeneralizedFock(T, &dF_gen, dh, dg, D_act, d_act, space, gpa);

            for (pairs, 0..) |pair_i, i| {
                A.ptr(i, j).* = 2 * (dF_gen.at(pair_i.p, pair_i.q) - dF_gen.at(pair_i.q, pair_i.p));
            }
        }
    }

    if (two_e.shape.len == 3) {
        const naux = two_e.shape[2];

        var ws = try GeneralizedFockRiWorkspace(T).init(two_e, D_act, d_act, space, gpa);
        defer ws.deinit(gpa);

        var dB = try Tensor(T, 3).initZero(.{ nmo, nmo, naux }, gpa);
        defer dB.deinit(gpa);

        for (pairs, 0..) |pair_j, j| {
            const p, const q = .{ pair_j.p, pair_j.q };

            @memset(dh.data, 0);

            for (0..nmo) |idx| {
                dh.ptr(p, idx).* += h_mo.at(q, idx);
                dh.ptr(q, idx).* -= h_mo.at(p, idx);
                dh.ptr(idx, p).* += h_mo.at(idx, q);
                dh.ptr(idx, q).* -= h_mo.at(idx, p);
            }

            @memset(dB.data, 0);

            for (0..nmo) |idx| for (0..naux) |P| {
                dB.ptr(.{ p, idx, P }).* += two_e.at(.{ q, idx, P });
                dB.ptr(.{ q, idx, P }).* -= two_e.at(.{ p, idx, P });
                dB.ptr(.{ idx, p, P }).* += two_e.at(.{ idx, q, P });
                dB.ptr(.{ idx, q, P }).* -= two_e.at(.{ idx, p, P });
            };

            try buildGeneralizedFockDerivRi(T, &dF_gen, dh, .{ .B = two_e, .dB = dB }, .{ D_act, d_act }, space, ws, gpa);

            for (pairs, 0..) |pair_i, i| {
                A.ptr(i, j).* = 2 * (dF_gen.at(pair_i.p, pair_i.q) - dF_gen.at(pair_i.q, pair_i.p));
            }
        }
    }

    A.symmetrize();
}

/// Evaluates the anti-symmetric orbital rotation generator matrix kappa using a diagonal Hessian approximation.
pub fn buildOrbitalRotationGenerator(comptime T: type, kappa: *Matrix(T), g_orb: Matrix(T), F_diag: Vector(T), D_act: Matrix(T), space: CasSpace) void {
    for (0..space.ncore) |j| for (space.nocc()..space.nmo) |b| {
        const denom = @max(0.1, @abs(4 * (F_diag.at(b) - F_diag.at(j))));

        const step = std.math.clamp(-g_orb.at(b, j) / denom, -0.2, 0.2);

        kappa.ptr(b, j).* = 0 + step;
        kappa.ptr(j, b).* = 0 - step;
    };

    for (0..space.nact) |u_idx| {
        const u = space.ncore + u_idx;

        const Duu = D_act.at(u_idx, u_idx);

        for (space.nocc()..space.nmo) |b| {
            const denom = @max(0.1, @abs(2 * Duu * (F_diag.at(b) - F_diag.at(u))));

            const step = std.math.clamp(-g_orb.at(b, u) / denom, -0.2, 0.2);

            kappa.ptr(b, u).* = 0 + step;
            kappa.ptr(u, b).* = 0 - step;
        }

        for (0..space.ncore) |j| {
            const denom = @max(0.1, @abs(2 * (2 - Duu) * (F_diag.at(u) - F_diag.at(j))));

            const step = std.math.clamp(-g_orb.at(u, j) / denom, -0.2, 0.2);

            kappa.ptr(u, j).* = 0 + step;
            kappa.ptr(j, u).* = 0 - step;
        }
    }
}

/// Canonicalizes core and virtual orbital blocks by diagonalizing their respective Fock matrices.
pub fn canonicalizeOrbitals(comptime T: type, C: *Matrix(T), ints: anytype, space: CasSpace, gpa: Allocator) !void {
    const nmo = space.nmo;

    var h_mo = try Matrix(T).init(nmo, nmo, gpa);
    defer h_mo.deinit(gpa);

    try ao2mo_1e(T, &h_mo, ints.H.?, C.*, gpa);

    var F_core = try Matrix(T).init(nmo, nmo, gpa);
    defer F_core.deinit(gpa);

    if (ints.B != null) {
        var B_mo = try Tensor(T, 3).initZero(.{ nmo, nmo, ints.B.?.shape[2] }, gpa);
        defer B_mo.deinit(gpa);

        try ao2mo_3c(T, &B_mo, ints.B.?, C.*, gpa);

        try buildCoreFock(T, &F_core, h_mo, B_mo, space, gpa);
    }

    if (ints.B == null) {
        var g_mo = try Tensor(T, 4).init(.{ nmo, nmo, nmo, nmo }, gpa);
        defer g_mo.deinit(gpa);

        try ao2mo_2e(T, &g_mo, ints.g.?, C.*, gpa);

        try buildCoreFock(T, &F_core, h_mo, g_mo, space, gpa);
    }

    try canonicalizeCoreSubspace(T, C, F_core, space.ncore, gpa);

    try canonicalizeVirtualSubspace(T, C, F_core, space, gpa);
}

/// Evaluates the orbital gradient matrix and returns the maximum gradient norm among non-redundant pairs.
pub fn computeOrbitalGradient(comptime T: type, g_orb: *Matrix(T), F_gen: Matrix(T), space: CasSpace) T {
    const nmo, var max_norm: T = .{ space.nmo, 0 };

    for (0..nmo) |p| for (0..nmo) |q| {
        g_orb.ptr(p, q).* = 2 * (F_gen.at(p, q) - F_gen.at(q, p));
    };

    for (0..space.ncore) |j| for (space.nocc()..space.nmo) |b| {
        max_norm = @max(max_norm, @abs(g_orb.at(b, j)));
    };

    for (space.ncore..space.nocc()) |u| for (space.nocc()..space.nmo) |b| {
        max_norm = @max(max_norm, @abs(g_orb.at(b, u)));
    };

    for (0..space.ncore) |j| for (space.ncore..space.nocc()) |u| {
        max_norm = @max(max_norm, @abs(g_orb.at(u, j)));
    };

    return max_norm;
}

/// Constructs active subspace columns of the generalized Fock matrix from active 1- and 2-RDM contractions.
fn buildActColumns(comptime T: type, F_gen: *Matrix(T), F_core: Matrix(T), g_mo: Tensor(T, 4), D_act: Matrix(T), d_act: Tensor(T, 4), space: CasSpace) void {
    const nmo, const nact = .{ space.nmo, space.nact };

    for (0..nmo) |p| for (0..nact) |u| {
        const U, var sum: T = .{ space.ncore + u, 0 };

        for (0..nact) |v| {
            sum += F_core.at(p, space.ncore + v) * D_act.at(v, u);
        }

        for (0..nact) |v| for (0..nact) |w| for (0..nact) |x| {
            const V, const W, const X = .{ space.ncore + v, space.ncore + w, space.ncore + x };

            sum += g_mo.at(.{ p, W, V, X }) * d_act.at(.{ u, v, w, x });
        };

        F_gen.ptr(p, U).* = sum;
    };
}

/// Constructs active subspace columns of generalized Fock matrix from active RDMs and fitting integrals.
fn buildActColumnsRi(comptime T: type, F_gen: *Matrix(T), F_core: Matrix(T), B_mo: Tensor(T, 3), D_act: Matrix(T), d_act: Tensor(T, 4), space: CasSpace, gpa: Allocator) !void {
    const nmo, const nact = .{ space.nmo, space.nact };

    const naux = B_mo.shape[2];

    var Q_aux = try Tensor(T, 3).initZero(.{ nact, nact, naux }, gpa);
    defer Q_aux.deinit(gpa);

    for (0..nact) |u| for (0..nact) |v| for (0..nact) |w| for (0..nact) |x| {
        const val_d = d_act.at(.{ u, v, w, x });

        if (val_d == 0) continue;

        const W, const X = .{ space.ncore + w, space.ncore + x };

        for (0..naux) |P| {
            Q_aux.ptr(.{ u, v, P }).* += val_d * B_mo.at(.{ W, X, P });
        }
    };

    for (0..nmo) |p| for (0..nact) |u| {
        const U, var sum: T = .{ space.ncore + u, 0 };

        for (0..nact) |v| {
            sum += F_core.at(p, space.ncore + v) * D_act.at(v, u);
        }

        for (0..nact) |v| for (0..naux) |P| {
            sum += B_mo.at(.{ p, space.ncore + v, P }) * Q_aux.at(.{ u, v, P });
        };

        F_gen.ptr(p, U).* = sum;
    };
}

/// Constructs core subspace columns of the generalized Fock matrix from active 1-RDM contractions.
fn buildCoreColumns(comptime T: type, F_gen: *Matrix(T), F_core: Matrix(T), g_mo: Tensor(T, 4), D_act: Matrix(T), space: CasSpace) void {
    const nmo, const nact = .{ space.nmo, space.nact };

    for (0..nmo) |p| for (0..space.ncore) |j| {
        var sum = F_core.at(p, j);

        for (0..nact) |u| for (0..nact) |v| {
            const U, const V = .{ space.ncore + u, space.ncore + v };

            sum += D_act.at(u, v) * (g_mo.at(.{ p, U, j, V }) - 0.5 * g_mo.at(.{ p, U, V, j }));
        };

        F_gen.ptr(p, j).* = 2 * sum;
    };
}

/// Constructs core subspace columns of generalized Fock matrix from active 1-RDM and fitting integrals.
fn buildCoreColumnsRi(comptime T: type, F_gen: *Matrix(T), F_core: Matrix(T), B_mo: Tensor(T, 3), D_act: Matrix(T), space: CasSpace, gpa: Allocator) !void {
    const nmo, const nact = .{ space.nmo, space.nact };

    const naux = B_mo.shape[2];

    var v_act = try Vector(T).initZero(naux, gpa);
    defer v_act.deinit(gpa);

    for (0..nact) |u| for (0..nact) |v| {
        const val_D = D_act.at(u, v);

        if (val_D == 0) continue;

        const U, const V = .{ space.ncore + u, space.ncore + v };

        for (0..naux) |P| {
            v_act.ptr(P).* += val_D * B_mo.at(.{ U, V, P });
        }
    };

    var W = try Tensor(T, 3).initZero(.{ nact, space.ncore, naux }, gpa);
    defer W.deinit(gpa);

    for (0..nact) |u| for (0..nact) |v| {
        const val_D = D_act.at(u, v);

        if (val_D == 0) continue;

        const U = space.ncore + u;

        for (0..space.ncore) |j| for (0..naux) |P| {
            W.ptr(.{ v, j, P }).* += val_D * B_mo.at(.{ U, j, P });
        };
    };

    for (0..nmo) |p| for (0..space.ncore) |j| {
        var sum = F_core.at(p, j);

        for (0..naux) |P| {
            sum += B_mo.at(.{ p, j, P }) * v_act.at(P);
        }

        for (0..nact) |v| for (0..naux) |P| {
            sum -= 0.5 * B_mo.at(.{ p, space.ncore + v, P }) * W.at(.{ v, j, P });
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
fn canonicalizeVirtualSubspace(comptime T: type, C: *Matrix(T), F_core: Matrix(T), space: CasSpace, gpa: Allocator) !void {
    const nvirt = space.nvir();

    if (nvirt <= 1) return;

    const nbf, const v_start = .{ C.nrow(), space.nocc() };

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

/// Evaluates derivative of core Fock matrix with respect to orbital rotation or nuclear displacement.
fn buildCoreFockDerivRi(comptime T: type, dF_core: *Matrix(T), dh: Matrix(T), b: FittingDeriv(T), v_core: Vector(T), space: CasSpace, gpa: Allocator) !void {
    const B_mo, const dB_mo = .{ b.B, b.dB };

    const nmo, const naux = .{ space.nmo, B_mo.shape[2] };

    var dv_core = try Vector(T).initZero(naux, gpa);
    defer dv_core.deinit(gpa);

    for (0..space.ncore) |i| for (0..naux) |P| {
        dv_core.ptr(P).* += 2 * dB_mo.at(.{ i, i, P });
    };

    for (0..nmo) |p| for (0..nmo) |q| {
        var sum = dh.at(p, q);

        for (0..naux) |P| {
            sum += dB_mo.at(.{ p, q, P }) * v_core.at(P) + B_mo.at(.{ p, q, P }) * dv_core.at(P);
        }

        for (0..space.ncore) |i| for (0..naux) |P| {
            sum -= dB_mo.at(.{ p, i, P }) * B_mo.at(.{ q, i, P }) + B_mo.at(.{ p, i, P }) * dB_mo.at(.{ q, i, P });
        };

        dF_core.ptr(p, q).* = sum;
    };
}

/// Constructs derivative of core subspace columns of generalized Fock matrix from 3-center fitting derivatives.
fn buildCoreColumnsDerivRi(comptime T: type, dF: *Matrix(T), dF_core: Matrix(T), b: FittingDeriv(T), D_act: Matrix(T), ws: anytype, sp: CasSpace, gpa: Allocator) !void {
    const B_mo, const dB_mo = .{ b.B, b.dB };

    const nmo, const nact, const naux = .{ sp.nmo, sp.nact, B_mo.shape[2] };

    const v_act, const W = .{ ws.v_act, ws.W };

    var dv_act = try Vector(T).initZero(naux, gpa);
    defer dv_act.deinit(gpa);

    for (0..nact) |u| for (0..nact) |v| {
        const val_D = D_act.at(u, v);

        if (val_D == 0) continue;

        const U, const V = .{ sp.ncore + u, sp.ncore + v };

        for (0..naux) |P| {
            dv_act.ptr(P).* += val_D * dB_mo.at(.{ U, V, P });
        }
    };

    var dW = try Tensor(T, 3).initZero(.{ nact, sp.ncore, naux }, gpa);
    defer dW.deinit(gpa);

    for (0..nact) |u| for (0..nact) |v| {
        const val_D = D_act.at(u, v);

        if (val_D == 0) continue;

        const U = sp.ncore + u;

        for (0..sp.ncore) |j| for (0..naux) |P| {
            dW.ptr(.{ v, j, P }).* += val_D * dB_mo.at(.{ U, j, P });
        };
    };

    for (0..nmo) |p| for (0..sp.ncore) |j| {
        var sum = dF_core.at(p, j);

        for (0..naux) |P| {
            sum += dB_mo.at(.{ p, j, P }) * v_act.at(P) + B_mo.at(.{ p, j, P }) * dv_act.at(P);
        }

        for (0..nact) |v| for (0..naux) |P| {
            sum -= 0.5 * (dB_mo.at(.{ p, sp.ncore + v, P }) * W.at(.{ v, j, P }) + B_mo.at(.{ p, sp.ncore + v, P }) * dW.at(.{ v, j, P }));
        };

        dF.ptr(p, j).* = 2 * sum;
    };
}

/// Constructs derivative of active subspace columns of generalized Fock matrix from 3-center fitting derivatives.
fn buildActColumnsDerivRi(comptime T: type, dF: *Matrix(T), dF_core: Matrix(T), b: FittingDeriv(T), rdm: anytype, Q_aux: Tensor(T, 3), space: CasSpace, gpa: Allocator) !void {
    const B_mo, const dB_mo = .{ b.B, b.dB };

    const D_act, const d_act = rdm;

    const nmo, const nact, const naux = .{ space.nmo, space.nact, B_mo.shape[2] };

    var dQ_aux = try Tensor(T, 3).initZero(.{ nact, nact, naux }, gpa);
    defer dQ_aux.deinit(gpa);

    for (0..nact) |u| for (0..nact) |v| for (0..nact) |w| for (0..nact) |x| {
        const val_d = d_act.at(.{ u, v, w, x });

        if (val_d == 0) continue;

        const W_idx, const X_idx = .{ space.ncore + w, space.ncore + x };

        for (0..naux) |P| {
            dQ_aux.ptr(.{ u, v, P }).* += val_d * dB_mo.at(.{ W_idx, X_idx, P });
        }
    };

    for (0..nmo) |p| for (0..nact) |u| {
        const U, var sum: T = .{ space.ncore + u, 0 };

        for (0..nact) |v| {
            sum += dF_core.at(p, space.ncore + v) * D_act.at(v, u);
        }

        for (0..nact) |v| for (0..naux) |P| {
            const q = space.ncore + v;

            sum += dB_mo.at(.{ p, q, P }) * Q_aux.at(.{ u, v, P }) + B_mo.at(.{ p, q, P }) * dQ_aux.at(.{ u, v, P });
        };

        dF.ptr(p, U).* = sum;
    };
}
