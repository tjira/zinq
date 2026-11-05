//! Solves the Complete Active Space Self-Consistent Field equations to optimize orbitals and CI states.

const std = @import("std");

const cblas = @import("cimport.zig").cblas;

const Allocator = std.mem.Allocator;

const FrequencyOptions = @import("frequency_analysis.zig").Options;
const HartreeFockOptions = @import("hartree_fock.zig").Options;
const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const ScalarDual = @import("dual.zig").ScalarDual;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

const ao2mo_1e = @import("integral_transform.zig").ao2mo_1e;
const ao2mo_2e = @import("integral_transform.zig").ao2mo_2e;
const buildActRdm = @import("rdm.zig").buildActRdm;
const buildStateAveragedActRdm = @import("rdm.zig").buildStateAveragedActRdm;
const computeCi = @import("configuration_interaction.zig").computeCi;
const gemm = @import("linear_algebra.zig").gemm;
const generateCasDets = @import("determinant.zig").generateCasDets;
const getrf = @import("linear_algebra.zig").getrf;
const getrs = @import("linear_algebra.zig").getrs;
const printf = @import("read_write.zig").printf;
const slater = @import("slater_condon.zig").slater;
const syevd = @import("linear_algebra.zig").syevd;
const transformInts = @import("integral_transform.zig").transformInts;

/// Configurations for CASSCF active space dimensions, iterations, and derivative calculations.
pub const Options = struct {
    active_electrons: u32,
    active_orbitals: u32,
    frequency: ?FrequencyOptions = null,
    gradient: ?GradientOptions = null,
    hartree_fock: HartreeFockOptions,
    hessian: ?HessianOptions = null,
    iterations: u32 = 100,
    nacv: ?NacvOptions = null,
    nstate: u32 = 1,
    nthreads: u32 = 1,
    optimize: ?OptimizeOptions = null,
    orbital_threshold: f64 = 1e-6,
    threshold: f64 = 1e-8,
    weights: ?[]const f64 = null,
    write: Write = .{},
};

/// Options for computing CASSCF energy gradients analytically or numerically for a target electronic state.
pub const GradientOptions = union(enum) {
    analytic: GradientAnalyticOptions,
    numeric: GradientNumericOptions,
};

/// Tagged union specifying CASSCF nuclear Hessian calculation methods.
pub const HessianOptions = union(enum) {
    numeric: HessianNumericOptions,
};

/// Tagged union specifying CASSCF nonadiabatic coupling vector calculation methods.
pub const NacvOptions = union(enum) {
    numeric: NacvNumericOptions,
};

/// Tagged union specifying molecular geometry optimization algorithms for CASSCF.
pub const OptimizeOptions = union(enum) {
    bfgs: BfgsOptions,
    steepest_descent: SteepestDescentOptions,
};

/// Configuration options for BFGS quasi-Newton molecular geometry optimization.
const BfgsOptions = struct {
    gradient: GradientOptions = .{ .analytic = .{} },
    iterations: u32 = 100,
    step: f64 = 1,
    threshold: f64 = 1e-4,
};

/// Analytical nuclear gradient evaluation settings for specific electronic states.
const GradientAnalyticOptions = struct {
    states: ?[]const u32 = null,
};

/// Finite difference displacement settings for numerical gradient evaluation.
const GradientNumericOptions = struct {
    states: ?[]const u32 = null,
    step: f64 = 1e-5,
};

/// Finite difference displacement step for numerical Hessian evaluation.
const HessianNumericOptions = struct {
    states: ?[]const u32 = null,
    step: f64 = 1e-5,
};

/// Finite difference displacement step for numerical nonadiabatic coupling vector evaluation.
const NacvNumericOptions = struct {
    pairs: ?[]const [2]u32 = null,
    step: f64 = 1e-5,
};

/// Orbital index pair representing a non-redundant rotation parameter between subspaces.
const OrbitalPair = struct {
    p: usize,
    q: usize,
};

/// Configuration options for steepest descent molecular geometry optimization.
const SteepestDescentOptions = struct {
    gradient: GradientOptions = .{ .analytic = .{} },
    iterations: u32 = 100,
    step: f64 = 1e-1,
    threshold: f64 = 1e-4,
};

/// Output file destinations for CASSCF geometry, gradients, Hessians, NACV, and orbitals.
const Write = struct {
    geometry: ?[]const u8 = null,
    gradient: ?[]const u8 = null,
    hessian: ?[]const u8 = null,
    nacv: ?[]const u8 = null,
    orbitals: ?[]const u8 = null,
};

/// Holds CASSCF converged results: reference SCF data, state energies, molecular orbitals, and derivatives.
pub fn Result(comptime T: type) type {
    return struct {
        hartree_fock: HartreeFockResult(T),

        C_mo: Matrix(T),
        C_ci: Matrix(T),

        energy: []T,

        grad: []Matrix(T) = &.{},
        hess: []Matrix(T) = &.{},
        nacv: []Matrix(T) = &.{},

        /// Frees all allocated memory associated with the CASSCF results.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.hartree_fock.deinit(gpa);

            self.C_mo.deinit(gpa);
            self.C_ci.deinit(gpa);
            gpa.free(self.energy);

            for (0..self.grad.len) |i| {
                self.grad[i].deinit(gpa);
            }

            if (self.grad.len > 0) gpa.free(self.grad);

            for (0..self.hess.len) |i| {
                self.hess[i].deinit(gpa);
            }

            if (self.hess.len > 0) gpa.free(self.hess);

            for (0..self.nacv.len) |i| {
                self.nacv[i].deinit(gpa);
            }

            if (self.nacv.len > 0) gpa.free(self.nacv);
        }
    };
}

/// Solves the Complete Active Space Self-Consistent Field equations to optimize orbitals and CI states.
pub fn solve(comptime T: type, io: std.Io, opt: Options, sys: *MolecularSystem(T), ref: HartreeFockResult(T), log: bool, gpa: Allocator) !Result(T) {
    _ = sys;

    var hfres = ref;
    errdefer hfres.deinit(gpa);

    cblas.openblas_set_num_threads(@intCast(opt.nthreads));

    var timer = std.Io.Timestamp.now(io, .real);

    const nel, const nact_el, const nact = .{ hfres.ints.sys.nel, opt.active_electrons, opt.active_orbitals };

    if (nact_el > nel or (nel - nact_el) % 2 != 0) {
        std.log.err("ACTIVE ELECTRONS MUST BE LESS THAN OR EQUAL TO TOTAL ELECTRONS WITH CLOSED SHELL CORE", .{});

        return error.InvalidInput;
    }

    const ncore, const nmo = .{ (nel - nact_el) / 2, hfres.C_mo.ncol() };

    if (ncore + nact > nmo) {
        std.log.err("ACTIVE PLUS CORE ORBITALS EXCEED MOLECULAR ORBITAL COUNT", .{});

        return error.InvalidInput;
    }

    if (log) {
        const fmt = "\nHARTREE-FOCK CYCLE TIME: {f}\n\nCASSCF SPACE: CAS({d}, {d})\n\n";

        try printf(io, fmt, .{ timer.untilNow(io, .real), nact_el, nact });
    }

    var dets = try generateCasDets(ncore, nact, nact_el, opt.hartree_fock.multiplicity, gpa);

    defer {
        for (0..dets.items.len) |i| gpa.free(dets.items[i]);

        dets.deinit(gpa);
    }

    if (opt.nstate > dets.items.len) {
        std.log.err("NUMBER OF REQUESTED CASSCF STATES EXCEEDS DETERMINANT COUNT", .{});

        return error.InvalidInput;
    }

    if (log) {
        try printf(io, "NUMBER OF CAS DETERMINANTS: {d}\n", .{dets.items.len});

        const fmt = "\nCASSCF MACRO-ITERATION CYCLE\n{s:4} {s:20} {s:9} {s:9} {s:9}\n";

        try printf(io, fmt, .{ "ITER", "TOTAL ENERGY (Eh)", "|DE| (Eh)", "MAX(G)", "TIME" });
    }

    var C = try Matrix(T).init(hfres.C_mo.nrow(), hfres.C_mo.ncol(), gpa);
    errdefer C.deinit(gpa);

    @memcpy(C.data, hfres.C_mo.data);

    const weights = try gpa.alloc(T, opt.nstate);
    defer gpa.free(weights);

    if (opt.weights) |user_w| {
        var w_sum: f64 = 0;

        for (user_w) |w| w_sum += w;

        for (0..opt.nstate) |s| {
            weights[s] = @floatCast(user_w[s] / w_sum);
        }
    }

    if (opt.weights == null) {
        const uniform_w = 1.0 / @as(T, @floatFromInt(opt.nstate));

        @memset(weights, uniform_w);
    }

    try casscf(T, io, &C, hfres.ints, dets.items, weights, opt, log, gpa);

    try canonicalizeOrbitals(T, &C, hfres.ints, ncore, nact, gpa);

    const nrep = try hfres.ints.sys.nrep();

    var final_energies, var final_C_ci = try computeCi(T, null, C, hfres.ints, dets.items, nrep, false, gpa);

    errdefer {
        final_energies.deinit(gpa);

        final_C_ci.deinit(gpa);
    }

    if (log) {
        const n_states = @min(opt.nstate, final_energies.length());

        if (n_states == 1) {
            try printf(io, "\nFINAL CASSCF ENERGY: {d:.14} Eh\n", .{final_energies.at(0)});
        }

        if (n_states > 1) {
            try printf(io, "\n", .{});

            var final_sa_energy: T = 0;

            for (0..n_states) |s| {
                final_sa_energy += weights[s] * final_energies.at(s);

                try printf(io, "STATE {d:2} FINAL CASSCF ENERGY: {d:.14} Eh\n", .{ s, final_energies.at(s) });
            }

            try printf(io, "\nSTATE-AVERAGED FINAL CASSCF ENERGY: {d:.14} Eh\n", .{final_sa_energy});
        }
    }

    var grad: []Matrix(T) = &.{};

    errdefer {
        for (0..grad.len) |k| grad[k].deinit(gpa);

        if (grad.len > 0) gpa.free(grad);
    }

    if (opt.gradient) |g| switch (g) {
        .analytic => |a| {
            const states = a.states orelse &[_]u32{0};

            grad = try calculateNuclearGradient(T, hfres, C, final_C_ci, dets.items, weights, ncore, nact, states, gpa);
        },
        .numeric => {},
    };

    return Result(T){
        .hartree_fock = hfres,
        .energy = final_energies.data,
        .C_mo = C,
        .grad = grad,
        .hess = &.{},
        .nacv = &.{},
        .C_ci = final_C_ci,
    };
}

/// Applies a unitary Cayley transform to update molecular orbital coefficients from rotation generator kappa.
fn applyCayleyTransform(comptime T: type, C: *Matrix(T), kappa: Matrix(T), gpa: Allocator) !void {
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

/// Constructs the core Fock matrix by contracting one-electron integrals with closed-shell core orbitals.
fn buildCoreFock(comptime T: type, F_core: *Matrix(T), h_mo: Matrix(T), g_mo: Tensor(T, 4), ncore: usize) void {
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
fn buildDiagonalFock(comptime T: type, F_diag: *Vector(T), h_mo: Matrix(T), g_mo: Tensor(T, 4), D_act: Matrix(T), ncore: usize, gpa: Allocator) !void {
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
fn buildGeneralizedFock(comptime T: type, F_gen: *Matrix(T), h_mo: Matrix(T), g_mo: Tensor(T, 4), D_act: Matrix(T), d_act: Tensor(T, 4), ncore: usize, gpa: Allocator) !void {
    const nmo = h_mo.nrow();

    var F_core = try Matrix(T).init(nmo, nmo, gpa);
    defer F_core.deinit(gpa);

    buildCoreFock(T, &F_core, h_mo, g_mo, ncore);

    @memset(F_gen.data, 0);

    buildCoreColumns(T, F_gen, F_core, g_mo, D_act, ncore);

    buildActColumns(T, F_gen, F_core, g_mo, D_act, d_act, ncore);
}

/// Enumerates non-redundant orbital pairs across core-virtual, active-virtual, and core-active subspaces.
fn buildNonRedundantPairs(ncore: usize, nact: usize, nmo: usize, gpa: Allocator) ![]OrbitalPair {
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
fn buildOrbitalHessian(comptime T: type, A: *Matrix(T), h_mo: Matrix(T), g_mo: Tensor(T, 4), D_act: Matrix(T), d_act: Tensor(T, 4), pairs: []const OrbitalPair, ncore: usize, gpa: Allocator) !void {
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
fn buildOrbitalRotationGenerator(comptime T: type, kappa: *Matrix(T), g_orb: Matrix(T), F_diag: Vector(T), D_act: Matrix(T), ncore: usize) void {
    const nmo, const nact = .{ kappa.ncol(), D_act.nrow() };

    for (0..ncore) |j| for (ncore + nact..nmo) |b| {
        var denom = 4 * (F_diag.at(b) - F_diag.at(j));

        if (@abs(denom) < 0.1) denom = if (denom < 0) -0.1 else 0.1;

        var step = -g_orb.at(b, j) / denom;

        step = std.math.clamp(step, -0.2, 0.2);

        kappa.ptr(b, j).* = 0 + step;
        kappa.ptr(j, b).* = 0 - step;
    };

    for (0..nact) |u_idx| {
        const u = ncore + u_idx;

        const Duu = D_act.at(u_idx, u_idx);

        for (ncore + nact..nmo) |b| {
            var denom = 2 * Duu * (F_diag.at(b) - F_diag.at(u));

            if (@abs(denom) < 0.1) denom = if (denom < 0) -0.1 else 0.1;

            var step = -g_orb.at(b, u) / denom;

            step = std.math.clamp(step, -0.2, 0.2);

            kappa.ptr(b, u).* = 0 + step;
            kappa.ptr(u, b).* = 0 - step;
        }

        for (0..ncore) |j| {
            var denom = 2 * (2 - Duu) * (F_diag.at(u) - F_diag.at(j));

            if (@abs(denom) < 0.1) denom = if (denom < 0) -0.1 else 0.1;

            var step = -g_orb.at(u, j) / denom;

            step = std.math.clamp(step, -0.2, 0.2);

            kappa.ptr(u, j).* = 0 + step;
            kappa.ptr(j, u).* = 0 - step;
        }
    }
}

/// Computes analytic nuclear gradients of CASSCF state energies using dual number differentiation.
fn calculateNuclearGradient(comptime T: type, hfres: HartreeFockResult(T), C: Matrix(T), C_ci: Matrix(T), dets: []const []const usize, weights: []const T, ncore: usize, nact: usize, states: []const u32, gpa: Allocator) ![]Matrix(T) {
    var grads = try gpa.alloc(Matrix(T), states.len);
    errdefer gpa.free(grads);

    var grads_done: usize = 0;
    errdefer for (0..grads_done) |k| grads[k].deinit(gpa);

    for (0..states.len) |idx| {
        grads[idx] = try hfres.ints.sys.nuclearRepulsionGradient(gpa);

        grads_done += 1;
    }

    const nbf = C.nrow();
    const nmo = C.ncol();

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

/// Canonicalizes core and virtual orbital blocks by diagonalizing their respective Fock matrices.
fn canonicalizeOrbitals(comptime T: type, C: *Matrix(T), ints: anytype, ncore: usize, nact: usize, gpa: Allocator) !void {
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

/// Canonicalizes the unoccupied virtual orbital subspace by diagonalizing its core Fock sub-block.
fn canonicalizeVirtualSubspace(comptime T: type, C: *Matrix(T), F_core: Matrix(T), ncore: usize, nact: usize, gpa: Allocator) !void {
    const nmo = C.ncol();

    const nvirt = nmo - ncore - nact;

    if (nvirt <= 1) return;

    const nbf = C.nrow();
    const v_start = ncore + nact;

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

    gemm(T, &C_new, C_virt, V_vv, 1.0, 0.0, false, false);

    for (0..nbf) |mu| for (0..nvirt) |a| {
        C.ptr(mu, v_start + a).* = C_new.at(mu, a);
    };
}

/// Executes the two-step CASSCF macro-iteration cycle alternating configuration and orbital optimizations.
fn casscf(comptime T: type, io: std.Io, C: *Matrix(T), ints: anytype, dets: []const []const usize, weights: []const T, opt: Options, log: bool, gpa: Allocator) !void {
    const ncore, const nact = .{ (ints.sys.nel - opt.active_electrons) / 2, opt.active_orbitals };

    const nmo, var prev_energy: T = .{ C.ncol(), 0 };

    for (0..opt.iterations) |iter| {
        var iter_timer = std.Io.Timestamp.now(io, .real);

        var all_E, var C_ci = try computeCi(T, null, C.*, ints, dets, try ints.sys.nrep(), false, gpa);

        defer {
            all_E.deinit(gpa);

            C_ci.deinit(gpa);
        }

        var sa_energy: T = 0;

        for (0..opt.nstate) |I| {
            sa_energy += weights[I] * all_E.at(I);
        }

        const e_diff = if (iter == 0) sa_energy else @abs(sa_energy - prev_energy);

        prev_energy = sa_energy;

        var D_act, var d_act = try buildStateAveragedActRdm(T, dets, C_ci, weights, ncore, nact, gpa);

        defer {
            D_act.deinit(gpa);
            d_act.deinit(gpa);
        }

        var h_mo = try Matrix(T).init(nmo, nmo, gpa);
        defer h_mo.deinit(gpa);

        ao2mo_1e(T, &h_mo, ints.H.?, C.*);

        var g_mo = try Tensor(T, 4).init(.{ nmo, nmo, nmo, nmo }, gpa);
        defer g_mo.deinit(gpa);

        try ao2mo_2e(T, &g_mo, ints.g.?, C.*, gpa);

        var F_gen = try Matrix(T).initZero(nmo, nmo, gpa);
        defer F_gen.deinit(gpa);

        try buildGeneralizedFock(T, &F_gen, h_mo, g_mo, D_act, d_act, ncore, gpa);

        var g_orb = try Matrix(T).initZero(nmo, nmo, gpa);
        defer g_orb.deinit(gpa);

        const grad_norm = computeOrbitalGradient(T, &g_orb, F_gen, ncore, nact);

        if (log) {
            const fmt = "{d:4} {d:20.14} {e:9.3} {e:9.3} {s:9}\n";

            const time_str = try std.fmt.allocPrint(gpa, "{f}", .{iter_timer.untilNow(io, .real)});
            defer gpa.free(time_str);

            const de = if (iter == 0) 0 else e_diff;

            try printf(io, fmt, .{ iter + 1, sa_energy, de, grad_norm, time_str });
        }

        if (iter > 0 and e_diff < opt.threshold and grad_norm < opt.orbital_threshold) {
            break;
        }

        if (iter < opt.iterations - 1) {
            var F_diag = try Vector(T).init(nmo, gpa);
            defer F_diag.deinit(gpa);

            try buildDiagonalFock(T, &F_diag, h_mo, g_mo, D_act, ncore, gpa);

            var kappa = try Matrix(T).initZero(nmo, nmo, gpa);
            defer kappa.deinit(gpa);

            buildOrbitalRotationGenerator(T, &kappa, g_orb, F_diag, D_act, ncore);

            try applyCayleyTransform(T, C, kappa, gpa);
        }
    } else return error.CasscfDidNotConverge;
}

/// Evaluates the orbital gradient matrix and returns the maximum gradient norm among non-redundant pairs.
fn computeOrbitalGradient(comptime T: type, g_orb: *Matrix(T), F_gen: Matrix(T), ncore: usize, nact: usize) T {
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