//! Solves the Configuration Interaction equations by diagonalizing the Hamiltonian over a Slater determinant basis.

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
const Value = @import("value.zig").Value;
const Vector = @import("tensor.zig").Vector;

const ao2mo_pp = @import("integral_transform.zig").ao2mo_pp;
const ao2mo_pppp = @import("integral_transform.zig").ao2mo_pppp;
const ao2so_coef = @import("integral_transform.zig").ao2so_coef;
const ao2so_pp = @import("integral_transform.zig").ao2so_pp;
const ao2so_pppp = @import("integral_transform.zig").ao2so_pppp;
const nuclearRepulsionGradient = @import("hartree_fock.zig").nuclearRepulsionGradient;
const primType = @import("value.zig").primType;
const printf = @import("read_write.zig").printf;
const syevd = @import("linear_algebra.zig").syevd;

/// Configurations for CI determinant subspaces, active spaces, optimization, and derivative settings.
pub const Options = struct {
    frequency: ?FrequencyOptions = null,
    gradient: ?GradientOptions = null,
    hartree_fock: HartreeFockOptions,
    hessian: ?HessianOptions = null,
    nacv: ?NacvOptions = null,
    nstate: u32 = 1,
    nthreads: u32 = 1,
    optimize: ?OptimizeOptions = null,
    space: SpaceOptions = .{ .excitations = &.{ 1, 2 } },
    write: Write = .{},
};

/// Active space dimensions specifying active electrons and active orbitals for CASCI calculations.
pub const CasSpaceOptions = struct {
    active_electrons: u32,
    active_orbitals: u32,
};

/// Options for computing CI energy gradients analytically or numerically for a specific electronic state.
pub const GradientOptions = union(enum) {
    analytic: GradientAnalyticOptions,
    numeric: GradientNumericOptions,
};

/// Tagged union specifying CI nuclear Hessian calculation methods.
pub const HessianOptions = union(enum) {
    numeric: HessianNumericOptions,
};

/// Tagged union specifying CI nonadiabatic coupling vector calculation methods.
pub const NacvOptions = union(enum) {
    numeric: NacvNumericOptions,
};

/// Tagged union specifying geometry optimization algorithms.
pub const OptimizeOptions = union(enum) {
    bfgs: BfgsOptions,
    steepest_descent: SteepestDescentOptions,
};

/// Tagged union specifying determinant subspace selection via CAS or excitation levels.
pub const SpaceOptions = union(enum) {
    cas: CasSpaceOptions,
    excitations: []const u32,
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

/// Configuration options for steepest descent molecular geometry optimization.
const SteepestDescentOptions = struct {
    gradient: GradientOptions = .{ .analytic = .{} },
    iterations: u32 = 100,
    step: f64 = 1e-1,
    threshold: f64 = 1e-4,
};

/// Output file destinations for CI geometry, gradients, Hessians, and nonadiabatic coupling vectors.
const Write = struct {
    geometry: ?[]const u8 = null,
    gradient: ?[]const u8 = null,
    hessian: ?[]const u8 = null,
    nacv: ?[]const u8 = null,
};

/// Holds CI results: reference SCF results, state energies, wavefunctions, and optional nuclear derivatives.
pub fn Result(comptime T: type) type {
    return struct {
        hartree_fock: HartreeFockResult(T),

        C_ci: Matrix(T),

        energy: []T,

        grad: []Matrix(T) = &.{},
        hess: []Matrix(T) = &.{},
        nacv: []Matrix(T) = &.{},

        /// Frees all allocated memory associated with the CI results.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.hartree_fock.deinit(gpa);

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

/// Aligns molecular orbital and CI state phases with a reference step to preserve wavefunction continuity.
pub fn alignPhases(comptime T: type, C_mo: *Matrix(T), C_ci: *Matrix(T), mo_prev: ?Matrix(T), ci_prev: ?Matrix(T), dets: []const []const usize, gen: bool, flipped: ?[]bool) void {
    if (mo_prev) |prev| {
        for (0..C_mo.ncol()) |p| {
            var mo_overlap: T = 0;

            for (0..C_mo.nrow()) |mu| {
                mo_overlap += C_mo.at(mu, p) * prev.at(mu, p);
            }

            if (mo_overlap < 0) {
                for (0..C_mo.nrow()) |mu| {
                    C_mo.ptr(mu, p).* = -C_mo.at(mu, p);
                }

                for (dets, 0..) |det, d| {
                    const indexOf = std.mem.indexOfScalar;

                    var flip_det = (indexOf(usize, det, 2 * p) != null) != (indexOf(usize, det, 2 * p + 1) != null);

                    if (gen) {
                        flip_det = indexOf(usize, det, p) != null;
                    }

                    if (flip_det) for (0..C_ci.ncol()) |k| {
                        C_ci.ptr(d, k).* = -C_ci.at(d, k);
                    };
                }
            }
        }
    }

    if (ci_prev) |prev| {
        for (0..C_ci.ncol()) |k| {
            var overlap: T = 0;

            for (0..C_ci.nrow()) |d| {
                overlap += C_ci.at(d, k) * prev.at(d, k);
            }

            if (overlap < 0) {
                if (flipped) |f| if (k < f.len) {
                    f[k] = true;
                };

                for (0..C_ci.nrow()) |d| {
                    C_ci.ptr(d, k).* = -C_ci.at(d, k);
                }
            }
        }
    }
}

/// Builds active-space 1-particle and 2-particle reduced density matrices from CAS-CI coefficients.
pub fn buildActiveRdm(comptime T: type, dets: []const []const usize, C_ci: Matrix(T), state: usize, ncore: usize, nact: usize, gpa: Allocator) !struct { Matrix(T), Tensor(T, 4) } {
    var D_act = try Matrix(T).initZero(nact, nact, gpa);
    errdefer D_act.deinit(gpa);

    var d_act = try Tensor(T, 4).initZero(.{ nact, nact, nact, nact }, gpa);
    errdefer d_act.deinit(gpa);

    for (0..dets.len) |i| for (0..dets.len) |j| {
        const c_prod = C_ci.at(i, state) * C_ci.at(j, state);

        if (@abs(c_prod) < 1e-15) continue;

        const diff = findDeterminantDifferences(dets[i], dets[j]) orelse continue;

        const sA = signOfExcitations(dets[i], diff.a[0..diff.len]);
        const sB = signOfExcitations(dets[j], diff.b[0..diff.len]);

        const phase: T = if (diff.len > 0) @as(T, @floatFromInt(sA * sB)) else 1;

        switch (diff.len) {
            0 => accumulateRdmDiagonal(T, &D_act, &d_act, dets[i], c_prod, ncore, nact),
            1 => accumulateRdmSingle(T, &D_act, &d_act, dets[i], diff.a[0], diff.b[0], phase, c_prod, ncore, nact),
            2 => accumulateRdmDouble(T, &d_act, diff.a, diff.b, phase, c_prod, ncore, nact),
            else => {},
        }
    };

    return .{ D_act, d_act };
}

/// Evaluates electronic Hamiltonian matrix elements across Slater determinants using Slater-Condon rules.
pub fn buildHamiltonianMatrix(comptime T: type, dets: []const []const usize, h_so: Matrix(T), g_so: Tensor(T, 4), gpa: Allocator) !Matrix(T) {
    var H_CI = try Matrix(T).init(dets.len, dets.len, gpa);
    errdefer H_CI.deinit(gpa);

    for (0..dets.len) |m| for (m..dets.len) |n_det| {
        const val = slater(T, dets[m], dets[n_det], h_so, g_so);

        H_CI.ptr(m, n_det).* = val;
        H_CI.ptr(n_det, m).* = val;
    };

    return H_CI;
}

/// Evaluates Hamiltonian matrix elements in the determinant basis and diagonalizes to get CI state energies.
pub fn computeCi(comptime T: type, io: ?std.Io, C: Matrix(T), ints: anytype, dets: []const []const usize, VN: T, generalized: bool, gpa: Allocator) !struct { Vector(T), Matrix(T) } {
    var timer: std.Io.Timestamp = if (io) |out| std.Io.Timestamp.now(out, .real) else undefined;

    var h_so, var g_so = try transformInts(T, C, ints.H.?, ints.g.?, generalized, gpa);

    defer {
        h_so.deinit(gpa);
        g_so.deinit(gpa);
    }

    if (io) |out| {
        try printf(out, " {f}\nCI MATRIX ASSEMBLY TIME:", .{timer.untilNow(out, .real)});

        timer = std.Io.Timestamp.now(out, .real);
    }

    var H_CI = try buildHamiltonianMatrix(T, dets, h_so, g_so, gpa);
    defer H_CI.deinit(gpa);

    if (io) |out| {
        try printf(out, " {f}\nCI DIAGONALIZATION TIME:", .{timer.untilNow(out, .real)});
        timer = std.Io.Timestamp.now(out, .real);
    }

    var E, var C_ci = try solveEigenvalueProblem(T, H_CI, VN, gpa);

    errdefer {
        E.deinit(gpa);

        C_ci.deinit(gpa);
    }

    if (io) |out| {
        try printf(out, " {f}\n", .{timer.untilNow(out, .real)});
    }

    return .{ E, C_ci };
}

/// Compares two Slater determinants and extracts orbital differences when excitation degree is at most two.
pub fn findDeterminantDifferences(A: []const usize, B: []const usize) ?struct { a: [2]usize, b: [2]usize, len: usize } {
    var diff_A: [2]usize = undefined;
    var diff_B: [2]usize = undefined;

    var diff_a_len: usize = 0;
    var diff_b_len: usize = 0;

    var idx_a: usize = 0;
    var idx_b: usize = 0;

    while (idx_a < A.len and idx_b < B.len) {
        switch (std.math.order(A[idx_a], B[idx_b])) {
            .lt => {
                if (diff_a_len == 2) return null;

                diff_A[diff_a_len] = A[idx_a];

                diff_a_len, idx_a = .{ diff_a_len + 1, idx_a + 1 };
            },
            .gt => {
                if (diff_b_len == 2) return null;

                diff_B[diff_b_len] = B[idx_b];

                diff_b_len, idx_b = .{ diff_b_len + 1, idx_b + 1 };
            },
            .eq => {
                idx_a += 1;
                idx_b += 1;
            },
        }
    }

    while (idx_a < A.len) : (idx_a += 1) {
        if (diff_a_len == 2) return null;

        diff_A[diff_a_len] = A[idx_a];

        diff_a_len += 1;
    }

    while (idx_b < B.len) : (idx_b += 1) {
        if (diff_b_len == 2) return null;

        diff_B[diff_b_len] = B[idx_b];

        diff_b_len += 1;
    }

    if (diff_a_len != diff_b_len) return null;

    return .{ .a = diff_A, .b = diff_B, .len = diff_a_len };
}

/// Generates complete active space determinants preserving closed-shell core occupations.
pub fn generateCasDets(ncore: usize, nact: usize, nact_el: usize, mult: u32, gpa: Allocator) !std.ArrayList([]const usize) {
    var dets = std.ArrayList([]const usize).empty;

    errdefer {
        for (0..dets.items.len) |i| gpa.free(dets.items[i]);

        dets.deinit(gpa);
    }

    if ((nact_el + mult - 1) % 2 != 0) return error.InvalidInput;

    const n_alpha = (nact_el + mult - 1) / 2;

    const n_beta = nact_el - n_alpha;

    if (n_alpha > nact or n_beta > nact) return error.InvalidInput;

    var alpha_combs = try generateCombinations(nact, n_alpha, 0, gpa);

    defer {
        for (0..alpha_combs.items.len) |i| gpa.free(alpha_combs.items[i]);

        alpha_combs.deinit(gpa);
    }

    var beta_combs = try generateCombinations(nact, n_beta, 0, gpa);

    defer {
        for (0..beta_combs.items.len) |i| gpa.free(beta_combs.items[i]);

        beta_combs.deinit(gpa);
    }

    const nel = 2 * ncore + nact_el;

    for (alpha_combs.items) |ac| for (beta_combs.items) |bc| {
        const det = try gpa.alloc(usize, nel);
        errdefer gpa.free(det);

        for (0..2 * ncore) |k| {
            det[k] = k;
        }

        for (ac, 0..) |a, k| {
            det[2 * ncore + k] = 2 * (ncore + a);
        }

        for (bc, 0..) |b, k| {
            det[2 * ncore + n_alpha + k] = 2 * (ncore + b) + 1;
        }

        std.mem.sort(usize, det, {}, std.sort.asc(usize));

        try dets.append(gpa, det);
    };

    return dets;
}

/// Generates all k-combinations from a set of size n, representing orbital indices.
pub fn generateCombinations(n: usize, k: usize, offset: usize, gpa: Allocator) !std.ArrayList([]const usize) {
    var results: std.ArrayList([]const usize) = .empty;

    errdefer {
        for (0..results.items.len) |i| gpa.free(results.items[i]);

        results.deinit(gpa);
    }

    if (k > n) return results;

    if (k == 0) {
        try results.append(gpa, try gpa.alloc(usize, 0));

        return results;
    }

    var current = try gpa.alloc(usize, k);
    defer gpa.free(current);

    for (0..k) |i| {
        current[i] = i;
    }

    while (true) {
        const copy = try gpa.alloc(usize, k);
        errdefer gpa.free(copy);

        for (0..k) |i| {
            copy[i] = current[i] + offset;
        }

        try results.append(gpa, copy);

        var i: usize = k - 1;

        while (current[i] == n - k + i) : (i -= 1) {
            if (i == 0) return results;
        }

        current[i] += 1;

        for (i + 1..k) |j| {
            current[j] = current[j - 1] + 1;
        }
    }

    return results;
}

/// Generates excited Slater determinants up to specified excitation levels from a Hartree-Fock reference.
pub fn generateDets(nel: usize, nsp: usize, excitations: []const u32, gpa: Allocator) !std.ArrayList([]const usize) {
    var dets = std.ArrayList([]const usize).empty;

    errdefer {
        for (0..dets.items.len) |i| gpa.free(dets.items[i]);

        dets.deinit(gpa);
    }

    {
        const ref_det = try gpa.alloc(usize, nel);

        for (0..nel) |i| {
            ref_det[i] = i;
        }

        try dets.append(gpa, ref_det);
    }

    for (0..excitations.len) |i| {
        const k = excitations[i];

        if (k == 0 or k > nel or k > nsp - nel) continue;

        var occ_combs = try generateCombinations(nel, k, 0, gpa);

        defer {
            for (0..occ_combs.items.len) |j| gpa.free(occ_combs.items[j]);

            occ_combs.deinit(gpa);
        }

        var vir_combs = try generateCombinations(nsp - nel, k, nel, gpa);

        defer {
            for (0..vir_combs.items.len) |j| gpa.free(vir_combs.items[j]);

            vir_combs.deinit(gpa);
        }

        for (0..occ_combs.items.len) |j| for (0..vir_combs.items.len) |l| {
            var list = try std.ArrayList(usize).initCapacity(gpa, nel);
            errdefer list.deinit(gpa);

            for (0..nel) |p| if (std.mem.indexOfScalar(usize, occ_combs.items[j], p) == null) {
                list.appendAssumeCapacity(p);
            };

            list.appendSliceAssumeCapacity(vir_combs.items[l]);

            std.mem.sort(usize, list.items, {}, std.sort.asc(usize));

            try dets.append(gpa, try list.toOwnedSlice(gpa));
        };
    }

    return dets;
}

/// Generates all unique upper-triangular state index pairs for nonadiabatic coupling evaluation.
pub fn generatePairs(n: usize, gpa: Allocator) ![][2]u32 {
    const pairs = try gpa.alloc([2]u32, n * (n - 1) / 2);

    var idx: usize = 0;

    for (0..n) |i| for (i + 1..n) |j| {
        pairs[idx] = .{ @intCast(i), @intCast(j) };

        idx += 1;
    };

    return pairs;
}

/// Computes the permutation sign (+1 or -1) associated with mapping two determinants via orbital excitations.
pub fn signOfExcitations(S: []const usize, U: []const usize) i2 {
    var exp: usize = 0;

    for (0..U.len) |j| {
        exp += std.mem.indexOfScalar(usize, S, U[j]).? - j;
    }

    return if (exp % 2 == 0) 1 else -1;
}

/// Computes the Hamiltonian matrix element between two Slater determinants using the Slater-Condon rules.
pub fn slater(comptime T: type, A: []const usize, B: []const usize, h_so: Matrix(T), g_so: Tensor(T, 4)) T {
    const diff = findDeterminantDifferences(A, B) orelse return Value(T).fromFloat(0).val;

    const s_A = signOfExcitations(A, diff.a[0..diff.len]);
    const s_B = signOfExcitations(B, diff.b[0..diff.len]);

    const phase = @as(primType(T), @floatFromInt(s_A * s_B));

    return evaluateSlaterCondon(T, A, diff.a, diff.b, diff.len, phase, h_so, g_so);
}

/// Solves the Configuration Interaction Hamiltonian eigenvalue problem over the determinant space.
pub fn solve(comptime T: type, io: std.Io, opt: Options, sys: *MolecularSystem(T), ref: HartreeFockResult(T), log: bool, gpa: Allocator) !Result(T) {
    _ = sys;

    var hfres = ref;
    errdefer hfres.deinit(gpa);

    cblas.openblas_set_num_threads(@intCast(opt.nthreads));

    const generalized = opt.hartree_fock.generalized;

    const nel = hfres.ints.sys.nel;

    const nmo = if (generalized) hfres.C_mo.ncol() / 2 else hfres.C_mo.ncol();

    try validateSpaceOptions(opt, nel, nmo);

    var timer = std.Io.Timestamp.now(io, .real);

    if (log) {
        switch (opt.space) {
            .cas => |c| {
                const fmt = "\nCASCI SPACE: CAS({d}, {d})\n\nSLATER DETERMINANT TIME:";

                try printf(io, fmt, .{ c.active_electrons, c.active_orbitals });
            },
            .excitations => |ex| {
                try printf(io, "\nCI EXCITATIONS CONSIDERED: [", .{});

                for (0..ex.len) |i| {
                    try printf(io, "{d}{s}", .{ ex[i], if (i == ex.len - 1) "]" else ", " });
                }

                try printf(io, "\n\nSLATER DETERMINANT TIME:", .{});
            },
        }
    }

    timer = std.Io.Timestamp.now(io, .real);

    var dets = switch (opt.space) {
        .cas => |c| blk: {
            const ncore = (hfres.ints.sys.nel - c.active_electrons) / 2;

            const mult = opt.hartree_fock.multiplicity;

            break :blk try generateCasDets(ncore, c.active_orbitals, c.active_electrons, mult, gpa);
        },
        .excitations => |ex| try generateDets(hfres.ints.sys.nel, 2 * hfres.ints.sys.nbf, ex, gpa),
    };

    defer {
        for (0..dets.items.len) |i| gpa.free(dets.items[i]);

        dets.deinit(gpa);
    }

    if (opt.nstate > dets.items.len) {
        std.log.err("NUMBER OF REQUESTED CI STATES EXCEEDS DETERMINANT COUNT", .{});

        return error.InvalidInput;
    }

    if (log) {
        const params = .{ timer.untilNow(io, .real), dets.items.len };

        try printf(io, " {f}\n\nNUMBER OF DETERMINANTS: {d}\n\nORBITAL TRANSFORMS TIME:", params);
    }

    const nrep = try hfres.ints.sys.nrep();

    const io_ci = if (log) io else null;

    var E, var C_ci = try computeCi(T, io_ci, hfres.C_mo, hfres.ints, dets.items, nrep, generalized, gpa);

    errdefer {
        E.deinit(gpa);

        C_ci.deinit(gpa);
    }

    if (log) {
        const n_states = @min(opt.nstate, E.length());

        if (n_states == 1) {
            try printf(io, "\nFINAL CI ENERGY: {d:.14} Eh\n", .{E.at(0)});
        }

        if (n_states > 1) {
            try printf(io, "\n", .{});

            for (0..n_states) |s| {
                try printf(io, "STATE {d:2} FINAL CI ENERGY: {d:.14} Eh\n", .{ s, E.at(s) });
            }
        }
    }

    var grad: []Matrix(T) = &.{};

    errdefer {
        for (0..grad.len) |k| grad[k].deinit(gpa);

        if (grad.len > 0) gpa.free(grad);
    }

    if (opt.gradient) |g| switch (g) {
        .analytic => |a| grad = try calculateNuclearGradient(T, hfres, C_ci, dets, a.states orelse &[_]u32{0}, gpa),
        .numeric => {},
    };

    return Result(T){
        .hartree_fock = hfres,

        .C_ci = C_ci,
        .grad = grad,
        .hess = &.{},
        .nacv = &.{},

        .energy = E.data,
    };
}

/// Diagonalizes the CI Hamiltonian matrix to yield state energies and coefficients.
pub fn solveEigenvalueProblem(comptime T: type, H_CI: Matrix(T), VN: T, gpa: Allocator) !struct { Vector(T), Matrix(T) } {
    var E = try Vector(T).init(H_CI.shape[0], gpa);
    errdefer E.deinit(gpa);

    var C_ci = try Matrix(T).initZero(H_CI.shape[0], H_CI.shape[0], gpa);
    errdefer C_ci.deinit(gpa);

    try syevd(T, &E, &C_ci, H_CI);

    for (0..E.length()) |i| {
        E.ptr(i).* += VN;
    }

    return .{ E, C_ci };
}

/// Transforms one- and two-electron integrals from the atomic orbital basis to the molecular spin-orbital basis.
pub fn transformInts(comptime T: type, C: Matrix(T), H: Matrix(T), g: Tensor(T, 4), generalized: bool, gpa: Allocator) !struct { Matrix(T), Tensor(T, 4) } {
    const nsp = if (generalized) C.shape[0] else 2 * C.shape[0];

    var C_SO = if (generalized) C else try Matrix(T).initZero(nsp, nsp, gpa);
    defer if (!generalized) C_SO.deinit(gpa);

    if (!generalized) {
        ao2so_coef(T, &C_SO, C);
    }

    var H_SO = if (generalized) H else try Matrix(T).initZero(nsp, nsp, gpa);
    defer if (!generalized) H_SO.deinit(gpa);

    if (!generalized) {
        ao2so_pp(T, &H_SO, H);
    }

    var g_SO = if (generalized) g else try Tensor(T, 4).initZero(.{ nsp, nsp, nsp, nsp }, gpa);
    defer if (!generalized) g_SO.deinit(gpa);

    if (!generalized) {
        ao2so_pppp(T, &g_SO, g);
    }

    var h_so = try Matrix(T).init(nsp, nsp, gpa);
    errdefer h_so.deinit(gpa);

    ao2mo_pp(T, &h_so, H_SO, C_SO);

    var g_so = try Tensor(T, 4).init(.{ nsp, nsp, nsp, nsp }, gpa);
    errdefer g_so.deinit(gpa);

    try ao2mo_pppp(T, &g_so, g_SO, C_SO, gpa);

    return .{ h_so, g_so };
}

/// Accumulates diagonal 1-RDM and 2-RDM contributions between identical configurations.
fn accumulateRdmDiagonal(comptime T: type, D_act: *Matrix(T), d_act: *Tensor(T, 4), det: []const usize, c_prod: T, ncore: usize, nact: usize) void {
    for (0..nact) |u| {
        const pa = 2 * (ncore + u) + 0;
        const pb = 2 * (ncore + u) + 1;

        const occ_a = std.mem.indexOfScalar(usize, det, pa) != null;
        const occ_b = std.mem.indexOfScalar(usize, det, pb) != null;

        if (occ_a) D_act.ptr(u, u).* += c_prod;
        if (occ_b) D_act.ptr(u, u).* += c_prod;

        if (occ_a and occ_b) d_act.ptr(.{ u, u, u, u }).* += 2 * c_prod;

        for (0..nact) |w| if (u != w) {
            const qa = 2 * (ncore + w) + 0;
            const qb = 2 * (ncore + w) + 1;

            const w_a = std.mem.indexOfScalar(usize, det, qa) != null;
            const w_b = std.mem.indexOfScalar(usize, det, qb) != null;

            const nu: T = (if (occ_a) @as(T, 1) else @as(T, 0)) + (if (occ_b) @as(T, 1) else @as(T, 0));
            const nw: T = (if (w_a) @as(T, 1) else @as(T, 0)) + (if (w_b) @as(T, 1) else @as(T, 0));

            d_act.ptr(.{ u, u, w, w }).* += nu * nw * c_prod;

            var exch: T = 0;

            if (occ_a and w_a) exch += 1;
            if (occ_b and w_b) exch += 1;

            d_act.ptr(.{ u, w, w, u }).* -= exch * c_prod;
        };
    }
}

/// Accumulates 2-RDM contributions between configurations differing by two spin-orbitals.
fn accumulateRdmDouble(comptime T: type, d_act: *Tensor(T, 4), diff_a: [2]usize, diff_b: [2]usize, phase: T, c_prod: T, ncore: usize, nact: usize) void {
    const p, const q = .{ diff_a[0], diff_a[1] };
    const r, const s = .{ diff_b[0], diff_b[1] };

    if (p >= 2 * ncore and q >= 2 * ncore and r >= 2 * ncore and s >= 2 * ncore) {
        const u = (p / 2) - ncore;
        const v = (q / 2) - ncore;
        const w = (r / 2) - ncore;
        const x = (s / 2) - ncore;

        if (u < nact and v < nact and w < nact and x < nact) {
            if ((p % 2 == r % 2) and (q % 2 == s % 2)) {
                d_act.ptr(.{ u, w, v, x }).* += c_prod * phase;
                d_act.ptr(.{ v, x, u, w }).* += c_prod * phase;
            }

            if ((p % 2 == s % 2) and (q % 2 == r % 2)) {
                d_act.ptr(.{ u, x, v, w }).* -= c_prod * phase;
                d_act.ptr(.{ v, w, u, x }).* -= c_prod * phase;
            }
        }
    }
}

/// Accumulates 1-RDM and 2-RDM contributions between configurations differing by one spin-orbital.
fn accumulateRdmSingle(comptime T: type, D_act: *Matrix(T), d_act: *Tensor(T, 4), det: []const usize, p: usize, q: usize, phase: T, c_prod: T, ncore: usize, nact: usize) void {
    if (p >= 2 * ncore and q >= 2 * ncore) {
        const u = (p / 2) - ncore;
        const v = (q / 2) - ncore;

        if (u < nact and v < nact and (p % 2 == q % 2)) {
            D_act.ptr(u, v).* += c_prod * phase;

            for (0..nact) |w| {
                const ka = 2 * (ncore + w) + 0;
                const kb = 2 * (ncore + w) + 1;

                const w_a = std.mem.indexOfScalar(usize, det, ka) != null;
                const w_b = std.mem.indexOfScalar(usize, det, kb) != null;

                const occ_ka = w_a and ka != p and ka != q;
                const occ_kb = w_b and kb != p and kb != q;

                const nw: T = (if (occ_ka) @as(T, 1) else @as(T, 0)) + (if (occ_kb) @as(T, 1) else @as(T, 0));

                d_act.ptr(.{ u, v, w, w }).* += nw * c_prod * phase;
                d_act.ptr(.{ w, w, u, v }).* += nw * c_prod * phase;

                var exch: T = 0;

                if ((p % 2 == 0) and occ_ka) exch += 1;
                if ((p % 2 == 1) and occ_kb) exch += 1;

                d_act.ptr(.{ u, w, w, v }).* -= exch * c_prod * phase;
                d_act.ptr(.{ w, v, u, w }).* -= exch * c_prod * phase;
            }
        }
    }
}

/// Computes analytic nuclear gradients of CI state energies using dual number differentiation.
fn calculateNuclearGradient(comptime T: type, hfres: HartreeFockResult(T), C_ci: Matrix(T), dets: std.ArrayList([]const usize), states: []const u32, gpa: Allocator) ![]Matrix(T) {
    const generalized = hfres.ints.sys.nbf != hfres.C_mo.shape[0];

    var grads = try gpa.alloc(Matrix(T), states.len);
    errdefer gpa.free(grads);

    var grads_done: usize = 0;

    errdefer for (0..grads_done) |k| grads[k].deinit(gpa);

    for (0..states.len) |idx| {
        grads[idx] = try nuclearRepulsionGradient(T, hfres.ints.sys, gpa);

        grads_done += 1;
    }

    var H_d = try Matrix(ScalarDual(T)).init(hfres.ints.H.?.nrow(), hfres.ints.H.?.ncol(), gpa);
    defer H_d.deinit(gpa);

    var g_d = try Tensor(ScalarDual(T), 4).init(hfres.ints.g.?.shape, gpa);
    defer g_d.deinit(gpa);

    var C_d = try Matrix(ScalarDual(T)).init(hfres.C_mo.nrow(), hfres.C_mo.ncol(), gpa);
    defer C_d.deinit(gpa);

    for (0..3 * hfres.ints.sys.atoms.len) |i| {
        for (0..H_d.nrow()) |mu| for (0..H_d.ncol()) |nu| {
            H_d.ptr(mu, nu).* = ScalarDual(T).init(hfres.ints.H.?.at(mu, nu), hfres.ints.dH.?.at(.{ i, mu, nu }));
        };

        const n = g_d.shape[0];

        for (0..n) |mu| for (0..n) |lambda| for (0..n) |nu| for (0..n) |sigma| {
            const indices = .{ mu, lambda, nu, sigma };

            const val = hfres.ints.g.?.at(indices);

            g_d.ptr(indices).* = ScalarDual(T).init(val, hfres.ints.dg.?.at(.{ i, mu, lambda, nu, sigma }));
        };

        for (0..C_d.nrow()) |mu| for (0..C_d.ncol()) |p| {
            C_d.ptr(mu, p).* = ScalarDual(T).init(hfres.C_mo.at(mu, p), hfres.dC_mo.?.at(.{ i, mu, p }));
        };

        var h_so_d, var g_so_d = try transformInts(ScalarDual(T), C_d, H_d, g_d, generalized, gpa);

        defer {
            h_so_d.deinit(gpa);
            g_so_d.deinit(gpa);
        }

        for (0..dets.items.len) |m| for (m..dets.items.len) |n_det| {
            const factor: T = if (m == n_det) 1 else 2;

            const val = slater(ScalarDual(T), dets.items[m], dets.items[n_det], h_so_d, g_so_d);

            for (states, 0..) |s, s_idx| {
                grads[s_idx].ptr(i / 3, i % 3).* += factor * C_ci.at(m, s) * C_ci.at(n_det, s) * val.der;
            }
        };
    }

    return grads;
}

/// Evaluates Hamiltonian matrix elements between Slater determinants using the Slater-Condon excitation rules.
fn evaluateSlaterCondon(comptime T: type, A: []const usize, diff_a: [2]usize, diff_b: [2]usize, len: usize, phase: primType(T), h_so: Matrix(T), g_so: Tensor(T, 4)) T {
    return switch (len) {
        0 => slaterDiagonal(T, h_so, g_so, A),
        1 => slaterSingle(T, h_so, g_so, A, diff_a[0], diff_b[0], phase),
        2 => slaterDouble(T, g_so, diff_a[0], diff_a[1], diff_b[0], diff_b[1], phase),
        else => Value(T).fromFloat(0).val,
    };
}

/// Evaluates diagonal Hamiltonian matrix elements between identical Slater determinants via Slater-Condon rules.
fn slaterDiagonal(comptime T: type, h_so: Matrix(T), g_so: Tensor(T, 4), det: []const usize) T {
    var sum = Value(T).fromFloat(0);

    for (0..det.len) |i| {
        sum = sum.add(Value(T).init(h_so.at(det[i], det[i])));
    }

    for (0..det.len) |i| for (i + 1..det.len) |j| {
        const p, const q = .{ det[i], det[j] };

        const coulomb = Value(T).init(g_so.at(.{ p, q, p, q }));

        const exchange = Value(T).init(g_so.at(.{ p, q, q, p }));

        sum = sum.add(coulomb.sub(exchange));
    };

    return sum.val;
}

/// Evaluates Hamiltonian matrix elements between Slater determinants differing by two spin-orbitals.
fn slaterDouble(comptime T: type, g_so: Tensor(T, 4), p: usize, q: usize, r: usize, s: usize, phase: primType(T)) T {
    const direct = Value(T).init(g_so.at(.{ p, q, r, s }));
    const exchange = Value(T).init(g_so.at(.{ p, q, s, r }));

    return direct.sub(exchange).muls(phase).val;
}

/// Evaluates Hamiltonian matrix elements between Slater determinants differing by a single spin-orbital.
fn slaterSingle(comptime T: type, h_so: Matrix(T), g_so: Tensor(T, 4), det: []const usize, p: usize, q: usize, phase: primType(T)) T {
    var elem = Value(T).init(h_so.at(p, q));

    for (0..det.len) |idx| if (det[idx] != p) {
        const k = det[idx];

        const coulomb = Value(T).init(g_so.at(.{ p, k, q, k }));
        const exchange = Value(T).init(g_so.at(.{ p, k, k, q }));

        elem = elem.add(coulomb.sub(exchange));
    };

    return elem.muls(phase).val;
}

/// Validates orbital space partitioning consistency for Complete Active Space calculations.
fn validateSpaceOptions(opt: Options, nel: usize, nmo: usize) !void {
    if (opt.space == .cas) {
        const c = opt.space.cas;

        const nact_el, const nact = .{ c.active_electrons, c.active_orbitals };

        if (nact_el > nel or (nel - nact_el) % 2 != 0) {
            std.log.err("ACTIVE ELECTRONS MUST BE LESS THAN OR EQUAL TO TOTAL ELECTRONS WITH CLOSED SHELL CORE", .{});

            return error.InvalidInput;
        }

        const ncore = (nel - nact_el) / 2;

        if (ncore + nact > nmo) {
            std.log.err("ACTIVE PLUS CORE ORBITALS EXCEED MOLECULAR ORBITAL COUNT", .{});

            return error.InvalidInput;
        }
    }
}
