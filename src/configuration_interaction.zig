//! Solves the Configuration Interaction equations by diagonalizing the Hamiltonian over a Slater determinant basis.

const std = @import("std");

const cblas = @import("cimport.zig").cblas;

const Allocator = std.mem.Allocator;

const FrequencyOptions = @import("frequency_analysis.zig").Options;
const HartreeFockOptions = @import("hartree_fock.zig").Options;
const HartreeFockResult = @import("hartree_fock.zig").Result;
const Integrals = @import("molecular_integrals.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const OptimizeOptions = @import("molecular_optimization.zig").OptimizeOptions(GradientOptions);
const ScalarDual = @import("dual.zig").ScalarDual;
const Tensor = @import("tensor.zig").Tensor;
const Value = @import("value.zig").Value;
const Vector = @import("tensor.zig").Vector;

const ao2mo_1e = @import("integral_transform.zig").ao2mo_1e;
const applySpinPenalty = @import("slater_condon.zig").applySpinPenalty;
const buildFittingIntegrals = @import("density_fitting.zig").buildFittingIntegrals;
const buildHamiltonianMatrix = @import("slater_condon.zig").buildHamiltonianMatrix;
const calculateCanonicalCiGradient = @import("ci_gradient.zig").calculateCanonicalCiGradient;
const calculateRiCiGradient = @import("ci_gradient.zig").calculateRiCiGradient;
const generateCasDets = @import("determinant.zig").generateCasDets;
const generateDets = @import("determinant.zig").generateDets;
const printDeterminants = @import("determinant.zig").printDeterminants;
const printf = @import("read_write.zig").printf;
const slater = @import("slater_condon.zig").slater;
const syevd = @import("linear_algebra.zig").syevd;
const transformInts = @import("integral_transform.zig").transformInts;
const transformIntsRi = @import("integral_transform.zig").transformIntsRi;

/// Configurations for CI determinant subspaces, active spaces, optimization, and derivative settings.
pub const Options = struct {
    aux_basis: ?[]const u8 = null,
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
    multiplicity: u32 = 1,
    spin_penalty: ?f64 = null,
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

/// Tagged union specifying determinant subspace selection via CAS or excitation levels.
pub const SpaceOptions = union(enum) {
    cas: CasSpaceOptions,
    excitations: []const u32,
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

/// Evaluates Hamiltonian matrix elements in the determinant basis and diagonalizes to get CI state energies.
pub fn computeCi(comptime T: type, io: ?std.Io, C: Matrix(T), ints: anytype, dets: []const []const usize, VN: T, generalized: bool, penalty: ?f64, s2_target: f64, gpa: Allocator) !struct { Vector(T), Matrix(T) } {
    var timer: std.Io.Timestamp = if (io) |out| std.Io.Timestamp.now(out, .real) else undefined;

    var h_ms, var g_ms = try transformInts(T, C, ints.H.?, ints.g.?, generalized, gpa);

    defer {
        h_ms.deinit(gpa);
        g_ms.deinit(gpa);
    }

    if (io) |out| {
        try printf(out, " {f}\nCI MATRIX ASSEMBLY TIME:", .{timer.untilNow(out, .real)});

        timer = std.Io.Timestamp.now(out, .real);
    }

    var H_ci = try buildHamiltonianMatrix(T, dets, h_ms, g_ms, gpa);
    defer H_ci.deinit(gpa);

    var H_orig = if (penalty != null) try Matrix(T).init(dets.len, dets.len, gpa) else undefined;
    defer if (penalty != null) H_orig.deinit(gpa);

    if (penalty) |pen| {
        @memcpy(H_orig.data, H_ci.data);

        applySpinPenalty(T, &H_ci, dets, pen, s2_target);
    }

    if (io) |out| {
        try printf(out, " {f}\nCI DIAGONALIZATION TIME:", .{timer.untilNow(out, .real)});
        timer = std.Io.Timestamp.now(out, .real);
    }

    var E, var C_ci = try solveEigenvalueProblem(T, H_ci, VN, gpa);

    errdefer {
        E.deinit(gpa);

        C_ci.deinit(gpa);
    }

    if (penalty != null) for (0..E.length()) |k| {
        var phys_e: T = 0;

        for (0..dets.len) |i| for (0..dets.len) |j| {
            const c_prod = Value(T).init(C_ci.at(i, k)).mul(Value(T).init(C_ci.at(j, k)));

            phys_e = Value(T).init(phys_e).add(c_prod.mul(Value(T).init(H_orig.at(i, j)))).val;
        };

        E.ptr(k).* = Value(T).init(phys_e).add(Value(T).init(VN)).val;
    };

    if (io) |out| {
        try printf(out, " {f}\n", .{timer.untilNow(out, .real)});
    }

    return .{ E, C_ci };
}

/// Computes CI state energies and coefficients using density fitting for active-space integrals.
pub fn computeCiRi(comptime T: type, dets: []const []const usize, C: Matrix(T), H_ao: Matrix(T), B_mo: Tensor(T, 3), VN: T, penalty: ?f64, s2_target: f64, ncore: usize, nact: usize, gpa: Allocator) !struct { Vector(T), Matrix(T) } {
    const nocc = ncore + nact;
    const nsp = 2 * nocc;
    const naux = B_mo.shape[2];

    var h_mo = try Matrix(T).init(C.ncol(), C.ncol(), gpa);
    defer h_mo.deinit(gpa);

    ao2mo_1e(T, &h_mo, H_ao, C);

    var h_ms = try Matrix(T).initZero(nsp, nsp, gpa);
    defer h_ms.deinit(gpa);

    for (0..nocc) |p| for (0..nocc) |q| {
        h_ms.ptr(2 * p, 2 * q).* = h_mo.at(p, q);
        h_ms.ptr(2 * p + 1, 2 * q + 1).* = h_mo.at(p, q);
    };

    var g_ms = try Tensor(T, 4).initZero(.{ nsp, nsp, nsp, nsp }, gpa);
    defer g_ms.deinit(gpa);

    for (0..nocc) |p| for (0..nocc) |r| {
        for (0..nocc) |q| for (0..nocc) |s| {
            var val: T = 0;

            for (0..naux) |P| {
                val += B_mo.at(.{ p, r, P }) * B_mo.at(.{ q, s, P });
            }

            g_ms.ptr(.{ 2 * p + 0, 2 * q + 0, 2 * r + 0, 2 * s + 0 }).* = val;
            g_ms.ptr(.{ 2 * p + 0, 2 * q + 1, 2 * r + 0, 2 * s + 1 }).* = val;
            g_ms.ptr(.{ 2 * p + 1, 2 * q + 0, 2 * r + 1, 2 * s + 0 }).* = val;
            g_ms.ptr(.{ 2 * p + 1, 2 * q + 1, 2 * r + 1, 2 * s + 1 }).* = val;
        };
    };

    var H_ci = try buildHamiltonianMatrix(T, dets, h_ms, g_ms, gpa);
    defer H_ci.deinit(gpa);

    var H_orig = if (penalty != null) try Matrix(T).init(dets.len, dets.len, gpa) else undefined;
    defer if (penalty != null) H_orig.deinit(gpa);

    if (penalty) |pen| {
        @memcpy(H_orig.data, H_ci.data);

        applySpinPenalty(T, &H_ci, dets, pen, s2_target);
    }

    var E, var C_ci = try solveEigenvalueProblem(T, H_ci, VN, gpa);

    errdefer {
        E.deinit(gpa);

        C_ci.deinit(gpa);
    }

    if (penalty != null) for (0..E.length()) |k| {
        var phys_e: T = 0;

        for (0..dets.len) |i| for (0..dets.len) |j| {
            const c_prod = Value(T).init(C_ci.at(i, k)).mul(Value(T).init(C_ci.at(j, k)));

            phys_e = Value(T).init(phys_e).add(c_prod.mul(Value(T).init(H_orig.at(i, j)))).val;
        };

        E.ptr(k).* = Value(T).init(phys_e).add(Value(T).init(VN)).val;
    };

    return .{ E, C_ci };
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

            break :blk try generateCasDets(ncore, c.active_orbitals, c.active_electrons, c.multiplicity, gpa);
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

    const penalty, const s2_target = switch (opt.space) {
        .cas => |c| blk: {
            const s_target = 0.5 * @as(f64, @floatFromInt(c.multiplicity - 1));

            break :blk .{ c.spin_penalty, s_target * (s_target + 1.0) };
        },
        .excitations => .{ null, 0.0 },
    };

    const nrep = try hfres.ints.sys.nrep();

    const io_ci = if (log) io else null;

    const analytic_grad = if (opt.gradient) |g| g == .analytic else false;

    if (opt.aux_basis) |aux| {
        try buildFittingIntegrals(T, io, &hfres.ints, aux, analytic_grad, opt.nthreads, gpa);
    }

    var E: Vector(T), var C_ci: Matrix(T) = .{ undefined, undefined };

    if (opt.aux_basis == null) {
        E, C_ci = try computeCi(T, io_ci, hfres.C_mo, hfres.ints, dets.items, nrep, generalized, penalty, s2_target, gpa);
    }

    if (opt.aux_basis != null) {
        E, C_ci = try computeCiWithRi(T, io_ci, hfres.C_mo, hfres.ints, dets.items, nrep, generalized, penalty, s2_target, gpa);
    }

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

            for (0..n_states) |k| {
                try printf(io, "STATE {d:2} FINAL CI ENERGY: {d:.14} Eh\n", .{ k, E.at(k) });
            }
        }

        const ncore = switch (opt.space) {
            .cas => |c| (hfres.ints.sys.nel - c.active_electrons) / 2,
            .excitations => 0,
        };

        const norbs = switch (opt.space) {
            .cas => |c| c.active_orbitals,
            .excitations => hfres.C_mo.ncol(),
        };

        try printDeterminants(T, io, C_ci, hfres.C_mo, hfres.ints.S.?, dets.items, n_states, ncore, norbs, generalized, gpa);
    }

    var grad: []Matrix(T) = &.{};

    errdefer {
        for (0..grad.len) |k| grad[k].deinit(gpa);

        if (grad.len > 0) gpa.free(grad);
    }

    if (opt.gradient) |g| switch (g) {
        .analytic => |a| {
            const states = a.states orelse &[_]u32{0};

            if (opt.aux_basis == null) {
                grad = try calculateCanonicalCiGradient(T, hfres, C_ci, dets.items, states, gpa);
            }

            if (opt.aux_basis != null) {
                grad = try calculateRiCiGradient(T, hfres, C_ci, dets.items, states, gpa);
            }
        },
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
pub fn solveEigenvalueProblem(comptime T: type, H_ci: Matrix(T), VN: T, gpa: Allocator) !struct { Vector(T), Matrix(T) } {
    var E = try Vector(T).init(H_ci.shape[0], gpa);
    errdefer E.deinit(gpa);

    var C_ci = try Matrix(T).initZero(H_ci.shape[0], H_ci.shape[0], gpa);
    errdefer C_ci.deinit(gpa);

    try syevd(T, &E, &C_ci, H_ci);

    for (0..E.length()) |i| {
        E.ptr(i).* += VN;
    }

    return .{ E, C_ci };
}

/// Computes CI state energies and coefficients using resolution of the identity three-center tensors.
fn computeCiWithRi(comptime T: type, io: ?std.Io, C: Matrix(T), ints: Integrals(T), dets: []const []const usize, VN: T, generalized: bool, penalty: ?f64, s2_target: f64, gpa: Allocator) !struct { Vector(T), Matrix(T) } {
    var timer = if (io) |out| std.Io.Timestamp.now(out, .real) else undefined;

    var h_ms, var g_ms = try transformIntsRi(T, C, ints.H.?, ints.B.?, generalized, gpa);

    defer {
        h_ms.deinit(gpa);
        g_ms.deinit(gpa);
    }

    if (io) |out| {
        try printf(out, " {f}\nCI MATRIX ASSEMBLY TIME:", .{timer.untilNow(out, .real)});

        timer = std.Io.Timestamp.now(out, .real);
    }

    var H_ci = try buildHamiltonianMatrix(T, dets, h_ms, g_ms, gpa);
    defer H_ci.deinit(gpa);

    var H_orig = if (penalty != null) try Matrix(T).init(dets.len, dets.len, gpa) else undefined;
    defer if (penalty != null) H_orig.deinit(gpa);

    if (penalty) |pen| {
        @memcpy(H_orig.data, H_ci.data);

        applySpinPenalty(T, &H_ci, dets, pen, s2_target);
    }

    if (io) |out| {
        try printf(out, " {f}\nCI DIAGONALIZATION TIME:", .{timer.untilNow(out, .real)});

        timer = std.Io.Timestamp.now(out, .real);
    }

    var E, var C_ci = try solveEigenvalueProblem(T, H_ci, VN, gpa);

    errdefer {
        E.deinit(gpa);

        C_ci.deinit(gpa);
    }

    if (penalty != null) for (0..E.length()) |k| {
        var phys_e: T = 0;

        for (0..dets.len) |i| for (0..dets.len) |j| {
            const c_prod = Value(T).init(C_ci.at(i, k)).mul(Value(T).init(C_ci.at(j, k)));

            phys_e = Value(T).init(phys_e).add(c_prod.mul(Value(T).init(H_orig.at(i, j)))).val;
        };

        E.ptr(k).* = Value(T).init(phys_e).add(Value(T).init(VN)).val;
    };

    if (io) |out| {
        try printf(out, " {f}\n", .{timer.untilNow(out, .real)});
    }

    return .{ E, C_ci };
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
