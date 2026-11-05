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
const Vector = @import("tensor.zig").Vector;

const buildHamiltonianMatrix = @import("slater_condon.zig").buildHamiltonianMatrix;
const generateCasDets = @import("determinant.zig").generateCasDets;
const generateDets = @import("determinant.zig").generateDets;
const printf = @import("read_write.zig").printf;
const slater = @import("slater_condon.zig").slater;
const syevd = @import("linear_algebra.zig").syevd;
const transformInts = @import("integral_transform.zig").transformInts;

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

/// Evaluates Hamiltonian matrix elements in the determinant basis and diagonalizes to get CI state energies.
pub fn computeCi(comptime T: type, io: ?std.Io, C: Matrix(T), ints: anytype, dets: []const []const usize, VN: T, generalized: bool, gpa: Allocator) !struct { Vector(T), Matrix(T) } {
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

    if (io) |out| {
        try printf(out, " {f}\nCI DIAGONALIZATION TIME:", .{timer.untilNow(out, .real)});
        timer = std.Io.Timestamp.now(out, .real);
    }

    var E, var C_ci = try solveEigenvalueProblem(T, H_ci, VN, gpa);

    errdefer {
        E.deinit(gpa);

        C_ci.deinit(gpa);
    }

    if (io) |out| {
        try printf(out, " {f}\n", .{timer.untilNow(out, .real)});
    }

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

            for (0..n_states) |k| {
                try printf(io, "STATE {d:2} FINAL CI ENERGY: {d:.14} Eh\n", .{ k, E.at(k) });
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

/// Computes analytic nuclear gradients of CI state energies using dual number differentiation.
fn calculateNuclearGradient(comptime T: type, hfres: HartreeFockResult(T), C_ci: Matrix(T), dets: std.ArrayList([]const usize), states: []const u32, gpa: Allocator) ![]Matrix(T) {
    const generalized = hfres.ints.sys.nbf != hfres.C_mo.shape[0];

    var grads = try gpa.alloc(Matrix(T), states.len);
    errdefer gpa.free(grads);

    var grads_done: usize = 0;

    errdefer for (0..grads_done) |k| grads[k].deinit(gpa);

    for (0..states.len) |idx| {
        grads[idx] = try hfres.ints.sys.nuclearRepulsionGradient(gpa);

        grads_done += 1;
    }

    var H_d = try Matrix(ScalarDual(T)).init(hfres.ints.H.?.nrow(), hfres.ints.H.?.ncol(), gpa);
    defer H_d.deinit(gpa);

    var g_d = try Tensor(ScalarDual(T), 4).init(hfres.ints.g.?.shape, gpa);
    defer g_d.deinit(gpa);

    var C_d = try Matrix(ScalarDual(T)).init(hfres.C_mo.nrow(), hfres.C_mo.ncol(), gpa);
    defer C_d.deinit(gpa);

    for (0..3 * hfres.ints.sys.atoms.len) |coord| {
        for (0..H_d.nrow()) |mu| for (0..H_d.ncol()) |nu| {
            H_d.ptr(mu, nu).* = ScalarDual(T).init(hfres.ints.H.?.at(mu, nu), hfres.ints.dH.?.at(.{ coord, mu, nu }));
        };

        const nao = g_d.shape[0];

        for (0..nao) |mu| for (0..nao) |lambda| for (0..nao) |nu| for (0..nao) |sigma| {
            const indices = .{ mu, lambda, nu, sigma };

            const val = hfres.ints.g.?.at(indices);

            g_d.ptr(indices).* = ScalarDual(T).init(val, hfres.ints.dg.?.at(.{ coord, mu, lambda, nu, sigma }));
        };

        for (0..C_d.nrow()) |mu| for (0..C_d.ncol()) |p| {
            C_d.ptr(mu, p).* = ScalarDual(T).init(hfres.C_mo.at(mu, p), hfres.dC_mo.?.at(.{ coord, mu, p }));
        };

        var h_ms_d, var g_ms_d = try transformInts(ScalarDual(T), C_d, H_d, g_d, generalized, gpa);

        defer {
            h_ms_d.deinit(gpa);
            g_ms_d.deinit(gpa);
        }

        for (0..dets.items.len) |m| for (m..dets.items.len) |n| {
            const factor: T = if (m == n) 1 else 2;

            const val = slater(ScalarDual(T), dets.items[m], dets.items[n], h_ms_d, g_ms_d);

            for (states, 0..) |k, k_idx| {
                grads[k_idx].ptr(coord / 3, coord % 3).* += factor * C_ci.at(m, k) * C_ci.at(n, k) * val.der;
            }
        };
    }

    return grads;
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