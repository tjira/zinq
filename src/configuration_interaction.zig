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
const bfgs = @import("molecular_optimization.zig").bfgs;
const calculateHarmonicFrequencies = @import("frequency_analysis.zig").calculateHarmonicFrequencies;
const calculateNumericalGradient = @import("nuclear_derivative.zig").calculateNumericalGradient;
const calculateNumericalHessian = @import("nuclear_derivative.zig").calculateNumericalHessian;
const calculateNumericalNacv = @import("nuclear_derivative.zig").calculateNumericalNacv;
const calculateThermochemistry = @import("frequency_analysis.zig").calculateThermochemistry;
const DerivativeContext = @import("nuclear_derivative.zig").DerivativeContext;
const exportIfBuiltin = @import("molecular_integrals.zig").exportIfBuiltin;
const hartree_fock_runFromSystem = @import("hartree_fock.zig").runFromSystem;
const nuclearRepulsionGradient = @import("hartree_fock.zig").nuclearRepulsionGradient;
const primType = @import("value.zig").primType;
const printHarmonicFrequencies = @import("frequency_analysis.zig").printHarmonicFrequencies;
const printThermochemistry = @import("frequency_analysis.zig").printThermochemistry;
const printf = @import("read_write.zig").printf;
const steepestDescent = @import("molecular_optimization.zig").steepestDescent;
const syevd = @import("linear_algebra.zig").syevd;
const writeMatrix = @import("read_write.zig").writeMatrix;
const writeXyzFile = @import("read_write.zig").writeXyzFile;

/// Configurations for the CI calculation, excitation levels, optimization, and derivative settings.
pub const Options = struct {
    excitations: []const u32 = &.{ 1, 2 },
    frequency: ?FrequencyOptions = null,
    gradient: ?GradientOptions = null,
    hartree_fock: HartreeFockOptions,
    hessian: ?HessianOptions = null,
    nacv: ?NacvOptions = null,
    nstate: u32 = 1,
    nthreads: u32 = 1,
    optimize: ?OptimizeOptions = null,
    write: Write = .{},
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

        C: Matrix(T),

        energy: []T,

        grad: []Matrix(T) = &.{},
        hess: []Matrix(T) = &.{},
        nacv: []Matrix(T) = &.{},

        /// Frees all allocated memory associated with the CI results.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.hartree_fock.deinit(gpa);

            self.C.deinit(gpa);

            gpa.free(self.energy);

            for (0..self.grad.len) |i| {
                self.grad[i].deinit(gpa);
            }

            gpa.free(self.grad);

            for (0..self.hess.len) |i| {
                self.hess[i].deinit(gpa);
            }

            gpa.free(self.hess);

            for (0..self.nacv.len) |i| {
                self.nacv[i].deinit(gpa);
            }

            gpa.free(self.nacv);
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

/// Performs a CI calculation on a molecular system specified by file paths.
pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator) !Result(T) {
    try checkInvalidInput(opt);

    cblas.openblas_set_num_threads(@intCast(opt.nthreads));

    const basis_path = try exportIfBuiltin(io, opt.hartree_fock.basis, gpa);

    defer if (std.mem.startsWith(u8, opt.hartree_fock.basis, "builtin:")) {
        std.Io.Dir.cwd().deleteFile(io, basis_path) catch {};
    };

    const charge, const multi = .{ opt.hartree_fock.charge, opt.hartree_fock.multiplicity };

    var sys = try MolecularSystem(T).init(io, opt.hartree_fock.system, basis_path, charge, multi, gpa);
    defer sys.deinit(gpa);

    if (std.mem.startsWith(u8, opt.hartree_fock.basis, "builtin:")) {
        try std.Io.Dir.cwd().deleteFile(io, basis_path);
    }

    return try runFromSystem(T, io, opt, &sys, null, log, gpa);
}

/// Runs a CI calculation starting from an initialized MolecularSystem structure.
pub fn runFromSystem(comptime T: type, io: std.Io, opt: Options, sys: *MolecularSystem(T), Pg: ?Matrix(T), log: bool, gpa: Allocator) !Result(T) {
    try checkInvalidInput(opt);

    var final_Pg: ?Matrix(T) = null;
    defer if (final_Pg) |*p| p.deinit(gpa);

    if (opt.optimize) |o| {
        switch (o) {
            .steepest_descent => final_Pg = try steepestDescent(T, io, runFromSystem, opt, sys, Pg, log, gpa),
            .bfgs => final_Pg = try bfgs(T, io, runFromSystem, opt, sys, Pg, log, gpa),
        }
    }

    const generalized, var hf_opt = .{ opt.hartree_fock.generalized, opt.hartree_fock };

    if (opt.gradient != null and opt.gradient.? == .analytic) {
        hf_opt.gradient = .{ .analytic = .{} };

        if (hf_opt.response == null) {
            hf_opt.response = .{};
        }
    }

    var timer = std.Io.Timestamp.now(io, .real);

    var hfres = try hartree_fock_runFromSystem(T, io, hf_opt, sys, final_Pg orelse Pg, log, gpa);
    errdefer hfres.deinit(gpa);

    cblas.openblas_set_num_threads(@intCast(opt.nthreads));

    if (log) {
        try printf(io, "\nHARTREE-FOCK CYCLE TIME: {f}\n\nCI EXCITATIONS CONSIDERED: [", .{timer.untilNow(io, .real)});

        for (0..opt.excitations.len) |i| {
            try printf(io, "{d}{s}", .{ opt.excitations[i], if (i == opt.excitations.len - 1) "]" else ", " });
        }

        try printf(io, "\n\nSLATER DETERMINANT TIME:", .{});
    }

    timer = std.Io.Timestamp.now(io, .real);

    var dets = try generateDets(hfres.ints.sys.nel, 2 * hfres.ints.sys.nbf, opt.excitations, gpa);

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

    var E, var C = blk: {
        const res = try computeCiStates(T, if (log) io else null, generalized, hfres, dets.items, gpa);

        break :blk .{ res[0], res[1] };
    };

    errdefer {
        E.deinit(gpa);
        C.deinit(gpa);
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

    const grad_states = if (opt.gradient) |gradopt| switch (gradopt) {
        inline else => |g| g.states orelse &[_]u32{0},
    } else &[_]u32{};

    const ctx: DerivativeContext(T) = .{ .gpa = gpa, .io = io, .log = log, .sys = sys, .Pg = hfres.P };

    var grad: []Matrix(T) = undefined;

    if (opt.gradient == null) {
        grad = try gpa.alloc(Matrix(T), 0);
    }

    if (opt.gradient) |gradopt| switch (gradopt) {
        .analytic => grad = try gradient(T, hfres, C, dets, grad_states, gpa),
        .numeric => grad = try calculateNumericalGradient(T, ctx, runFromSystem, opt, grad_states),
    };

    errdefer {
        for (0..grad.len) |k| grad[k].deinit(gpa);

        gpa.free(grad);
    }

    if (log and opt.gradient != null) {
        const grad_type_str = if (opt.gradient.? == .analytic) "ANALYTIC" else "NUMERIC";

        for (grad_states, 0..) |state, i| {
            try printf(io, "\nCI STATE {d} {s} NUCLEAR ENERGY GRADIENT (Eh/a0)\n", .{ state, grad_type_str });

            for (0..grad[i].nrow()) |j| for (0..grad[i].ncol()) |k| {
                try printf(io, "{d:20.14}{s}", .{ grad[i].at(j, k), if (k == 2) "\n" else " " });
            };
        }
    }

    const hess_states = if (opt.hessian) |hessopt| switch (hessopt) {
        inline else => |h| h.states orelse &[_]u32{0},
    } else &[_]u32{};

    const hess = try handleHessianAndFrequencies(T, io, opt, sys, E.data, log, gpa);

    errdefer {
        for (0..hess.len) |k| hess[k].deinit(gpa);

        gpa.free(hess);
    }

    const nacv_pairs = if (opt.nacv) |nacvopt| switch (nacvopt) {
        inline else => |n| n.pairs orelse &[_][2]u32{.{ 0, 1 }},
    } else &[_][2]u32{};

    var nacv: []Matrix(T) = undefined;

    if (opt.nacv == null) {
        nacv = try gpa.alloc(Matrix(T), 0);
    }

    if (opt.nacv) |nacvopt| switch (nacvopt) {
        .numeric => nacv = try calculateNumericalNacv(T, ctx, runFromSystem, opt, nacv_pairs),
    };

    errdefer {
        for (0..nacv.len) |k| nacv[k].deinit(gpa);

        gpa.free(nacv);
    }

    if (log and opt.nacv != null) {
        for (nacv_pairs, 0..) |pair, idx| {
            const fmt = "\nCI NONADIABATIC COUPLING VECTOR BETWEEN STATES {d} AND {d} (1/a0)\n";

            try printf(io, fmt, .{ pair[0], pair[1] });

            for (0..nacv[idx].nrow()) |j| for (0..nacv[idx].ncol()) |k| {
                try printf(io, "{d:20.14}{s}", .{ nacv[idx].at(j, k), if (k == 2) "\n" else " " });
            };
        }
    }

    if (opt.write.gradient) |fname| {
        if (grad.len == 1) {
            try writeMatrix(T, io, fname, grad[0]);
        }

        if (grad.len > 1) for (grad, 0..) |g, idx| {
            const out_name = try getStatePath(fname, grad_states[idx], gpa);
            defer gpa.free(out_name);

            try writeMatrix(T, io, out_name, g);
        };
    }

    if (opt.write.hessian) |fname| {
        if (hess.len == 1) {
            try writeMatrix(T, io, fname, hess[0]);
        }

        if (hess.len > 1) for (hess, 0..) |h, idx| {
            const out_name = try getStatePath(fname, hess_states[idx], gpa);
            defer gpa.free(out_name);

            try writeMatrix(T, io, out_name, h);
        };
    }

    if (opt.write.nacv) |fname| {
        if (nacv.len == 1) {
            try writeMatrix(T, io, fname, nacv[0]);
        }

        if (nacv.len > 1) for (nacv, 0..) |n, idx| {
            const p0, const p1 = .{ nacv_pairs[idx][0], nacv_pairs[idx][1] };

            const out_name = try getNacvPath(fname, p0, p1, gpa);
            defer gpa.free(out_name);

            try writeMatrix(T, io, out_name, n);
        };
    }

    if (opt.write.geometry) |fname| {
        try writeXyzFile(T, io, fname, sys.atoms, sys.coors);
    }

    return Result(T){ .hartree_fock = hfres, .energy = E.data, .C = C, .grad = grad, .hess = hess, .nacv = nacv };
}

/// Computes the Hamiltonian matrix element between two Slater determinants using the Slater-Condon rules.
pub fn slater(comptime T: type, A: []const usize, B: []const usize, H_MS: Matrix(T), g_MS: Tensor(T, 4)) T {
    const diff = findDeterminantDifferences(A, B) orelse return Value(T).fromFloat(0).val;

    const s_A = signOfExcitations(A, diff.a[0..diff.len]);
    const s_B = signOfExcitations(B, diff.b[0..diff.len]);

    const phase = @as(primType(T), @floatFromInt(s_A * s_B));

    return evaluateSlaterCondon(T, A, diff.a, diff.b, diff.len, phase, H_MS, g_MS);
}

/// Validates input options for physical consistency and method requirements.
fn checkInvalidInput(opt: Options) !void {
    if (opt.hartree_fock.integral_direct) {
        std.log.err("INTEGRAL DIRECT CALCULATION IS NOT SUPPORTED FOR CONFIGURATION INTERACTION", .{});

        return error.InvalidInput;
    }

    if (opt.write.gradient != null and opt.gradient == null) {
        std.log.err("GRADIENT WRITE REQUESTED BUT GRADIENT IS NOT CALCULATED", .{});

        return error.InvalidInput;
    }

    if (opt.write.hessian != null and opt.hessian == null) {
        std.log.err("HESSIAN WRITE REQUESTED BUT HESSIAN IS NOT CALCULATED", .{});

        return error.InvalidInput;
    }

    if (opt.write.nacv != null and opt.nacv == null) {
        std.log.err("NACV WRITE REQUESTED BUT NACV IS NOT CALCULATED", .{});

        return error.InvalidInput;
    }

    if (opt.frequency != null and opt.hessian == null) {
        std.log.err("FREQUENCY CALCULATION REQUESTED BUT HESSIAN IS NOT CALCULATED", .{});

        return error.InvalidInput;
    }

    if (opt.excitations.len == 0) {
        std.log.err("CI EXCITATIONS LIST MUST NOT BE EMPTY", .{});

        return error.InvalidInput;
    }

    for (opt.excitations) |k| if (k == 0) {
        std.log.err("CI EXCITATION LEVEL MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    };

    if (opt.nthreads == 0) {
        std.log.err("THREAD COUNT MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.nstate == 0) {
        std.log.err("NUMBER OF CI STATES MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.gradient) |gradopt| {
        const states = switch (gradopt) {
            inline else => |g| g.states orelse &[_]u32{0},
        };

        if (states.len == 0) {
            std.log.err("CI GRADIENT STATES LIST MUST NOT BE EMPTY", .{});

            return error.InvalidInput;
        }

        for (states) |s| if (s >= opt.nstate) {
            std.log.err("GRADIENT STATE MUST BE LESS THAN NUMBER OF STATES", .{});

            return error.InvalidInput;
        };
    }

    if (opt.hessian) |hessopt| {
        const states = switch (hessopt) {
            inline else => |h| h.states orelse &[_]u32{0},
        };

        if (states.len == 0) {
            std.log.err("CI HESSIAN STATES LIST MUST NOT BE EMPTY", .{});

            return error.InvalidInput;
        }

        for (states) |s| if (s >= opt.nstate) {
            std.log.err("HESSIAN STATE MUST BE LESS THAN NUMBER OF STATES", .{});

            return error.InvalidInput;
        };
    }

    if (opt.nacv) |nacvopt| {
        if (opt.nstate < 2) {
            std.log.err("AT LEAST 2 STATES ARE REQUIRED FOR NACV CALCULATION", .{});

            return error.InvalidInput;
        }

        const pairs = switch (nacvopt) {
            inline else => |n| n.pairs orelse &[_][2]u32{.{ 0, 1 }},
        };

        if (pairs.len == 0) {
            std.log.err("CI NACV PAIRS LIST MUST NOT BE EMPTY", .{});

            return error.InvalidInput;
        }

        for (pairs) |pair| {
            if (pair[0] == pair[1]) {
                std.log.err("NACV STATES MUST BE DISTINCT", .{});

                return error.InvalidInput;
            }

            if (pair[0] >= opt.nstate or pair[1] >= opt.nstate) {
                std.log.err("NACV STATES MUST BE LESS THAN NUMBER OF STATES", .{});

                return error.InvalidInput;
            }
        }
    }

    if (opt.optimize) |optopt| {
        const gradopt = switch (optopt) {
            .bfgs => |b| b.gradient,
            .steepest_descent => |s| s.gradient,
        };

        const states = switch (gradopt) {
            inline else => |g| g.states orelse &[_]u32{0},
        };

        if (states.len != 1) {
            std.log.err("EXACTLY ONE STATE MUST BE SPECIFIED FOR GEOMETRY OPTIMIZATION", .{});

            return error.InvalidInput;
        }

        if (states[0] >= opt.nstate) {
            std.log.err("OPTIMIZATION STATE MUST BE LESS THAN NUMBER OF STATES", .{});

            return error.InvalidInput;
        }
    }
}

/// Evaluates Hamiltonian matrix elements in the determinant basis and diagonalizes to get CI state energies.
fn computeCiStates(comptime T: type, io: ?std.Io, generalized: bool, hfres: HartreeFockResult(T), dets: []const []const usize, gpa: Allocator) !struct { Vector(T), Matrix(T) } {
    var timer: std.Io.Timestamp = if (io) |out| std.Io.Timestamp.now(out, .real) else undefined;

    var H_MS, var g_MS = try transformInts(T, hfres.C, hfres.ints.H.?, hfres.ints.g.?, generalized, gpa);

    defer {
        H_MS.deinit(gpa);
        g_MS.deinit(gpa);
    }

    if (io) |out| {
        try printf(out, " {f}\nCI MATRIX ASSEMBLY TIME:", .{timer.untilNow(out, .real)});
    }

    if (io) |out| {
        timer = std.Io.Timestamp.now(out, .real);
    }

    var H_CI = try Matrix(T).init(dets.len, dets.len, gpa);
    defer H_CI.deinit(gpa);

    for (0..dets.len) |i| for (i..dets.len) |j| {
        const val = slater(T, dets[i], dets[j], H_MS, g_MS);

        H_CI.ptr(i, j).* = val;
        H_CI.ptr(j, i).* = val;
    };

    if (io) |out| {
        try printf(out, " {f}\nCI DIAGONALIZATION TIME:", .{timer.untilNow(out, .real)});
    }

    if (io) |out| {
        timer = std.Io.Timestamp.now(out, .real);
    }

    var E, var C = try solveEigenvalueProblem(T, H_CI, try hfres.ints.sys.nrep(), gpa);

    errdefer {
        E.deinit(gpa);
        C.deinit(gpa);
    }

    if (io) |out| {
        try printf(out, " {f}\n", .{timer.untilNow(out, .real)});
    }

    return .{ E, C };
}

/// Evaluates Hamiltonian matrix elements between Slater determinants using the Slater-Condon excitation rules.
fn evaluateSlaterCondon(comptime T: type, A: []const usize, diff_a: [2]usize, diff_b: [2]usize, len: usize, phase: primType(T), H_MS: Matrix(T), g_MS: Tensor(T, 4)) T {
    switch (len) {
        0 => {
            var sum = Value(T).fromFloat(0);

            for (0..A.len) |i| {
                sum = sum.add(Value(T).init(H_MS.at(A[i], A[i])));
            }

            for (0..A.len) |i| for (i + 1..A.len) |j| {
                const idx_i, const idx_j = .{ A[i], A[j] };

                const term1 = Value(T).init(g_MS.at(.{ idx_i, idx_j, idx_i, idx_j }));
                const term2 = Value(T).init(g_MS.at(.{ idx_i, idx_j, idx_j, idx_i }));

                sum = sum.add(term1.sub(term2));
            };

            return sum.val;
        },
        1 => {
            const m, const p = .{ diff_a[0], diff_b[0] };

            var term = Value(T).init(H_MS.at(m, p));

            for (0..A.len) |i| if (A[i] != m) {
                const term1 = Value(T).init(g_MS.at(.{ m, A[i], p, A[i] }));
                const term2 = Value(T).init(g_MS.at(.{ m, A[i], A[i], p }));

                term = term.add(term1.sub(term2));
            };

            return term.muls(phase).val;
        },
        2 => {
            const m, const n = .{ diff_a[0], diff_a[1] };
            const p, const q = .{ diff_b[0], diff_b[1] };

            const term1 = Value(T).init(g_MS.at(.{ m, n, p, q }));
            const term2 = Value(T).init(g_MS.at(.{ m, n, q, p }));

            return term1.sub(term2).muls(phase).val;
        },
        else => return Value(T).fromFloat(0).val,
    }
}

/// Compares two Slater determinants and extracts orbital differences when excitation degree is at most two.
fn findDeterminantDifferences(A: []const usize, B: []const usize) ?struct { a: [2]usize, b: [2]usize, len: usize } {
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

/// Generates all k-combinations from a set of size n, representing orbital indices.
fn generateCombinations(n: usize, k: usize, offset: usize, gpa: Allocator) !std.ArrayList([]const usize) {
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

/// Derives output paths for nonadiabatic coupling vectors by inserting state pair indices before the file extension.
fn getNacvPath(path: []const u8, p0: u32, p1: u32, gpa: Allocator) ![]const u8 {
    const ext = std.fs.path.extension(path);

    return try std.fmt.allocPrint(gpa, "{s}_{d}_{d}{s}", .{ path[0 .. path.len - ext.len], p0, p1, ext });
}

/// Derives output paths for state nuclear derivatives by inserting the state index before the file extension.
fn getStatePath(path: []const u8, state: u32, gpa: Allocator) ![]const u8 {
    const ext = std.fs.path.extension(path);

    return try std.fmt.allocPrint(gpa, "{s}_{d}{s}", .{ path[0 .. path.len - ext.len], state, ext });
}

/// Computes analytic nuclear gradients of CI state energies using dual number differentiation.
fn gradient(comptime T: type, hfres: HartreeFockResult(T), C: Matrix(T), dets: std.ArrayList([]const usize), states: []const u32, gpa: Allocator) ![]Matrix(T) {
    const generalized = hfres.ints.sys.nbf != hfres.C.shape[0];

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

    var C_d = try Matrix(ScalarDual(T)).init(hfres.C.nrow(), hfres.C.ncol(), gpa);
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
            C_d.ptr(mu, p).* = ScalarDual(T).init(hfres.C.at(mu, p), hfres.dC.?.at(.{ i, mu, p }));
        };

        var H_MS_d, var g_MS_d = try transformInts(ScalarDual(T), C_d, H_d, g_d, generalized, gpa);

        defer {
            H_MS_d.deinit(gpa);
            g_MS_d.deinit(gpa);
        }

        var H_CI_d = try Matrix(ScalarDual(T)).init(dets.items.len, dets.items.len, gpa);
        defer H_CI_d.deinit(gpa);

        for (0..dets.items.len) |a| for (a..dets.items.len) |b| {
            const val = slater(ScalarDual(T), dets.items[a], dets.items[b], H_MS_d, g_MS_d);

            H_CI_d.ptr(a, b).* = val;
            H_CI_d.ptr(b, a).* = val;
        };

        for (0..dets.items.len) |a| for (0..dets.items.len) |b| {
            const der = H_CI_d.at(a, b).der;

            for (states, 0..) |s, s_idx| {
                grads[s_idx].ptr(i / 3, i % 3).* += C.at(a, s) * C.at(b, s) * der;
            }
        };
    }

    return grads;
}

/// Computes nuclear Hessians of CI states numerically and performs frequency and thermochemical analyses.
fn handleHessianAndFrequencies(comptime T: type, io: std.Io, opt: Options, sys: *MolecularSystem(T), energies: []const T, log: bool, gpa: Allocator) ![]Matrix(T) {
    if (opt.hessian == null) return try gpa.alloc(Matrix(T), 0);

    const states = switch (opt.hessian.?) {
        inline else => |h| h.states orelse &[_]u32{0},
    };

    const ctx: DerivativeContext(T) = .{ .gpa = gpa, .io = io, .log = log, .sys = sys };

    var hess = try calculateNumericalHessian(T, ctx, runFromSystem, opt, states);

    errdefer {
        for (0..hess.len) |k| hess[k].deinit(gpa);

        gpa.free(hess);
    }

    if (log and opt.frequency != null) {
        for (states, 0..) |state, idx| {
            var freqs = try calculateHarmonicFrequencies(T, hess[idx], sys.*, gpa);
            defer freqs.deinit(gpa);

            const method_str = try std.fmt.allocPrint(gpa, "CI STATE {d}", .{state});
            defer gpa.free(method_str);

            try printHarmonicFrequencies(T, io, freqs, method_str);

            const mult = opt.hartree_fock.multiplicity;

            const thermo = try calculateThermochemistry(T, opt.frequency.?, sys.*, freqs, mult, gpa);

            try printThermochemistry(T, io, opt.frequency.?, thermo, energies[state], method_str);
        }
    }

    return hess;
}

/// Computes the permutation sign (+1 or -1) associated with mapping two determinants via orbital excitations.
fn signOfExcitations(S: []const usize, U: []const usize) i2 {
    var exp: usize = 0;

    for (0..U.len) |j| {
        exp += std.mem.indexOfScalar(usize, S, U[j]).? - j;
    }

    return if (exp % 2 == 0) 1 else -1;
}

/// Diagonalizes the CI Hamiltonian matrix to yield state energies and coefficients.
fn solveEigenvalueProblem(comptime T: type, H_CI: Matrix(T), VN: T, gpa: Allocator) !struct { Vector(T), Matrix(T) } {
    var E = try Vector(T).init(H_CI.shape[0], gpa);
    errdefer E.deinit(gpa);

    var C = try Matrix(T).initZero(H_CI.shape[0], H_CI.shape[0], gpa);
    errdefer C.deinit(gpa);

    try syevd(T, &E, &C, H_CI);

    for (0..E.length()) |i| {
        E.ptr(i).* += VN;
    }

    return .{ E, C };
}

/// Transforms one- and two-electron integrals from the atomic orbital basis to the molecular spin-orbital basis.
fn transformInts(comptime T: type, C: Matrix(T), H: Matrix(T), g: Tensor(T, 4), generalized: bool, gpa: Allocator) !struct { Matrix(T), Tensor(T, 4) } {
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

    var H_MS = try Matrix(T).init(nsp, nsp, gpa);
    errdefer H_MS.deinit(gpa);

    ao2mo_pp(T, &H_MS, H_SO, C_SO);

    var g_MS = try Tensor(T, 4).init(.{ nsp, nsp, nsp, nsp }, gpa);
    errdefer g_MS.deinit(gpa);

    try ao2mo_pppp(T, &g_MS, g_SO, C_SO, gpa);

    return .{ H_MS, g_MS };
}
