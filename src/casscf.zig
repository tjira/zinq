//! Solves the Complete Active Space Self-Consistent Field equations to optimize orbitals and CI states.

const std = @import("std");

const cblas = @import("cimport.zig").cblas;

const Allocator = std.mem.Allocator;

const DerivativeContext = @import("nuclear_derivative.zig").DerivativeContext;
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
const bfgs = @import("molecular_optimization.zig").bfgs;
const buildActiveRdm = @import("configuration_interaction.zig").buildActiveRdm;
const calculateHarmonicFrequencies = @import("frequency_analysis.zig").calculateHarmonicFrequencies;
const calculateNumericalGradient = @import("nuclear_derivative.zig").calculateNumericalGradient;
const calculateNumericalHessian = @import("nuclear_derivative.zig").calculateNumericalHessian;
const calculateNumericalNacv = @import("nuclear_derivative.zig").calculateNumericalNacv;
const calculateThermochemistry = @import("frequency_analysis.zig").calculateThermochemistry;
const computeCi = @import("configuration_interaction.zig").computeCi;
const exportIfBuiltin = @import("molecular_integrals.zig").exportIfBuiltin;
const generateCasDets = @import("configuration_interaction.zig").generateCasDets;
const gemm = @import("linear_algebra.zig").gemm;
const getNacvPath = @import("configuration_interaction.zig").getNacvPath;
const getrf = @import("linear_algebra.zig").getrf;
const getrs = @import("linear_algebra.zig").getrs;
const getStatePath = @import("configuration_interaction.zig").getStatePath;
const hartree_fock_runFromSystem = @import("hartree_fock.zig").runFromSystem;
const nuclearRepulsionGradient = @import("hartree_fock.zig").nuclearRepulsionGradient;
const printHarmonicFrequencies = @import("frequency_analysis.zig").printHarmonicFrequencies;
const printThermochemistry = @import("frequency_analysis.zig").printThermochemistry;
const printf = @import("read_write.zig").printf;
const slater = @import("configuration_interaction.zig").slater;
const steepestDescent = @import("molecular_optimization.zig").steepestDescent;
const syevd = @import("linear_algebra.zig").syevd;
const transformInts = @import("configuration_interaction.zig").transformInts;
const writeMatrix = @import("read_write.zig").writeMatrix;
const writeXyzFile = @import("read_write.zig").writeXyzFile;

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
    threshold: f64 = 1e-8,
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

        C: Matrix(T),

        energy: []T,

        grad: []Matrix(T) = &.{},
        hess: []Matrix(T) = &.{},
        nacv: []Matrix(T) = &.{},

        C_ci: Matrix(T),

        /// Frees all allocated memory associated with the CASSCF results.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.hartree_fock.deinit(gpa);

            self.C.deinit(gpa);

            self.C_ci.deinit(gpa);

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

/// Performs a CASSCF calculation on a molecular system specified by file paths.
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

/// Executes the CASSCF optimization algorithm on a provided MolecularSystem instance.
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

    var timer = std.Io.Timestamp.now(io, .real);

    var hf_opt = opt.hartree_fock;

    if (opt.gradient != null and opt.gradient.? == .analytic) {
        hf_opt.gradient = .{ .analytic = .{} };
    }

    var hfres = try hartree_fock_runFromSystem(T, io, hf_opt, sys, final_Pg orelse Pg, log, gpa);
    errdefer hfres.deinit(gpa);

    cblas.openblas_set_num_threads(@intCast(opt.nthreads));

    const nel, const nact_el, const nact = .{ hfres.ints.sys.nel, opt.active_electrons, opt.active_orbitals };

    if (nact_el > nel or (nel - nact_el) % 2 != 0) {
        std.log.err("ACTIVE ELECTRONS MUST BE LESS THAN OR EQUAL TO TOTAL ELECTRONS WITH CLOSED SHELL CORE", .{});

        return error.InvalidInput;
    }

    const ncore, const nmo = .{ (nel - nact_el) / 2, hfres.C.ncol() };

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

    var C = try Matrix(T).init(hfres.C.nrow(), hfres.C.ncol(), gpa);
    errdefer C.deinit(gpa);

    @memcpy(C.data, hfres.C.data);

    var prev_energy: T = 0;

    for (0..opt.iterations) |iter| {
        var iter_timer = std.Io.Timestamp.now(io, .real);

        var all_E, var C_ci = try computeCi(T, null, C, hfres.ints, dets.items, try hfres.ints.sys.nrep(), false, gpa);

        defer {
            all_E.deinit(gpa);

            C_ci.deinit(gpa);
        }

        const ci_E = all_E.at(0);

        const e_diff = if (iter == 0) ci_E else @abs(ci_E - prev_energy);

        prev_energy = ci_E;

        var D_act, var d_act = try buildActiveRdm(T, dets.items, C_ci, 0, ncore, nact, gpa);

        defer {
            D_act.deinit(gpa);
            d_act.deinit(gpa);
        }

        var H_mo = try Matrix(T).init(nmo, nmo, gpa);
        defer H_mo.deinit(gpa);

        ao2mo_pp(T, &H_mo, hfres.ints.H.?, C);

        var g_mo = try Tensor(T, 4).init(.{ nmo, nmo, nmo, nmo }, gpa);
        defer g_mo.deinit(gpa);

        try ao2mo_pppp(T, &g_mo, hfres.ints.g.?, C, gpa);

        var F_gen = try Matrix(T).initZero(nmo, nmo, gpa);
        defer F_gen.deinit(gpa);

        var F_diag = try Vector(T).init(nmo, gpa);
        defer F_diag.deinit(gpa);

        try buildGeneralizedFock(T, &F_gen, &F_diag, H_mo, g_mo, D_act, d_act, ncore, nact, gpa);

        var g_orb = try Matrix(T).initZero(nmo, nmo, gpa);
        defer g_orb.deinit(gpa);

        const grad_norm = computeOrbitalGradient(T, &g_orb, F_gen, ncore, nact);

        if (log) {
            const fmt = "{d:4} {d:20.14} {e:9.3} {e:9.3} {s:9}\n";

            const time_str = try std.fmt.allocPrint(gpa, "{f}", .{iter_timer.untilNow(io, .real)});
            defer gpa.free(time_str);

            const de = if (iter == 0) 0 else e_diff;

            try printf(io, fmt, .{ iter + 1, ci_E, de, grad_norm, time_str });
        }

        if (iter > 0 and e_diff < opt.threshold and grad_norm < @max(10 * opt.threshold, 1e-5)) {
            break;
        }

        if (iter < opt.iterations - 1) {
            try computeOrbitalStep(T, &C, g_orb, F_diag, D_act, ncore, nact, gpa);
        }
    }

    try canonicalizeOrbitals(T, &C, hfres.ints, ncore, nact, gpa);

    const nrep = try hfres.ints.sys.nrep();

    var final_energies, var final_C_ci = try computeCi(T, null, C, hfres.ints, dets.items, nrep, false, gpa);

    errdefer {
        final_energies.deinit(gpa);

        final_C_ci.deinit(gpa);
    }
    const final_energy = final_energies.at(0);

    if (log) {
        const n_states = @min(opt.nstate, final_energies.length());

        if (n_states == 1) {
            try printf(io, "\nFINAL CASSCF ENERGY: {d:.14} Eh\n", .{final_energy});
        }

        if (n_states > 1) {
            try printf(io, "\n", .{});

            for (0..n_states) |s| {
                try printf(io, "STATE {d:2} FINAL CASSCF ENERGY: {d:.14} Eh\n", .{ s, final_energies.at(s) });
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
        .analytic => grad = try gradient(T, hfres, C, final_C_ci, dets.items, grad_states, gpa),
        .numeric => grad = try calculateNumericalGradient(T, ctx, runFromSystem, opt, grad_states),
    };

    errdefer {
        for (0..grad.len) |k| grad[k].deinit(gpa);

        gpa.free(grad);
    }

    if (log and opt.gradient != null) {
        const grad_type_str = if (opt.gradient.? == .analytic) "ANALYTIC" else "NUMERIC";

        for (grad_states, 0..) |s, idx| {
            if (grad_states.len > 1) {
                const fmt = "\nCASSCF {s} GRADIENT FOR STATE {d} (Eh/a0)\n";

                try printf(io, fmt, .{ grad_type_str, s });
            }

            if (grad_states.len == 1) {
                const fmt = "\nCASSCF {s} GRADIENT (Eh/a0)\n";

                try printf(io, fmt, .{grad_type_str});
            }

            for (0..grad[idx].nrow()) |j| for (0..grad[idx].ncol()) |k| {
                try printf(io, "{d:20.14}{s}", .{ grad[idx].at(j, k), if (k == 2) "\n" else " " });
            };
        }
    }

    const hess = try handleHessianAndFrequencies(T, io, opt, sys, final_energies.data, log, gpa);

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
            const fmt = "\nCASSCF NONADIABATIC COUPLING VECTOR BETWEEN STATES {d} AND {d} (1/a0)\n";
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
            const path = try getStatePath(fname, grad_states[idx], gpa);
            defer gpa.free(path);

            try writeMatrix(T, io, path, g);
        };
    }

    if (opt.write.nacv) |fname| {
        if (nacv.len == 1) {
            try writeMatrix(T, io, fname, nacv[0]);
        }

        if (nacv.len > 1) for (nacv, 0..) |n, idx| {
            const p0, const p1 = .{ nacv_pairs[idx][0], nacv_pairs[idx][1] };

            const path = try getNacvPath(fname, p0, p1, gpa);
            defer gpa.free(path);

            try writeMatrix(T, io, path, n);
        };
    }

    if (opt.write.orbitals) |fname| {
        try writeMatrix(T, io, fname, C);
    }

    if (opt.write.geometry) |fname| {
        try writeXyzFile(T, io, fname, sys.atoms, sys.coors);
    }

    return Result(T){
        .hartree_fock = hfres,
        .energy = final_energies.data,
        .C = C,
        .grad = grad,
        .hess = hess,
        .nacv = nacv,
        .C_ci = final_C_ci,
    };
}

/// Constructs the core Fock matrix by contracting one-electron integrals with closed-shell core orbitals.
fn buildCoreFock(comptime T: type, F_core: *Matrix(T), H_mo: Matrix(T), g_mo: Tensor(T, 4), ncore: usize) void {
    const nmo = H_mo.nrow();

    for (0..nmo) |p| for (0..nmo) |q| {
        var sum = H_mo.at(p, q);

        for (0..ncore) |k| {
            sum += 2 * g_mo.at(.{ p, k, q, k }) - g_mo.at(.{ p, k, k, q });
        }

        F_core.ptr(p, q).* = sum;
    };
}

/// Constructs the generalized Fock matrix and diagonal orbital energies using active reduced density matrices.
fn buildGeneralizedFock(comptime T: type, F_gen: *Matrix(T), F_diag: *Vector(T), H_mo: Matrix(T), g_mo: Tensor(T, 4), D_act: Matrix(T), d_act: Tensor(T, 4), ncore: usize, nact: usize, gpa: Allocator) !void {
    const nmo = H_mo.nrow();

    var F_core = try Matrix(T).init(nmo, nmo, gpa);
    defer F_core.deinit(gpa);

    buildCoreFock(T, &F_core, H_mo, g_mo, ncore);

    @memset(F_gen.data, 0);

    for (0..nmo) |p| {
        var sum_diag = F_core.at(p, p);

        for (0..nact) |u| for (0..nact) |v| {
            const U, const V = .{ ncore + u, ncore + v };

            sum_diag += D_act.at(u, v) * (g_mo.at(.{ p, U, p, V }) - 0.5 * g_mo.at(.{ p, U, V, p }));
        };

        F_diag.ptr(p).* = sum_diag;

        for (0..ncore) |q| {
            var sum = F_core.at(p, q);

            for (0..nact) |u| for (0..nact) |v| {
                const U, const V = .{ ncore + u, ncore + v };

                sum += D_act.at(u, v) * (g_mo.at(.{ p, U, q, V }) - 0.5 * g_mo.at(.{ p, U, V, q }));
            };

            F_gen.ptr(p, q).* = 2 * sum;
        }

        for (0..nact) |u| {
            const q, var sum: T = .{ ncore + u, 0 };

            for (0..nact) |v| {
                sum += F_core.at(p, ncore + v) * D_act.at(v, u);
            }

            for (0..nact) |v| for (0..nact) |w| for (0..nact) |x| {
                const V, const W, const X = .{ ncore + v, ncore + w, ncore + x };

                sum += g_mo.at(.{ p, W, V, X }) * d_act.at(.{ u, v, w, x });
            };

            F_gen.ptr(p, q).* = sum;
        }
    }
}

/// Canonicalizes core and virtual orbital blocks by diagonalizing their respective Fock matrices.
fn canonicalizeOrbitals(comptime T: type, C: *Matrix(T), ints: anytype, ncore: usize, nact: usize, gpa: Allocator) !void {
    const nmo = C.ncol();
    const nbf = C.nrow();

    var H_mo = try Matrix(T).init(nmo, nmo, gpa);
    defer H_mo.deinit(gpa);

    ao2mo_pp(T, &H_mo, ints.H.?, C.*);

    var g_mo = try Tensor(T, 4).init(.{ nmo, nmo, nmo, nmo }, gpa);
    defer g_mo.deinit(gpa);

    try ao2mo_pppp(T, &g_mo, ints.g.?, C.*, gpa);

    var F_core = try Matrix(T).init(nmo, nmo, gpa);
    defer F_core.deinit(gpa);

    buildCoreFock(T, &F_core, H_mo, g_mo, ncore);

    if (ncore > 1) {
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

    const nvirt = nmo - ncore - nact;

    if (nvirt > 1) {
        var F_vv = try Matrix(T).init(nvirt, nvirt, gpa);
        defer F_vv.deinit(gpa);

        const v_start = ncore + nact;

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
}

/// Validates input configuration for consistency, orbital ranges, and active space limits.
fn checkInvalidInput(opt: Options) !void {
    if (opt.hartree_fock.integral_direct) {
        std.log.err("INTEGRAL DIRECT CALCULATION IS NOT SUPPORTED FOR CASSCF", .{});

        return error.InvalidInput;
    }

    if (opt.active_electrons == 0 or opt.active_orbitals == 0) {
        std.log.err("CASSCF ACTIVE ELECTRONS AND ACTIVE ORBITALS MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.nthreads == 0 or opt.iterations == 0 or opt.threshold <= 0 or opt.nstate == 0) {
        std.log.err("INVALID THREADS, ITERATIONS, THRESHOLD, OR NSTATE", .{});

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

    if (opt.gradient) |gradopt| {
        const states = switch (gradopt) {
            inline else => |g| g.states orelse &[_]u32{0},
        };

        if (states.len == 0) {
            std.log.err("CASSCF GRADIENT STATES LIST MUST NOT BE EMPTY", .{});

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
            std.log.err("CASSCF HESSIAN STATES LIST MUST NOT BE EMPTY", .{});

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
            std.log.err("CASSCF NACV PAIRS LIST MUST NOT BE EMPTY", .{});

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
}

/// Evaluates the orbital gradient matrix and returns the maximum gradient norm among non-redundant pairs.
fn computeOrbitalGradient(comptime T: type, g_orb: *Matrix(T), F_gen: Matrix(T), ncore: usize, nact: usize) T {
    const nmo = F_gen.nrow();
    var max_norm: T = 0;

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

/// Updates molecular orbitals via Cayley transform using the orbital gradient.
fn computeOrbitalStep(comptime T: type, C: *Matrix(T), g_orb: Matrix(T), F_diag: Vector(T), D_act: Matrix(T), ncore: usize, nact: usize, gpa: Allocator) !void {
    const nmo = C.ncol();
    const nbf = C.nrow();

    var kappa = try Matrix(T).initZero(nmo, nmo, gpa);
    defer kappa.deinit(gpa);

    for (0..ncore) |j| for (ncore + nact..nmo) |b| {
        var denom = 4 * (F_diag.at(b) - F_diag.at(j));

        if (denom < 0.1) denom = 0.1;

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

            if (denom < 0.1) denom = 0.1;

            var step = -g_orb.at(b, u) / denom;

            step = std.math.clamp(step, -0.2, 0.2);

            kappa.ptr(b, u).* = 0 + step;
            kappa.ptr(u, b).* = 0 - step;
        }

        for (0..ncore) |j| {
            var denom = 2 * (2 - Duu) * (F_diag.at(u) - F_diag.at(j));

            if (denom < 0.1) denom = 0.1;

            var step = -g_orb.at(u, j) / denom;

            step = std.math.clamp(step, -0.2, 0.2);

            kappa.ptr(u, j).* = 0 + step;
            kappa.ptr(j, u).* = 0 - step;
        }
    }

    var A = try Matrix(T).initZero(nmo, nmo, gpa);
    defer A.deinit(gpa);

    var U = try Matrix(T).initZero(nmo, nmo, gpa);
    defer U.deinit(gpa);

    for (0..nmo) |i| for (0..nmo) |j| {
        const delta: T = if (i == j) 1 else 0;

        A.ptr(i, j).* = delta - 0.5 * kappa.at(i, j);
        U.ptr(i, j).* = delta + 0.5 * kappa.at(i, j);
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

/// Computes analytic nuclear gradients of CASSCF state energies using dual number differentiation.
fn gradient(comptime T: type, hfres: HartreeFockResult(T), C: Matrix(T), C_ci: Matrix(T), dets: []const []const usize, states: []const u32, gpa: Allocator) ![]Matrix(T) {
    var grads = try gpa.alloc(Matrix(T), states.len);
    errdefer gpa.free(grads);

    var grads_done: usize = 0;
    errdefer for (0..grads_done) |k| grads[k].deinit(gpa);

    for (0..states.len) |idx| {
        grads[idx] = try nuclearRepulsionGradient(T, hfres.ints.sys, gpa);

        grads_done += 1;
    }

    const nbf = C.nrow();
    const nmo = C.ncol();

    const natoms = hfres.ints.sys.atoms.len;

    var H_d = try Matrix(ScalarDual(T)).init(nbf, nbf, gpa);
    defer H_d.deinit(gpa);

    var g_d = try Tensor(ScalarDual(T), 4).init(.{ nbf, nbf, nbf, nbf }, gpa);
    defer g_d.deinit(gpa);

    var C_d = try Matrix(ScalarDual(T)).init(nbf, nmo, gpa);
    defer C_d.deinit(gpa);

    var dS_AO = try Matrix(T).init(nbf, nbf, gpa);
    defer dS_AO.deinit(gpa);

    var tmp = try Matrix(T).init(nbf, nmo, gpa);
    defer tmp.deinit(gpa);

    var S_mo = try Matrix(T).init(nmo, nmo, gpa);
    defer S_mo.deinit(gpa);

    var dC = try Matrix(T).init(nbf, nmo, gpa);
    defer dC.deinit(gpa);

    for (0..3 * natoms) |coord| {
        for (0..nbf) |mu| for (0..nbf) |nu| {
            dS_AO.ptr(mu, nu).* = hfres.ints.dS.?.at(.{ coord, mu, nu });
        };

        gemm(T, &tmp, dS_AO, C, 1, 0, false, false);

        gemm(T, &S_mo, C, tmp, 1, 0, true, false);

        gemm(T, &dC, C, S_mo, -0.5, 0, false, false);

        for (0..nbf) |mu| for (0..nmo) |p| {
            C_d.ptr(mu, p).* = ScalarDual(T).init(C.at(mu, p), dC.at(mu, p));
        };

        for (0..nbf) |mu| for (0..nbf) |nu| {
            H_d.ptr(mu, nu).* = ScalarDual(T).init(hfres.ints.H.?.at(mu, nu), hfres.ints.dH.?.at(.{ coord, mu, nu }));
        };

        const n = g_d.shape[0];

        for (0..n) |mu| for (0..n) |lambda| for (0..n) |nu| for (0..n) |sigma| {
            const indices = .{ mu, lambda, nu, sigma };

            const val = hfres.ints.g.?.at(indices);

            const der = hfres.ints.dg.?.at(.{ coord, mu, lambda, nu, sigma });

            g_d.ptr(indices).* = ScalarDual(T).init(val, der);
        };

        var H_MS_d, var g_MS_d = try transformInts(ScalarDual(T), C_d, H_d, g_d, false, gpa);

        defer {
            H_MS_d.deinit(gpa);
            g_MS_d.deinit(gpa);
        }

        for (0..dets.len) |d1| for (d1..dets.len) |d2| {
            const factor: T = if (d1 == d2) 1 else 2;

            const val = slater(ScalarDual(T), dets[d1], dets[d2], H_MS_d, g_MS_d);

            for (states, 0..) |s, s_idx| {
                grads[s_idx].ptr(coord / 3, coord % 3).* += factor * C_ci.at(d1, s) * C_ci.at(d2, s) * val.der;
            }
        };
    }

    return grads;
}

/// Computes numerical nuclear Hessians for CASSCF states and conducts frequency and thermochemical analysis.
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

            const method_str = try std.fmt.allocPrint(gpa, "CASSCF STATE {d}", .{state});
            defer gpa.free(method_str);

            try printHarmonicFrequencies(T, io, freqs, method_str);

            const mult = opt.hartree_fock.multiplicity;

            const thermo = try calculateThermochemistry(T, opt.frequency.?, sys.*, freqs, mult, gpa);

            try printThermochemistry(T, io, opt.frequency.?, thermo, energies[state], method_str);
        }
    }

    return hess;
}
