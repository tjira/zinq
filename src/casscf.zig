//! Solves the Complete Active Space Self-Consistent Field equations to optimize orbitals and CI states.

const std = @import("std");

const cblas = @import("cimport.zig").cblas;

const Allocator = std.mem.Allocator;

const FrequencyOptions = @import("frequency_analysis.zig").Options;
const HartreeFockOptions = @import("hartree_fock.zig").Options;
const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const OptimizeOptions = @import("molecular_optimization.zig").OptimizeOptions(GradientOptions);
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

const ao2mo_1e = @import("integral_transform.zig").ao2mo_1e;
const ao2mo_2e = @import("integral_transform.zig").ao2mo_2e;
const ao2mo_3c = @import("integral_transform.zig").ao2mo_3c;
const applyCayleyTransform = @import("casscf_orbital.zig").applyCayleyTransform;
const buildDiagonalFock = @import("casscf_orbital.zig").buildDiagonalFock;
const buildDiagonalFockRi = @import("casscf_orbital.zig").buildDiagonalFockRi;
const buildFittingIntegrals = @import("density_fitting.zig").buildFittingIntegrals;
const buildGeneralizedFock = @import("casscf_orbital.zig").buildGeneralizedFock;
const buildGeneralizedFockRi = @import("casscf_orbital.zig").buildGeneralizedFockRi;
const buildOrbitalRotationGenerator = @import("casscf_orbital.zig").buildOrbitalRotationGenerator;
const buildStateAveragedActRdm = @import("reduced_density_matrix.zig").buildStateAveragedActRdm;
const calculateCanonicalCasscfGradient = @import("casscf_gradient.zig").calculateCanonicalCasscfGradient;
const calculateRiCasscfGradient = @import("casscf_gradient.zig").calculateRiCasscfGradient;
const canonicalizeOrbitals = @import("casscf_orbital.zig").canonicalizeOrbitals;
const canonicalizeOrbitalsRi = @import("casscf_orbital.zig").canonicalizeOrbitalsRi;
const computeCi = @import("configuration_interaction.zig").computeCi;
const computeCiRi = @import("configuration_interaction.zig").computeCiRi;
const computeOrbitalGradient = @import("casscf_orbital.zig").computeOrbitalGradient;
const gemm = @import("linear_algebra.zig").gemm;
const generateCasDets = @import("determinant.zig").generateCasDets;
const printDeterminants = @import("determinant.zig").printDeterminants;
const printf = @import("read_write.zig").printf;
const syevd = @import("linear_algebra.zig").syevd;

/// Configurations for CASSCF active space dimensions, iterations, and derivative calculations.
pub const Options = struct {
    active_electrons: u32,
    active_orbitals: u32,
    aux_basis: ?[]const u8 = null,
    frequency: ?FrequencyOptions = null,
    gradient: ?GradientOptions = null,
    hartree_fock: HartreeFockOptions,
    hessian: ?HessianOptions = null,
    iterations: u32 = 100,
    multiplicity: u32 = 1,
    nacv: ?NacvOptions = null,
    nstate: u32 = 1,
    nthreads: u32 = 1,
    optimize: ?OptimizeOptions = null,
    orbital_threshold: f64 = 1e-6,
    spin_penalty: ?f64 = null,
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
pub fn solve(comptime T: type, io: std.Io, opt: Options, sys: *MolecularSystem(T), ref: HartreeFockResult(T), guess: ?Matrix(T), log: bool, gpa: Allocator) !Result(T) {
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

    var dets = try generateCasDets(ncore, nact, nact_el, opt.multiplicity, gpa);

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

    if (guess) |g| {
        try projectOrbitals(T, &C, g, hfres.ints.S.?, gpa);
    }

    if (guess == null) {
        @memcpy(C.data, hfres.C_mo.data);
    }

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

    if (opt.aux_basis) |aux| {
        const analytic_grad = if (opt.gradient) |g| g == .analytic else false;

        try buildFittingIntegrals(T, io, &hfres.ints, aux, analytic_grad, opt.nthreads, gpa);
    }

    if (opt.aux_basis != null) {
        try evaluateRiCasscf(T, io, &C, hfres.ints, dets.items, weights, opt, log, gpa);

        try canonicalizeOrbitalsRi(T, &C, hfres.ints, ncore, nact, gpa);
    }

    if (opt.aux_basis == null) {
        try evaluateCanonicalCasscf(T, io, &C, hfres.ints, dets.items, weights, opt, log, gpa);

        try canonicalizeOrbitals(T, &C, hfres.ints, ncore, nact, gpa);
    }

    const nrep = try hfres.ints.sys.nrep();

    const s_target = 0.5 * @as(f64, @floatFromInt(opt.multiplicity - 1));

    const s2_target = s_target * (s_target + 1.0);

    var final_energies: Vector(T), var final_C_ci: Matrix(T) = .{ undefined, undefined };

    if (opt.aux_basis != null) {
        var B_mo = try Tensor(T, 3).initZero(.{ nmo, nmo, hfres.ints.B.?.shape[2] }, gpa);
        defer B_mo.deinit(gpa);

        try ao2mo_3c(T, &B_mo, hfres.ints.B.?, C, gpa);

        final_energies, final_C_ci = try computeCiRi(T, dets.items, C, hfres.ints.H.?, B_mo, nrep, opt.spin_penalty, s2_target, ncore, nact, gpa);
    }

    if (opt.aux_basis == null) {
        final_energies, final_C_ci = try computeCi(T, null, C, hfres.ints, dets.items, nrep, false, opt.spin_penalty, s2_target, gpa);
    }

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

        try printDeterminants(T, io, final_C_ci, C, hfres.ints.S.?, dets.items, n_states, ncore, nact, false, gpa);
    }

    var grad: []Matrix(T) = &.{};

    errdefer {
        for (0..grad.len) |k| grad[k].deinit(gpa);

        if (grad.len > 0) gpa.free(grad);
    }

    if (opt.gradient) |g| switch (g) {
        .analytic => |a| {
            const states = a.states orelse &[_]u32{0};

            grad = if (opt.aux_basis != null)
                try calculateRiCasscfGradient(T, hfres, C, final_C_ci, dets.items, weights, ncore, nact, states, gpa)
            else
                try calculateCanonicalCasscfGradient(T, hfres, C, final_C_ci, dets.items, weights, ncore, nact, states, gpa);
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

/// Executes the two-step CASSCF macro-iteration cycle alternating configuration and orbital optimizations.
fn evaluateCanonicalCasscf(comptime T: type, io: std.Io, C: *Matrix(T), ints: anytype, dets: []const []const usize, weights: []const T, opt: Options, log: bool, gpa: Allocator) !void {
    const ncore, const nact = .{ (ints.sys.nel - opt.active_electrons) / 2, opt.active_orbitals };

    const nmo, var prev_energy: T = .{ C.ncol(), 0 };

    const s_target = 0.5 * @as(f64, @floatFromInt(opt.multiplicity - 1));

    const s2_target = s_target * (s_target + 1.0);

    for (0..opt.iterations) |iter| {
        var iter_timer = std.Io.Timestamp.now(io, .real);

        var all_E, var C_ci = try computeCi(T, null, C.*, ints, dets, try ints.sys.nrep(), false, opt.spin_penalty, s2_target, gpa);

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

/// Executes two-step CASSCF macro-iterations using density fitting for integral representations.
fn evaluateRiCasscf(comptime T: type, io: std.Io, C: *Matrix(T), ints: anytype, dets: []const []const usize, weights: []const T, opt: Options, log: bool, gpa: Allocator) !void {
    const ncore, const nact = .{ (ints.sys.nel - opt.active_electrons) / 2, opt.active_orbitals };

    const nmo, var prev_energy: T = .{ C.ncol(), 0 };

    const s_target = 0.5 * @as(f64, @floatFromInt(opt.multiplicity - 1));

    const s2_target = s_target * (s_target + 1.0);

    for (0..opt.iterations) |iter| {
        var iter_timer = std.Io.Timestamp.now(io, .real);

        var B_mo = try Tensor(T, 3).initZero(.{ nmo, nmo, ints.B.?.shape[2] }, gpa);
        defer B_mo.deinit(gpa);

        try ao2mo_3c(T, &B_mo, ints.B.?, C.*, gpa);

        var all_E, var C_ci = try computeCiRi(T, dets, C.*, ints.H.?, B_mo, try ints.sys.nrep(), opt.spin_penalty, s2_target, ncore, nact, gpa);

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

        var F_gen = try Matrix(T).initZero(nmo, nmo, gpa);
        defer F_gen.deinit(gpa);

        try buildGeneralizedFockRi(T, &F_gen, h_mo, B_mo, D_act, d_act, ncore, gpa);

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

            try buildDiagonalFockRi(T, &F_diag, h_mo, B_mo, D_act, ncore, gpa);

            var kappa = try Matrix(T).initZero(nmo, nmo, gpa);
            defer kappa.deinit(gpa);

            buildOrbitalRotationGenerator(T, &kappa, g_orb, F_diag, D_act, ncore);

            try applyCayleyTransform(T, C, kappa, gpa);
        }
    } else return error.CasscfDidNotConverge;
}

/// Projects and Lowdin-orthogonalizes an orbital guess against the metric overlap matrix S.
fn projectOrbitals(comptime T: type, C: *Matrix(T), guess: Matrix(T), S: Matrix(T), gpa: Allocator) !void {
    const nbf, const nmo = .{ C.nrow(), C.ncol() };

    var tmp = try Matrix(T).init(nbf, nmo, gpa);
    defer tmp.deinit(gpa);

    gemm(T, &tmp, S, guess, 1, 0, false, false);

    var S_mo = try Matrix(T).init(nmo, nmo, gpa);
    defer S_mo.deinit(gpa);

    gemm(T, &S_mo, guess, tmp, 1, 0, true, false);

    var eigvals = try Vector(T).init(nmo, gpa);
    defer eigvals.deinit(gpa);

    var U = try Matrix(T).init(nmo, nmo, gpa);
    defer U.deinit(gpa);

    try syevd(T, &eigvals, &U, S_mo);

    var U_scaled = try Matrix(T).init(nmo, nmo, gpa);
    defer U_scaled.deinit(gpa);

    for (0..nmo) |i| for (0..nmo) |j| {
        const ev = @max(1e-12, eigvals.at(j));

        U_scaled.ptr(i, j).* = U.at(i, j) / @sqrt(ev);
    };

    var S_inv_half = try Matrix(T).init(nmo, nmo, gpa);
    defer S_inv_half.deinit(gpa);

    gemm(T, &S_inv_half, U_scaled, U, 1, 0, false, true);

    gemm(T, C, guess, S_inv_half, 1, 0, false, false);
}
