//! Solves the Hartree-Fock self-consistent field equations to approximate the electronic wavefunction of a system.

const std = @import("std");

const cblas = @import("cimport.zig").cblas;

const Allocator = std.mem.Allocator;

const DftPotential = @import("density_functional_theory.zig").DftPotential;
const FrequencyOptions = @import("frequency_analysis.zig").Options;
const Integrals = @import("molecular_integrals.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const MolecularIntegralsOptions = @import("molecular_integrals.zig").Options;
const MolecularIntegralsWrite = @import("molecular_integrals.zig").Write;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

const ao2mo_pp = @import("integral_transform.zig").ao2mo_pp;
const axpy = @import("linear_algebra.zig").axpy;
const calculateTotalSpin = @import("spin_analysis.zig").calculateTotalSpin;
const dot = @import("linear_algebra.zig").dot;
const gemm = @import("linear_algebra.zig").gemm;
const getrf = @import("linear_algebra.zig").getrf;
const getrs = @import("linear_algebra.zig").getrs;
const getSymbol = @import("constant.zig").getSymbol;
const lowdin = @import("population_analysis.zig").lowdin;
const mayer = @import("population_analysis.zig").mayer;
const mo2ao_xx = @import("integral_transform.zig").mo2ao_xx;
const molecular_integrals_run = @import("molecular_integrals.zig").run;
const molecular_integrals_runFromSystem = @import("molecular_integrals.zig").runFromSystem;
const mulliken = @import("population_analysis.zig").mulliken;
const orbitalResponse = @import("cphf.zig").orbitalResponse;
const printLowdinCharges = @import("population_analysis.zig").printLowdinCharges;
const printMayerBondOrders = @import("population_analysis.zig").printMayerBondOrders;
const printMullikenCharges = @import("population_analysis.zig").printMullikenCharges;
const printTotalSpin = @import("spin_analysis.zig").printTotalSpin;
const printWibergBondOrders = @import("population_analysis.zig").printWibergBondOrders;
const printf = @import("read_write.zig").printf;
const sygvd = @import("linear_algebra.zig").sygvd;
const wiberg = @import("population_analysis.zig").wiberg;
const writeMatrix = @import("read_write.zig").writeMatrix;

const AN2SM = @import("constant.zig").AN2SM;
const AU2CM = @import("constant.zig").AU2CM;

/// Parameters governing the self-consistent field (SCF) calculation convergence and method options.
pub const Options = struct {
    basis: []const u8,
    charge: i32 = 0,
    dft: ?DftOptions = null,
    diis: ?u32 = 8,
    frequency: ?FrequencyOptions = null,
    generalized: bool = false,
    gradient: ?GradientOptions = null,
    hessian: ?HessianOptions = null,
    integral_direct: bool = false,
    iterations: u32 = 100,
    lowdin: bool = false,
    mayer: bool = false,
    mulliken: bool = false,
    multiplicity: u32 = 1,
    nthreads: u32 = 1,
    optimize: ?OptimizeOptions = null,
    response: ?ResponseOptions = null,
    system: []const u8,
    threshold: f64 = 1e-8,
    wiberg: bool = false,
    write: Write = .{},
};

/// Options for computing the nuclear gradient analytically or numerically.
pub const GradientOptions = union(enum) {
    analytic: GradientAnalyticOptions,
    numeric: GradientNumericOptions,
};

/// Tagged union specifying nuclear Hessian calculation methods.
pub const HessianOptions = union(enum) {
    numeric: HessianNumericOptions,
};

/// Tagged union specifying geometry optimization algorithms.
pub const OptimizeOptions = union(enum) {
    bfgs: BfgsOptions,
    steepest_descent: SteepestDescentOptions,
};

/// Configuration options for BFGS quasi-Newton molecular geometry optimization.
const BfgsOptions = struct {
    gradient: GradientOptions = .analytic,
    iterations: u32 = 100,
    step: f64 = 1,
    threshold: f64 = 1e-4,
};

/// Angular and radial grid point specifications for numerical DFT integration.
const DftGridOptions = struct {
    angular: usize = 302,
    radial: usize = 50,
};

/// Exchange-correlation functionals and integration grid settings for DFT.
const DftOptions = struct {
    correlation: ?[]const u8 = null,
    exchange: ?[]const u8 = null,
    exchange_correlation: ?[]const u8 = null,
    grid: DftGridOptions = .{},
};

/// Analytical nuclear gradient evaluation settings.
const GradientAnalyticOptions = struct {};

/// Finite difference displacement settings for numerical gradient evaluation.
const GradientNumericOptions = struct {
    step: f64 = 1e-5,
};

/// Finite difference displacement step for numerical Hessian evaluation.
const HessianNumericOptions = struct {
    step: f64 = 1e-5,
};

/// Convergence and iterative accelerator settings for coupled-perturbed HF response.
const ResponseOptions = struct {
    diis: ?u32 = 8,
    iterations: u32 = 100,
    threshold: f64 = 1e-8,
};

/// Configuration options for steepest descent molecular geometry optimization.
const SteepestDescentOptions = struct {
    gradient: GradientOptions = .analytic,
    iterations: u32 = 100,
    step: f64 = 1e-1,
    threshold: f64 = 1e-4,
};

/// File paths for exporting computed SCF matrices, molecular integrals, and geometries.
const Write = struct {
    coefficients: ?[]const u8 = null,
    density: ?[]const u8 = null,
    fock: ?[]const u8 = null,
    geometry: ?[]const u8 = null,
    gradient: ?[]const u8 = null,
    hessian: ?[]const u8 = null,
    integrals: MolecularIntegralsWrite = .{},
};

/// Output molecular orbitals, density, Fock matrix, orbital energies, and gradients from an SCF calculation.
pub fn Result(comptime T: type) type {
    return struct {
        ints: Integrals(T),

        C_mo: Matrix(T),
        e_mo: Vector(T),

        P: Matrix(T),
        F: Matrix(T),

        energy: []T,

        grad: []Matrix(T) = &.{},
        hess: []Matrix(T) = &.{},

        dC_mo: ?Tensor(T, 3) = null,

        de_mo: ?Matrix(T) = null,

        /// Frees allocated memory associated with the Hartree-Fock or DFT result structure.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.ints.deinit(gpa);

            self.C_mo.deinit(gpa);
            self.e_mo.deinit(gpa);

            self.P.deinit(gpa);
            self.F.deinit(gpa);

            gpa.free(self.energy);

            for (0..self.grad.len) |i| {
                self.grad[i].deinit(gpa);
            }

            if (self.grad.len > 0) gpa.free(self.grad);

            for (0..self.hess.len) |i| {
                self.hess[i].deinit(gpa);
            }

            if (self.hess.len > 0) gpa.free(self.hess);

            if (self.dC_mo) |*dC_mo| dC_mo.deinit(gpa);
            if (self.de_mo) |*de_mo| de_mo.deinit(gpa);
        }
    };
}

/// Container for reference pointers to density, Fock, and orbital matrices updated during SCF iterations.
fn ScfWorkspace(comptime T: type) type {
    return struct {
        P: *Matrix(T),
        F: *Matrix(T),
        C: *Matrix(T),
        e: *Vector(T),
    };
}

/// Evaluates the nuclear gradient of the total molecular energy from Hellmann-Feynman, two-electron, and Pulay forces.
pub fn calculateNuclearGradient(comptime T: type, ints: Integrals(T), ws: ScfWorkspace(T), generalized: bool, dft: ?*DftPotential(T), nthreads: usize, gpa: Allocator) !Matrix(T) {
    const dS = ints.dS orelse unreachable;
    const dH = ints.dH orelse unreachable;

    const nocc = if (generalized) ints.sys.nel else ints.sys.nel / 2;

    var G = try nuclearRepulsionGradient(T, ints.sys, gpa);
    errdefer G.deinit(gpa);

    var exch_factor: T = if (generalized) 1 else 0.5;

    if (dft) |pot| {
        exch_factor *= pot.exx_coef;
    }

    var W = try Matrix(T).init(ws.P.shape[0], ws.P.shape[0], gpa);
    defer W.deinit(gpa);

    buildEnergyWeightedDensity(T, &W, ws.C.*, ws.e.*, nocc, generalized);

    for (0..ints.sys.atoms.len) |i| for (0..3) |j| {
        const offset = (3 * i + j) * dS.shape[1] * dS.shape[2];

        const h_G = dot(T, ws.P.asVector(), Vector(T).fromSlice(dH.data[offset .. offset + ws.P.data.len]));
        const s_G = dot(T, W.asVector(), Vector(T).fromSlice(dS.data[offset .. offset + W.data.len]));

        var g_G: T = 0;

        if (ints.dg) |dg| {
            for (0..dg.shape[1]) |p| for (0..dg.shape[2]) |q| for (0..dg.shape[3]) |r| for (0..dg.shape[4]) |s| {
                const dg1 = dg.at(.{ 3 * i + j, p, r, q, s });
                const dg2 = dg.at(.{ 3 * i + j, p, q, r, s });

                g_G += 0.5 * ws.P.at(p, q) * ws.P.at(r, s) * (dg1 - exch_factor * dg2);
            };
        }

        G.ptr(i, j).* += h_G + g_G - s_G;
    };

    if (ints.dg == null) {
        if (generalized) {
            ints.sys.coulombGradientGhf(&G, ws.P.*, exch_factor, nthreads);
        }

        if (!generalized) {
            ints.sys.coulombGradientRhf(&G, ws.P.*, exch_factor, nthreads);
        }
    }

    return G;
}

/// Direct Inversion in the Iterative Subspace (DIIS) to accelerate Fock matrix convergence via error minimization.
pub fn diis(comptime T: type, fck_hist: []const Matrix(T), err_hist: []const Matrix(T), F: *Matrix(T), symmetric: bool, gpa: Allocator) !void {
    if (fck_hist.len < 2) return;

    var B = try Matrix(T).initZero(fck_hist.len + 1, fck_hist.len + 1, gpa);
    defer B.deinit(gpa);

    var b = try Matrix(T).initZero(fck_hist.len + 1, 1, gpa);
    defer b.deinit(gpa);

    for (0..fck_hist.len) |i| for (i..fck_hist.len) |j| {
        const sum = dot(T, err_hist[i].asVector(), err_hist[j].asVector());

        B.ptr(i, j).* = sum;
        B.ptr(j, i).* = sum;
    };

    for (0..fck_hist.len) |i| {
        B.ptr(i, fck_hist.len).* = -1.0;
        B.ptr(fck_hist.len, i).* = -1.0;
    }

    for (0..fck_hist.len) |i| {
        b.ptr(i, 0).* = 0.0;
    }

    b.ptr(fck_hist.len, 0).* = -1.0;

    const ipiv = try gpa.alloc(i32, fck_hist.len + 1);
    defer gpa.free(ipiv);

    try getrf(T, &B, ipiv);

    var c = try Matrix(T).init(fck_hist.len + 1, 1, gpa);
    defer c.deinit(gpa);

    try getrs(T, &c, B, ipiv, b);

    F.zero();

    for (0..fck_hist.len) |i| {
        var f_vec = F.asVector();

        axpy(T, &f_vec, fck_hist[i].asVector(), c.at(i, 0));
    }

    if (symmetric) F.symmetrize();
}

/// Computes the classical Coulomb repulsion gradient between nuclei.
pub fn nuclearRepulsionGradient(comptime T: type, sys: MolecularSystem(T), gpa: Allocator) !Matrix(T) {
    var dVN = try Matrix(T).initZero(sys.atoms.len, 3, gpa);
    errdefer dVN.deinit(gpa);

    for (0..sys.atoms.len) |i| {
        const Zi = @as(T, @floatFromInt(sys.atoms[i]));

        const xi = sys.coors[3 * i + 0];
        const yi = sys.coors[3 * i + 1];
        const zi = sys.coors[3 * i + 2];

        for (0..sys.atoms.len) |k| {
            if (k == i) continue;

            const Zk = @as(T, @floatFromInt(sys.atoms[k]));

            const xk = sys.coors[3 * k + 0];
            const yk = sys.coors[3 * k + 1];
            const zk = sys.coors[3 * k + 2];

            const dx = xi - xk;
            const dy = yi - yk;
            const dz = zi - zk;

            const dist = std.math.sqrt(dx * dx + dy * dy + dz * dz);

            dVN.ptr(i, 0).* -= Zi * Zk * dx / (dist * dist * dist);
            dVN.ptr(i, 1).* -= Zi * Zk * dy / (dist * dist * dist);
            dVN.ptr(i, 2).* -= Zi * Zk * dz / (dist * dist * dist);
        }
    }

    return dVN;
}

/// Computes the self-consistent field molecular orbitals, density matrix, and ground-state electronic energy.
pub fn solve(comptime T: type, io: std.Io, opt: Options, sys: *MolecularSystem(T), guess: ?Matrix(T), log: bool, gpa: Allocator) !Result(T) {
    const analytic_grad = if (opt.gradient) |g| g == .analytic else false;

    const molopts = MolecularIntegralsOptions{
        .system = opt.system,
        .basis = opt.basis,
        .generalized = opt.generalized,
        .charge = opt.charge,
        .multiplicity = opt.multiplicity,
        .nthreads = opt.nthreads,
        .calculate = .{
            .coulomb = !opt.integral_direct,
            .coulomb_d1 = !opt.integral_direct and analytic_grad,
            .hamiltonian_d1 = analytic_grad,
            .kinetic_d1 = analytic_grad,
            .nuclear_d1 = analytic_grad,
            .overlap_d1 = analytic_grad,
        },
        .write = opt.write.integrals,
    };

    var ints = try molecular_integrals_runFromSystem(T, io, molopts, sys.*, log, gpa);
    errdefer ints.deinit(gpa);

    var energy = try gpa.alloc(T, 1);
    errdefer gpa.free(energy);

    var grad: []Matrix(T) = &.{};

    errdefer {
        for (0..grad.len) |k| grad[k].deinit(gpa);

        if (grad.len > 0) gpa.free(grad);
    }

    const VN = try sys.nrep();

    const nocc = if (opt.generalized) sys.nel else blk: {
        if (sys.nel % 2 != 0) {
            return error.OnlyClosedShellSupported;
        }

        break :blk sys.nel / 2;
    };

    const nbf = if (opt.generalized) 2 * sys.nbf else sys.nbf;

    if (log) {
        try printf(io, "\nNUMBER OF BASIS FUNCTIONS: {d}, NUMBER OF OCCUPIED ORBITALS: {d}\n", .{ nbf, nocc });
    }

    var B = try Matrix(T).init(nbf, nbf, gpa);
    defer B.deinit(gpa);

    var C = try Matrix(T).init(nbf, nbf, gpa);
    errdefer C.deinit(gpa);

    var P = try Matrix(T).initZero(nbf, nbf, gpa);
    errdefer P.deinit(gpa);

    var F = try Matrix(T).init(nbf, nbf, gpa);
    errdefer F.deinit(gpa);

    var e = try Vector(T).init(nbf, gpa);
    errdefer e.deinit(gpa);

    var dft: ?DftPotential(T) = null;

    if (opt.dft) |dft_opt| {
        const n_rad, const n_leb = .{ dft_opt.grid.radial, dft_opt.grid.angular };

        const funcs = .{ dft_opt.exchange, dft_opt.correlation, dft_opt.exchange_correlation };

        dft = try DftPotential(T).init(sys.*, funcs, n_rad, n_leb, opt.generalized, gpa);
    }

    defer if (dft) |*pot| {
        pot.deinit(gpa);
    };

    if (guess) |p| {
        @memcpy(P.data, p.data);
    }

    if (guess == null) {
        try solveRoothaanHall(T, &e, &C, ints.H.?, ints.S.?, &B);

        buildDensityMatrix(T, &P, C, nocc, opt.generalized);
    }

    const ws: ScfWorkspace(T) = .{ .P = &P, .F = &F, .C = &C, .e = &e };

    energy[0] = try scf(T, io, opt, ints, ws, if (dft) |*d| d else null, log, gpa);

    if (log and opt.generalized) {
        const s2 = try calculateTotalSpin(T, P, ints.S.?, gpa);

        const method_str = if (dft) |_| "DFT" else "HARTREE-FOCK";

        try printTotalSpin(T, io, s2, opt.multiplicity, method_str);
    }

    if (log and opt.lowdin) {
        var charges = try lowdin(T, sys.*, P, ints.S.?, gpa);
        defer charges.deinit(gpa);

        const method_str = if (dft) |_| "DFT" else "HARTREE-FOCK";

        try printLowdinCharges(T, io, sys.*, charges, method_str);
    }

    if (log and opt.mayer) {
        var bo = try mayer(T, sys.*, P, ints.S.?, gpa);
        defer bo.deinit(gpa);

        const method_str = if (dft) |_| "DFT" else "HARTREE-FOCK";

        try printMayerBondOrders(T, io, sys.*, bo, method_str);
    }

    if (log and opt.mulliken) {
        var charges = try mulliken(T, sys.*, P, ints.S.?, gpa);
        defer charges.deinit(gpa);

        const method_str = if (dft) |_| "DFT" else "HARTREE-FOCK";

        try printMullikenCharges(T, io, sys.*, charges, method_str);
    }

    if (log and opt.wiberg) {
        var bo = try wiberg(T, sys.*, P, ints.S.?, gpa);
        defer bo.deinit(gpa);

        const method_str = if (dft) |_| "DFT" else "HARTREE-FOCK";

        try printWibergBondOrders(T, io, sys.*, bo, method_str);
    }

    if (log) {
        try printf(io, "\nNUCLEAR REPULSION ENERGY: {d:.14} Eh\n", .{VN});
    }

    if (log) if (dft) |*pot| {
        const names = try pot.getFunctionalNames(gpa);
        defer gpa.free(names);

        try printf(io, "\nFINAL DFT ENERGY ({s}): {d:.14} Eh\n", .{ names, energy[0] });
    };

    if (log and dft == null) {
        try printf(io, "\nFINAL HARTREE-FOCK ENERGY: {d:.14} Eh\n", .{energy[0]});
    }

    try exportMatrices(T, io, opt.write, C, P, F);

    if (analytic_grad) {
        const g0 = try calculateNuclearGradient(T, ints, ws, opt.generalized, if (dft) |*d| d else null, opt.nthreads, gpa);

        grad = try gpa.alloc(Matrix(T), 1);

        grad[0] = g0;
    }

    var result: Result(T) = .{
        .ints = ints,

        .P = P,
        .F = F,
        .e_mo = e,

        .energy = energy,

        .grad = grad,
        .hess = &.{},

        .C_mo = C,
    };

    if (opt.response) |response| {
        result.dC_mo, result.de_mo = try orbitalResponse(T, io, result, response, log, gpa);
    }

    return result;
}

/// Constructs the one-particle reduced density matrix in the atomic basis from molecular orbital coefficients.
fn buildDensityMatrix(comptime T: type, P: *Matrix(T), C: Matrix(T), nocc: usize, generalized: bool) void {
    std.debug.assert(C.shape[0] == P.shape[0]);
    std.debug.assert(C.shape[1] == P.shape[1]);

    const factor: T = if (generalized) 1 else 2;

    for (0..P.shape[0]) |i| for (0..P.shape[1]) |j| {
        var sum: T = 0;

        for (0..nocc) |k| {
            sum += factor * C.at(i, k) * C.at(j, k);
        }

        P.ptr(i, j).* = sum;
    };
}

/// Evaluates the Pulay energy-weighted density matrix from orbital energies and molecular orbital coefficients.
fn buildEnergyWeightedDensity(comptime T: type, W: *Matrix(T), C: Matrix(T), e: Vector(T), nocc: usize, generalized: bool) void {
    const factor: T = if (generalized) 1 else 2;

    for (0..W.nrow()) |i| for (0..W.ncol()) |j| {
        var sum: T = 0;

        for (0..nocc) |k| {
            sum += factor * e.at(k) * C.at(i, k) * C.at(j, k);
        }

        W.ptr(i, j).* = sum;
    };
}

/// Constructs the effective one-electron Fock matrix from core Hamiltonian, Coulomb, exchange, and DFT potentials.
fn buildFockMatrix(comptime T: type, F: *Matrix(T), ints: Integrals(T), P: Matrix(T), opt: Options, dft: ?*DftPotential(T), gpa: Allocator) !void {
    const nbf = ints.H.?.shape[0];

    std.debug.assert(F.shape[0] == nbf);
    std.debug.assert(F.shape[1] == nbf);
    std.debug.assert(P.shape[0] == nbf);
    std.debug.assert(P.shape[1] == nbf);

    for (0..nbf) |i| for (0..nbf) |j| {
        F.ptr(i, j).* = ints.H.?.at(i, j);
    };

    var exch_factor: T = if (opt.generalized) 1 else 0.5;

    if (dft) |pot| {
        exch_factor *= pot.exx_coef;
    }

    if (opt.integral_direct) {
        if (opt.generalized) {
            ints.sys.fockGhf(F, P, exch_factor, opt.nthreads);
        }

        if (!opt.generalized) {
            ints.sys.fockRhf(F, P, exch_factor, opt.nthreads);
        }
    }

    if (!opt.integral_direct) {
        contractTwoElectronFock(T, F, ints.g.?, P, exch_factor);
    }

    if (dft) |pot| {
        try pot.evaluate(ints.sys, P, gpa);

        for (0..nbf) |j| for (0..nbf) |k| {
            F.ptr(j, k).* += pot.Vxc.at(j, k);
        };
    }

    F.symmetrize();
}

/// Computes the electronic energy expectation value from density, core Hamiltonian, and Fock matrices.
fn calculateElectronicEnergy(comptime T: type, ints: Integrals(T), F: Matrix(T), P: Matrix(T), dft: ?*DftPotential(T)) T {
    var energy: T = if (dft) |pot| pot.Exc else 0;

    if (dft) |pot| for (0..P.shape[0]) |j| for (0..P.shape[1]) |k| {
        energy += 0.5 * P.at(j, k) * (ints.H.?.at(j, k) + F.at(j, k) - pot.Vxc.at(j, k));
    };

    if (dft == null) for (0..P.shape[0]) |i| for (0..P.shape[1]) |j| {
        energy += 0.5 * P.at(i, j) * (ints.H.?.at(i, j) + F.at(i, j));
    };

    return energy;
}

/// Computes the root-mean-square difference between current and previous density matrices.
fn computeDensityRmsDifference(comptime T: type, P_new: Matrix(T), P_old: Matrix(T)) T {
    std.debug.assert(P_new.shape[0] == P_old.shape[0]);
    std.debug.assert(P_new.shape[1] == P_old.shape[1]);

    var sum_sq: T = 0;

    for (0..P_new.data.len) |i| {
        const diff = P_new.data[i] - P_old.data[i];

        sum_sq += diff * diff;
    }

    return @sqrt(sum_sq / @as(T, @floatFromInt(P_new.data.len)));
}

/// Evaluates the DIIS orbital gradient error matrix representing the commutator [F, P] in the non-orthogonal AO basis.
fn computeOrbitalGradient(comptime T: type, err: *Matrix(T), F: Matrix(T), P: Matrix(T), S: Matrix(T), gpa: Allocator) !void {
    const nbf = F.shape[0];

    std.debug.assert(F.shape[1] == nbf);
    std.debug.assert(P.shape[0] == nbf);
    std.debug.assert(P.shape[1] == nbf);
    std.debug.assert(S.shape[0] == nbf);
    std.debug.assert(S.shape[1] == nbf);

    std.debug.assert(err.shape[0] == nbf);
    std.debug.assert(err.shape[1] == nbf);

    var FP = try Matrix(T).init(nbf, nbf, gpa);
    defer FP.deinit(gpa);

    gemm(T, &FP, F, P, 1, 0, false, false);

    var FPS = try Matrix(T).init(nbf, nbf, gpa);
    defer FPS.deinit(gpa);

    gemm(T, &FPS, FP, S, 1, 0, false, false);

    for (0..nbf) |i| for (0..nbf) |j| {
        err.ptr(i, j).* = FPS.at(i, j) - FPS.at(j, i);
    };
}

/// Contracts four-center two-electron repulsion integrals with the density matrix to form Coulomb and exchange terms.
fn contractTwoElectronFock(comptime T: type, F: *Matrix(T), g: Tensor(T, 4), P: Matrix(T), exch_factor: T) void {
    const nbf = F.nrow();

    std.debug.assert(g.shape[0] == nbf);
    std.debug.assert(g.shape[1] == nbf);
    std.debug.assert(g.shape[2] == nbf);
    std.debug.assert(g.shape[3] == nbf);
    std.debug.assert(P.shape[0] == nbf);
    std.debug.assert(P.shape[1] == nbf);

    for (0..nbf) |i| for (0..nbf) |j| for (0..nbf) |k| for (0..nbf) |l| {
        F.ptr(k, l).* += P.at(i, j) * (g.at(.{ i, k, j, l }) - exch_factor * g.at(.{ i, j, k, l }));
    };
}

/// Exports SCF result matrices to files if corresponding paths are provided.
fn exportMatrices(comptime T: type, io: std.Io, write: Write, C: Matrix(T), P: Matrix(T), F: Matrix(T)) !void {
    if (write.coefficients) |fname| {
        try writeMatrix(T, io, fname, C);
    }

    if (write.density) |fname| {
        try writeMatrix(T, io, fname, P);
    }

    if (write.fock) |fname| {
        try writeMatrix(T, io, fname, F);
    }
}

/// Solves the self-consistent field equations iteratively using a density-driven approach.
fn scf(comptime T: type, io: std.Io, opt: Options, ints: Integrals(T), ws: ScfWorkspace(T), dft: ?*DftPotential(T), log: bool, gpa: Allocator) !T {
    const VN = try ints.sys.nrep();

    const nocc = if (opt.generalized) ints.sys.nel else blk: {
        if (ints.sys.nel % 2 != 0) {
            return error.OnlyClosedShellSupported;
        }

        break :blk ints.sys.nel / 2;
    };

    const nbf = if (opt.generalized) 2 * ints.sys.nbf else ints.sys.nbf;

    var e_old: T = VN;
    var e_new: T = VN;

    var B = try Matrix(T).init(nbf, nbf, gpa);
    defer B.deinit(gpa);

    var P_old = try Matrix(T).init(nbf, nbf, gpa);
    defer P_old.deinit(gpa);

    if (opt.iterations > 0 and log) {
        const fmt = "\nSELF CONSISTENT FIELD\n{s:4} {s:20} {s:9} {s:9} {s:9}\n";

        try printf(io, fmt, .{ "ITER", "TOTAL ENERGY (Eh)", "|DE| (Eh)", "RMS(DP)", "TIME" });
    }

    var fck_hist = std.ArrayList(Matrix(T)).empty;

    defer {
        for (0..fck_hist.items.len) |i| fck_hist.items[i].deinit(gpa);

        fck_hist.deinit(gpa);
    }

    var err_hist = std.ArrayList(Matrix(T)).empty;

    defer {
        for (0..err_hist.items.len) |i| err_hist.items[i].deinit(gpa);

        err_hist.deinit(gpa);
    }

    for (0..opt.iterations) |i| {
        var timer = std.Io.Timestamp.now(io, .real);

        try buildFockMatrix(T, ws.F, ints, ws.P.*, opt, dft, gpa);

        e_new = calculateElectronicEnergy(T, ints, ws.F.*, ws.P.*, dft) + VN;

        if (opt.diis) |d| if (d > 0 and i > 0) {
            try updateDiis(T, ws, ints.S.?, &fck_hist, &err_hist, d, gpa);
        };

        try solveRoothaanHall(T, ws.e, ws.C, ws.F.*, ints.S.?, &B);

        @memcpy(P_old.data, ws.P.data);

        buildDensityMatrix(T, ws.P, ws.C.*, nocc, opt.generalized);

        const delta_energy, const p_rm = .{ @abs(e_new - e_old), computeDensityRmsDifference(T, ws.P.*, P_old) };

        const elapsed = timer.untilNow(io, .real);

        if (log) {
            const fmt = "{d:4} {d:20.14} {e:9.3} {e:9.3} {s:9} {s}\n";

            const time_str = try std.fmt.allocPrint(gpa, "{f}", .{elapsed});
            defer gpa.free(time_str);

            const de = if (i == 0) 0 else delta_energy;

            try printf(io, fmt, .{ i + 1, e_new, de, p_rm, time_str, if (err_hist.items.len >= 2) "DIIS" else "" });
        }

        e_old = e_new;

        if (i > 0 and delta_energy < opt.threshold and p_rm < opt.threshold) {
            break;
        }
    } else return error.ScfDidNotConverge;

    return e_new;
}

/// Solves the generalized Roothaan-Hall eigenvalue problem FC = SCE for molecular orbital coefficients and energies.
fn solveRoothaanHall(comptime T: type, e: *Vector(T), C: *Matrix(T), F: Matrix(T), S: Matrix(T), B: *Matrix(T)) !void {
    @memcpy(B.data, S.data);

    try sygvd(T, e, C, F, B);
}

/// Accelerates SCF convergence by updating Fock and error matrix histories and extrapolating via DIIS.
fn updateDiis(comptime T: type, ws: ScfWorkspace(T), S: Matrix(T), fck_hist: *std.ArrayList(Matrix(T)), err_hist: *std.ArrayList(Matrix(T)), max_diis: usize, gpa: Allocator) !void {
    const nbf = ws.F.nrow();

    var e_diis = try Matrix(T).init(nbf, nbf, gpa);
    errdefer e_diis.deinit(gpa);

    try computeOrbitalGradient(T, &e_diis, ws.F.*, ws.P.*, S, gpa);

    try fck_hist.ensureUnusedCapacity(gpa, 1);
    try err_hist.ensureUnusedCapacity(gpa, 1);

    var f_diis = try Matrix(T).init(nbf, nbf, gpa);
    errdefer f_diis.deinit(gpa);

    @memcpy(f_diis.data, ws.F.data);

    if (fck_hist.items.len >= max_diis) {
        var old_f = fck_hist.orderedRemove(0);
        var old_e = err_hist.orderedRemove(0);

        old_f.deinit(gpa);
        old_e.deinit(gpa);
    }

    fck_hist.appendAssumeCapacity(f_diis);
    err_hist.appendAssumeCapacity(e_diis);

    diis(T, fck_hist.items, err_hist.items, ws.F, true, gpa) catch {
        for (0..fck_hist.items.len) |j| fck_hist.items[j].deinit(gpa);
        for (0..err_hist.items.len) |j| err_hist.items[j].deinit(gpa);

        fck_hist.clearRetainingCapacity();
        err_hist.clearRetainingCapacity();

        try fck_hist.ensureUnusedCapacity(gpa, 1);
        try err_hist.ensureUnusedCapacity(gpa, 1);

        var f_retry = try Matrix(T).init(nbf, nbf, gpa);
        errdefer f_retry.deinit(gpa);

        var e_retry = try Matrix(T).init(nbf, nbf, gpa);
        errdefer e_retry.deinit(gpa);

        @memcpy(f_retry.data, ws.F.data);

        try computeOrbitalGradient(T, &e_retry, ws.F.*, ws.P.*, S, gpa);

        fck_hist.appendAssumeCapacity(f_retry);
        err_hist.appendAssumeCapacity(e_retry);
    };
}
