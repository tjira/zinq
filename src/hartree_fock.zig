//! Solves the Hartree-Fock self-consistent field equations to approximate the electronic wavefunction of a system.

const std = @import("std");

const Allocator = std.mem.Allocator;

const DftPotential = @import("density_functional_theory.zig").DftPotential;
const Diis = @import("diis.zig").Diis;
const FrequencyOptions = @import("frequency_analysis.zig").Options;
const Integrals = @import("molecular_integrals.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const MolecularIntegralsOptions = @import("molecular_integrals.zig").Options;
const MolecularIntegralsWrite = @import("molecular_integrals.zig").Write;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

const calculateTotalSpin = @import("spin_analysis.zig").calculateTotalSpin;
const contractTwoElectronFock = @import("two_electron_contraction.zig").contractTwoElectronFock;
const contractTwoElectronGradient = @import("two_electron_contraction.zig").contractTwoElectronGradient;
const dot = @import("linear_algebra.zig").dot;
const gemm = @import("linear_algebra.zig").gemm;
const getSymbol = @import("constants.zig").getSymbol;
const lowdin = @import("population_analysis.zig").lowdin;
const mayer = @import("population_analysis.zig").mayer;
const molecular_integrals_run = @import("molecular_integrals.zig").run;
const molecular_integrals_runFromSystem = @import("molecular_integrals.zig").runFromSystem;
const mulliken = @import("population_analysis.zig").mulliken;
const orbitalResponse = @import("orbital_response.zig").orbitalResponse;
const printLowdinCharges = @import("population_analysis.zig").printLowdinCharges;
const printMayerBondOrders = @import("population_analysis.zig").printMayerBondOrders;
const printMullikenCharges = @import("population_analysis.zig").printMullikenCharges;
const printTotalSpin = @import("spin_analysis.zig").printTotalSpin;
const printWibergBondOrders = @import("population_analysis.zig").printWibergBondOrders;
const printf = @import("read_write.zig").printf;
const sygvd = @import("linear_algebra.zig").sygvd;
const wiberg = @import("population_analysis.zig").wiberg;
const writeMatrix = @import("read_write.zig").writeMatrix;

const AN2SM = @import("constants.zig").AN2SM;
const AU2CM = @import("constants.zig").AU2CM;

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

        P_ao: Matrix(T),
        F_ao: Matrix(T),

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

            self.P_ao.deinit(gpa);
            self.F_ao.deinit(gpa);

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

    var G = try ints.sys.nuclearRepulsionGradient(gpa);
    errdefer G.deinit(gpa);

    var exch_factor: T = if (generalized) 1 else 0.5;

    if (dft) |pot| {
        exch_factor *= pot.a_x;
    }

    var W = try Matrix(T).init(ws.P.shape[0], ws.P.shape[0], gpa);
    defer W.deinit(gpa);

    buildEnergyWeightedDensity(T, &W, ws.C.*, ws.e.*, nocc, generalized);

    for (0..ints.sys.atoms.len) |A| for (0..3) |alpha| {
        const offset = (3 * A + alpha) * dS.shape[1] * dS.shape[2];

        const h_G = dot(T, ws.P.asVector(), Vector(T).fromSlice(dH.data[offset .. offset + ws.P.data.len]));

        const s_G = dot(T, W.asVector(), Vector(T).fromSlice(dS.data[offset .. offset + W.data.len]));

        var g_G: T = 0;

        if (ints.dg) |dg| {
            g_G += contractTwoElectronGradient(T, dg, 3 * A + alpha, ws.P.*, exch_factor);
        }

        G.ptr(A, alpha).* += h_G + g_G - s_G;
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

    if (log) {
        try computeElectronicProperties(T, io, opt, sys.*, P, ints.S.?, if (dft) |*d| d else null, gpa);

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
        const pot = if (dft) |*d| d else null;

        const g0 = try calculateNuclearGradient(T, ints, ws, opt.generalized, pot, opt.nthreads, gpa);

        grad = try gpa.alloc(Matrix(T), 1);

        grad[0] = g0;
    }

    var result: Result(T) = .{
        .ints = ints,

        .P_ao = P,
        .F_ao = F,
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

/// Evaluates the exchange-correlation potential and adds its contribution to the Fock matrix.
fn addDftPotential(comptime T: type, F: *Matrix(T), dft: *DftPotential(T), P: Matrix(T), sys: MolecularSystem(T), gpa: Allocator) !void {
    try dft.evaluate(sys, P, gpa);

    const nbf = F.nrow();

    for (0..nbf) |mu| for (0..nbf) |nu| {
        F.ptr(mu, nu).* += dft.V_xc.at(mu, nu);
    };
}

/// Constructs the one-particle reduced density matrix in the atomic basis from molecular orbital coefficients.
fn buildDensityMatrix(comptime T: type, P: *Matrix(T), C: Matrix(T), nocc: usize, generalized: bool) void {
    std.debug.assert(C.shape[0] == P.shape[0]);
    std.debug.assert(C.shape[1] == P.shape[1]);

    const factor: T = if (generalized) 1 else 2;

    for (0..P.shape[0]) |mu| for (0..P.shape[1]) |nu| {
        var sum: T = 0;

        for (0..nocc) |i| {
            sum += factor * C.at(mu, i) * C.at(nu, i);
        }

        P.ptr(mu, nu).* = sum;
    };
}

/// Evaluates the Pulay energy-weighted density matrix from orbital energies and molecular orbital coefficients.
fn buildEnergyWeightedDensity(comptime T: type, W: *Matrix(T), C: Matrix(T), e: Vector(T), nocc: usize, generalized: bool) void {
    const factor: T = if (generalized) 1 else 2;

    for (0..W.nrow()) |mu| for (0..W.ncol()) |nu| {
        var sum: T = 0;

        for (0..nocc) |i| {
            sum += factor * e.at(i) * C.at(mu, i) * C.at(nu, i);
        }

        W.ptr(mu, nu).* = sum;
    };
}

/// Constructs the effective one-electron Fock matrix from core Hamiltonian, Coulomb, exchange, and DFT potentials.
fn buildFockMatrix(comptime T: type, F: *Matrix(T), ints: Integrals(T), P: Matrix(T), opt: Options, dft: ?*DftPotential(T), gpa: Allocator) !void {
    const nbf = ints.H.?.shape[0];

    std.debug.assert(F.shape[0] == nbf);
    std.debug.assert(F.shape[1] == nbf);
    std.debug.assert(P.shape[0] == nbf);
    std.debug.assert(P.shape[1] == nbf);

    for (0..nbf) |mu| for (0..nbf) |nu| {
        F.ptr(mu, nu).* = ints.H.?.at(mu, nu);
    };

    var exch_factor: T = if (opt.generalized) 1 else 0.5;

    if (dft) |pot| {
        exch_factor *= pot.a_x;
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
        try addDftPotential(T, F, pot, P, ints.sys, gpa);
    }

    F.symmetrize();
}

/// Computes the electronic energy expectation value from density, core Hamiltonian, and Fock matrices.
fn calculateElectronicEnergy(comptime T: type, ints: Integrals(T), F: Matrix(T), P: Matrix(T), dft: ?*DftPotential(T)) T {
    var energy: T = if (dft) |pot| pot.E_xc else 0;

    if (dft) |pot| for (0..P.shape[0]) |mu| for (0..P.shape[1]) |nu| {
        energy += 0.5 * P.at(mu, nu) * (ints.H.?.at(mu, nu) + F.at(mu, nu) - pot.V_xc.at(mu, nu));
    };

    if (dft == null) for (0..P.shape[0]) |mu| for (0..P.shape[1]) |nu| {
        energy += 0.5 * P.at(mu, nu) * (ints.H.?.at(mu, nu) + F.at(mu, nu));
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

/// Computes and prints requested population and spin properties from the converged density matrix.
fn computeElectronicProperties(comptime T: type, io: std.Io, opt: Options, sys: MolecularSystem(T), P: Matrix(T), S: Matrix(T), dft: ?*DftPotential(T), gpa: Allocator) !void {
    const method_str = if (dft) |_| "DFT" else "HARTREE-FOCK";

    if (opt.generalized) {
        const s2 = try calculateTotalSpin(T, P, S, gpa);

        try printTotalSpin(T, io, s2, opt.multiplicity, method_str);
    }

    if (opt.lowdin) {
        var charges = try lowdin(T, sys, P, S, gpa);
        defer charges.deinit(gpa);

        try printLowdinCharges(T, io, sys, charges, method_str);
    }

    if (opt.mayer) {
        var bo = try mayer(T, sys, P, S, gpa);
        defer bo.deinit(gpa);

        try printMayerBondOrders(T, io, sys, bo, method_str);
    }

    if (opt.mulliken) {
        var charges = try mulliken(T, sys, P, S, gpa);
        defer charges.deinit(gpa);

        try printMullikenCharges(T, io, sys, charges, method_str);
    }

    if (opt.wiberg) {
        var bo = try wiberg(T, sys, P, S, gpa);
        defer bo.deinit(gpa);

        try printWibergBondOrders(T, io, sys, bo, method_str);
    }
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

    for (0..nbf) |mu| for (0..nbf) |nu| {
        err.ptr(mu, nu).* = FPS.at(mu, nu) - FPS.at(nu, mu);
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

    var diis_acc: ?Diis(T) = if (opt.diis) |d| if (d > 0) Diis(T).init(d) else null else null;
    defer if (diis_acc) |*d| d.deinit(gpa);

    for (0..opt.iterations) |i| {
        var timer = std.Io.Timestamp.now(io, .real);

        try buildFockMatrix(T, ws.F, ints, ws.P.*, opt, dft, gpa);

        e_new = calculateElectronicEnergy(T, ints, ws.F.*, ws.P.*, dft) + VN;

        if (diis_acc) |*d| if (i > 0) {
            var e_diis = try Matrix(T).init(nbf, nbf, gpa);
            errdefer e_diis.deinit(gpa);

            try computeOrbitalGradient(T, &e_diis, ws.F.*, ws.P.*, ints.S.?, gpa);

            try d.update(ws.F, e_diis, true, gpa);
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

            const diis_tag = if (diis_acc) |d| (if (d.err_hist.items.len >= 2) "DIIS" else "") else "";

            try printf(io, fmt, .{ i + 1, e_new, de, p_rm, time_str, diis_tag });
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
