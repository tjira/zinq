//! Solves the Hartree-Fock self-consistent field equations to approximate the electronic wavefunction of a system.

const std = @import("std");

const cblas = @import("cimport.zig").cblas;

const Allocator = std.mem.Allocator;

const DftPotential = @import("density_functional_theory.zig").DftPotential;
const FrequencyOptions = @import("frequency_analysis.zig").Options;
const Integrals = @import("molecular_integrals.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const MolecularIntegralsOptions = @import("molecular_integrals.zig").Options;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

const addScaled = @import("linear_algebra.zig").addScaled;
const ao2mo_pp = @import("integral_transform.zig").ao2mo_pp;
const bfgs = @import("molecular_optimization.zig").bfgs;
const calculateHarmonicFrequencies = @import("frequency_analysis.zig").calculateHarmonicFrequencies;
const calculateNumericalGradient = @import("nuclear_derivative.zig").calculateNumericalGradient;
const calculateNumericalHessian = @import("nuclear_derivative.zig").calculateNumericalHessian;
const calculateThermochemistry = @import("frequency_analysis.zig").calculateThermochemistry;
const calculateTotalSpin = @import("spin_analysis.zig").calculateTotalSpin;
const dot = @import("linear_algebra.zig").dot;
const exportIfBuiltin = @import("molecular_integrals.zig").exportIfBuiltin;
const geigh = @import("linear_algebra.zig").geigh;
const getSymbol = @import("constant.zig").getSymbol;
const lowdin = @import("population_analysis.zig").lowdin;
const luFactorize = @import("linear_algebra.zig").luFactorize;
const luSolve = @import("linear_algebra.zig").luSolve;
const mayer = @import("population_analysis.zig").mayer;
const mm = @import("linear_algebra.zig").mm;
const mo2ao_xx = @import("integral_transform.zig").mo2ao_xx;
const molecular_integrals_run = @import("molecular_integrals.zig").run;
const molecular_integrals_runFromSystem = @import("molecular_integrals.zig").runFromSystem;
const mulliken = @import("population_analysis.zig").mulliken;
const orbitalResponse = @import("cphf.zig").orbitalResponse;
const printHarmonicFrequencies = @import("frequency_analysis.zig").printHarmonicFrequencies;
const printLowdinCharges = @import("population_analysis.zig").printLowdinCharges;
const printMayerBondOrders = @import("population_analysis.zig").printMayerBondOrders;
const printMullikenCharges = @import("population_analysis.zig").printMullikenCharges;
const printThermochemistry = @import("frequency_analysis.zig").printThermochemistry;
const printTotalSpin = @import("spin_analysis.zig").printTotalSpin;
const printWibergBondOrders = @import("population_analysis.zig").printWibergBondOrders;
const printf = @import("read_write.zig").printf;
const steepestDescent = @import("molecular_optimization.zig").steepestDescent;
const wiberg = @import("population_analysis.zig").wiberg;
const writeMatrix = @import("read_write.zig").writeMatrix;
const writeXyzFile = @import("read_write.zig").writeXyzFile;

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

/// File paths for exporting computed SCF matrices and geometries.
const Write = struct {
    coefficients: ?[]const u8 = null,
    density: ?[]const u8 = null,
    fock: ?[]const u8 = null,
    geometry: ?[]const u8 = null,
    gradient: ?[]const u8 = null,
    hessian: ?[]const u8 = null,
};

/// Output molecular orbitals, density, Fock matrix, orbital energies, and gradients from an SCF calculation.
pub fn Result(comptime T: type) type {
    return struct {
        ints: Integrals(T),

        C: Matrix(T),
        P: Matrix(T),
        F: Matrix(T),
        e: Vector(T),

        energy: []T,

        grad: []Matrix(T) = &.{},
        hess: []Matrix(T) = &.{},

        dC: ?Tensor(T, 3) = null,

        de: ?Matrix(T) = null,

        /// Frees allocated memory associated with the Hartree-Fock or DFT result structure.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.ints.deinit(gpa);

            self.C.deinit(gpa);
            self.P.deinit(gpa);
            self.F.deinit(gpa);
            self.e.deinit(gpa);

            gpa.free(self.energy);

            for (0..self.grad.len) |i| {
                self.grad[i].deinit(gpa);
            }

            gpa.free(self.grad);

            for (0..self.hess.len) |i| {
                self.hess[i].deinit(gpa);
            }

            gpa.free(self.hess);

            if (self.dC) |*dC| dC.deinit(gpa);
            if (self.de) |*de| de.deinit(gpa);
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

    try luFactorize(T, &B, ipiv);

    var c = try Matrix(T).init(fck_hist.len + 1, 1, gpa);
    defer c.deinit(gpa);

    try luSolve(T, &c, B, ipiv, b);

    F.zero();

    for (0..fck_hist.len) |i| {
        var f_vec = F.asVector();

        addScaled(T, c.at(i, 0), fck_hist[i].asVector(), &f_vec);
    }

    if (symmetric) for (0..F.shape[0]) |i| for (i + 1..F.shape[1]) |j| {
        const avg = (F.at(i, j) + F.at(j, i)) / 2;

        F.ptr(i, j).* = avg;
        F.ptr(j, i).* = avg;
    };
}

/// Computes the nuclear gradient of the total energy with respect to atomic coordinates.
pub fn gradient(comptime T: type, ints: Integrals(T), ws: ScfWorkspace(T), generalized: bool, dft: ?*DftPotential(T), nthreads: usize, gpa: Allocator) !Matrix(T) {
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

    const factor: T = if (generalized) 1 else 2;

    for (0..W.nrow()) |i| for (0..W.ncol()) |j| {
        var sum: T = 0;

        for (0..nocc) |k| {
            sum += factor * ws.e.at(k) * ws.C.at(i, k) * ws.C.at(j, k);
        }

        W.ptr(i, j).* = sum;
    };

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

/// Executes a Hartree-Fock or DFT calculation on a molecular system specified by file paths.
pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator) !Result(T) {
    try checkInvalidInput(opt);

    cblas.openblas_set_num_threads(@intCast(opt.nthreads));

    const basis_path = try exportIfBuiltin(io, opt.basis, gpa);

    defer if (std.mem.startsWith(u8, opt.basis, "builtin:")) {
        std.Io.Dir.cwd().deleteFile(io, basis_path) catch {};
    };

    var sys = try MolecularSystem(T).init(io, opt.system, basis_path, opt.charge, opt.multiplicity, gpa);
    defer sys.deinit(gpa);

    if (std.mem.startsWith(u8, opt.basis, "builtin:")) {
        try std.Io.Dir.cwd().deleteFile(io, basis_path);
    }

    return try runFromSystem(T, io, opt, &sys, null, log, gpa);
}

/// Runs Hartree-Fock or DFT SCF starting from an initialized MolecularSystem structure.
pub fn runFromSystem(comptime T: type, io: std.Io, opt: Options, sys: *MolecularSystem(T), Pg: ?Matrix(T), log: bool, gpa: Allocator) !Result(T) {
    try checkInvalidInput(opt);

    cblas.openblas_set_num_threads(@intCast(opt.nthreads));

    var final_Pg: ?Matrix(T) = null;
    defer if (final_Pg) |*p| p.deinit(gpa);

    if (opt.optimize) |o| {
        switch (o) {
            .steepest_descent => final_Pg = try steepestDescent(T, io, runFromSystem, opt, sys, Pg, log, gpa),
            .bfgs => final_Pg = try bfgs(T, io, runFromSystem, opt, sys, Pg, log, gpa),
        }
    }

    const molopts = MolecularIntegralsOptions{
        .system = opt.system,
        .basis = opt.basis,
        .spin = opt.generalized,
        .charge = opt.charge,
        .multiplicity = opt.multiplicity,
        .nthreads = opt.nthreads,
        .calculate = .{
            .coulomb = !opt.integral_direct,
            .kinetic_d1 = opt.gradient != null and opt.gradient.? == .analytic,
            .overlap_d1 = opt.gradient != null and opt.gradient.? == .analytic,
            .nuclear_d1 = opt.gradient != null and opt.gradient.? == .analytic,
            .hmatrix_d1 = opt.gradient != null and opt.gradient.? == .analytic,
            .coulomb_d1 = !opt.integral_direct and opt.gradient != null and opt.gradient.? == .analytic,
        },
    };

    var ints = try molecular_integrals_runFromSystem(T, io, molopts, sys.*, log, gpa);
    errdefer ints.deinit(gpa);

    var energy = try gpa.alloc(T, 1);
    errdefer gpa.free(energy);

    var grad = try gpa.alloc(Matrix(T), if (opt.gradient) |_| 1 else 0);
    errdefer gpa.free(grad);

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

    const guess_p = final_Pg orelse Pg;

    if (guess_p) |guess| {
        @memcpy(P.data, guess.data);
    }

    if (guess_p == null) {
        @memcpy(B.data, ints.S.?.data);

        try geigh(T, &e, &C, ints.H.?, &B);

        _ = getDensity(T, &P, C, nocc, opt.generalized);
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

    if (opt.write.geometry) |fname| {
        try writeXyzFile(T, io, fname, sys.atoms, sys.coors);
    }

    if (opt.gradient) |gradopt| switch (gradopt) {
        .analytic => grad[0] = try gradient(T, ints, ws, opt.generalized, if (dft) |*d| d else null, opt.nthreads, gpa),
        .numeric => grad[0] = try calculateNumericalGradient(T, io, runFromSystem, opt, sys, log, gpa),
    };

    errdefer if (opt.gradient) |_| grad[0].deinit(gpa);

    if (opt.write.gradient) |fname| if (grad.len > 0) {
        try writeMatrix(T, io, fname, grad[0]);
    };

    if (log) for (0..grad.len) |i| {
        const grad_type_str = if (opt.gradient.? == .analytic) "ANALYTICAL" else "NUMERICAL";

        const method_str = if (dft) |_| "DFT" else "HARTREE-FOCK";

        try printf(io, "\n{s} {s} NUCLEAR ENERGY GRADIENT (Eh/a0)\n", .{ method_str, grad_type_str });

        for (0..grad[i].shape[0]) |j| for (0..grad[i].shape[1]) |k| {
            try printf(io, "{d:20.14}{s}", .{ grad[i].at(j, k), if (k == 2) "\n" else " " });
        };
    };

    const hess = try handleHessianAndFrequencies(T, io, opt, runFromSystem, sys, energy[0], log, gpa);

    errdefer {
        if (opt.hessian) |_| hess[0].deinit(gpa);

        gpa.free(hess);
    }

    if (opt.write.hessian) |fname| if (hess.len > 0) {
        try writeMatrix(T, io, fname, hess[0]);
    };

    var result: Result(T) = .{
        .ints = ints,
        .P = P,
        .C = C,
        .F = F,
        .e = e,
        .energy = energy,
        .grad = grad,
        .hess = hess,
    };

    if (opt.response) |response| {
        result.dC, result.de = try orbitalResponse(T, io, result, response, log, gpa);
    }

    return result;
}

/// Validates input options for physical consistency and method compatibility.
fn checkInvalidInput(opt: Options) !void {
    if (opt.integral_direct) {
        if (opt.response != null) {
            std.log.err("CPHF RESPONSE PROPERTIES ARE NOT SUPPORTED FOR INTEGRAL DIRECT HARTREE-FOCK", .{});

            return error.InvalidInput;
        }
    }

    if (opt.nthreads == 0) {
        std.log.err("THREAD COUNT MUST BE GREATER THAN 0", .{});

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

    if (opt.frequency != null and opt.hessian == null) {
        std.log.err("FREQUENCY CALCULATION REQUESTED BUT HESSIAN IS NOT CALCULATED", .{});

        return error.InvalidInput;
    }

    if (opt.multiplicity == 0) {
        std.log.err("MULTIPLICITY MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.multiplicity > 1 and !opt.generalized) {
        std.log.err("OPEN SHELL IS NOT SUPPORTED UNLESS GENERALIZED HARTREE-FOCK IS ENABLED", .{});

        return error.OnlyClosedShellSupported;
    }

    if (opt.system.len == 0) {
        std.log.err("MOLECULAR SYSTEM XYZ PATH IS EMPTY", .{});

        return error.InvalidInput;
    }

    if (opt.basis.len == 0) {
        std.log.err("BASIS SET G94 PATH IS EMPTY", .{});

        return error.InvalidInput;
    }

    if (opt.iterations == 0) {
        std.log.err("SCF ITERATIONS MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.threshold <= 0) {
        std.log.err("SCF CONVERGENCE THRESHOLD MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.dft) |d| {
        if (opt.gradient != null and opt.gradient.? == .analytic) {
            std.log.err("ANALYTIC GRADIENTS ARE NOT SUPPORTED FOR DFT CALCULATIONS", .{});

            return error.InvalidInput;
        }

        if (opt.response != null) {
            std.log.err("CPKS RESPONSE PROPERTIES ARE NOT SUPPORTED FOR DFT CALCULATIONS", .{});

            return error.InvalidInput;
        }

        if (d.grid.radial == 0) {
            std.log.err("DFT RADIAL GRID POINTS MUST BE GREATER THAN 0", .{});

            return error.InvalidInput;
        }

        const valid_angular_points = switch (d.grid.angular) {
            50, 74, 110, 302, 590, 974 => true,

            else => false,
        };

        if (!valid_angular_points) {
            std.log.err("DFT ANGULAR GRID SIZE IS NOT SUPPORTED. USE 50, 74, 110, 302, 590, OR 974", .{});

            return error.InvalidInput;
        }

        if (d.exchange_correlation != null and (d.exchange != null or d.correlation != null)) {
            std.log.err("CANNOT SPECIFY BOTH EXCHANGE_CORRELATION AND INDIVIDUAL FUNCTIONALS", .{});

            return error.InvalidInput;
        }

        if (d.exchange_correlation == null and d.exchange == null and d.correlation == null) {
            std.log.err("MUST SPECIFY EITHER EXCHANGE_CORRELATION OR EXCHANGE AND/OR CORRELATION FUNCTIONALS", .{});

            return error.InvalidInput;
        }
    }

    if (opt.response) |r| {
        if (r.iterations == 0) {
            std.log.err("CPHF RESPONSE ITERATIONS MUST BE GREATER THAN 0", .{});

            return error.InvalidInput;
        }

        if (r.threshold <= 0) {
            std.log.err("CPHF RESPONSE THRESHOLD MUST BE GREATER THAN 0", .{});

            return error.InvalidInput;
        }
    }
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

/// Computes the electronic density matrix from molecular orbital coefficients.
fn getDensity(comptime T: type, P: *Matrix(T), C: Matrix(T), nocc: usize, generalized: bool) T {
    std.debug.assert(C.shape[0] == P.shape[0]);
    std.debug.assert(C.shape[1] == P.shape[1]);

    const factor: T, var sum_sq_dp: T = .{ if (generalized) 1 else 2, 0 };

    for (0..P.shape[0]) |i| for (0..P.shape[1]) |j| {
        var sum: T = 0;

        for (0..nocc) |k| {
            sum += factor * C.at(i, k) * C.at(j, k);
        }

        sum_sq_dp += (sum - P.at(i, j)) * (sum - P.at(i, j));

        P.ptr(i, j).* = sum;
    };

    return @sqrt(sum_sq_dp / @as(T, @floatFromInt(P.shape[0] * P.shape[1])));
}

/// Computes the total electronic energy from the Fock and density matrices.
fn getEnergy(comptime T: type, ints: Integrals(T), F: Matrix(T), P: Matrix(T), dft: ?*DftPotential(T)) T {
    var energy: T = if (dft) |pot| pot.Exc else 0;

    if (dft) |pot| for (0..P.shape[0]) |j| for (0..P.shape[1]) |k| {
        energy += 0.5 * P.at(j, k) * (ints.H.?.at(j, k) + F.at(j, k) - pot.Vxc.at(j, k));
    };

    if (dft == null) for (0..P.shape[0]) |i| for (0..P.shape[1]) |j| {
        energy += 0.5 * P.at(i, j) * (ints.H.?.at(i, j) + F.at(i, j));
    };

    return energy;
}

/// Computes the DIIS error matrix representing the commutator [F, P] in the overlap representation.
fn getError(comptime T: type, err: *Matrix(T), F: Matrix(T), P: Matrix(T), S: Matrix(T), gpa: Allocator) !void {
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

    mm(T, &FP, F, P, 1, 0, false, false);

    var FPS = try Matrix(T).init(nbf, nbf, gpa);
    defer FPS.deinit(gpa);

    mm(T, &FPS, FP, S, 1, 0, false, false);

    for (0..nbf) |i| for (0..nbf) |j| {
        err.ptr(i, j).* = FPS.at(i, j) - FPS.at(j, i);
    };
}

/// Constructs the Fock (or Kohn-Sham) matrix from core Hamiltonian and two-electron interactions.
fn getFock(comptime T: type, F: *Matrix(T), ints: Integrals(T), P: Matrix(T), opt: Options, dft: ?*DftPotential(T), gpa: Allocator) !void {
    const nbf = ints.H.?.shape[0];

    std.debug.assert(F.shape[0] == nbf);
    std.debug.assert(F.shape[1] == nbf);
    std.debug.assert(P.shape[0] == nbf);
    std.debug.assert(P.shape[1] == nbf);

    for (0..nbf) |i| for (0..nbf) |j| {
        F.ptr(i, j).* = ints.H.?.at(i, j);
    };

    if (opt.integral_direct) {
        var exch_factor: T = if (opt.generalized) 1 else 0.5;

        if (dft) |pot| {
            exch_factor *= pot.exx_coef;
        }

        if (opt.generalized) {
            ints.sys.fockGhf(F, P, exch_factor, opt.nthreads);
        }

        if (!opt.generalized) {
            ints.sys.fockRhf(F, P, exch_factor, opt.nthreads);
        }

        if (dft) |pot| {
            try pot.evaluate(ints.sys, P, gpa);

            for (0..nbf) |j| for (0..nbf) |k| {
                F.ptr(j, k).* += pot.Vxc.at(j, k);
            };
        }
    }

    if (!opt.integral_direct) {
        std.debug.assert(ints.g.?.shape[0] == nbf);
        std.debug.assert(ints.g.?.shape[1] == nbf);
        std.debug.assert(ints.g.?.shape[2] == nbf);
        std.debug.assert(ints.g.?.shape[3] == nbf);

        if (dft) |pot| {
            for (0..nbf) |j| for (0..nbf) |k| for (0..nbf) |l| for (0..nbf) |m| {
                F.ptr(l, m).* += P.at(j, k) * ints.g.?.at(.{ j, l, k, m });
            };

            if (pot.exx_coef > 0) {
                const factor: T = if (opt.generalized) 1 else 0.5;

                for (0..nbf) |i| for (0..nbf) |j| for (0..nbf) |k| for (0..nbf) |l| {
                    F.ptr(k, l).* -= factor * pot.exx_coef * P.at(i, j) * ints.g.?.at(.{ i, j, k, l });
                };
            }

            try pot.evaluate(ints.sys, P, gpa);

            for (0..nbf) |j| for (0..nbf) |k| {
                F.ptr(j, k).* += pot.Vxc.at(j, k);
            };
        }

        if (dft == null) {
            const exch_factor: T = if (opt.generalized) 1 else 0.5;

            for (0..nbf) |i| for (0..nbf) |j| for (0..nbf) |k| for (0..nbf) |l| {
                F.ptr(k, l).* += P.at(i, j) * (ints.g.?.at(.{ i, k, j, l }) - exch_factor * ints.g.?.at(.{ i, j, k, l }));
            };
        }
    }

    for (0..nbf) |i| for (i + 1..nbf) |j| {
        const avg = (F.at(i, j) + F.at(j, i)) / 2;

        F.ptr(i, j).* = avg;
        F.ptr(j, i).* = avg;
    };
}

/// Computes the nuclear Hessian and performs harmonic frequency and thermochemical analyses.
fn handleHessianAndFrequencies(comptime T: type, io: std.Io, opt: Options, runFn: anytype, sys: *MolecularSystem(T), energy: T, log: bool, gpa: Allocator) ![]Matrix(T) {
    var hess = try gpa.alloc(Matrix(T), if (opt.hessian) |_| 1 else 0);
    errdefer if (opt.hessian) |_| gpa.free(hess);

    if (opt.hessian) |hessopt| switch (hessopt) {
        .numeric => hess[0] = try calculateNumericalHessian(T, io, runFn, opt, sys, log, gpa),
    };

    errdefer if (opt.hessian) |_| hess[0].deinit(gpa);

    if (log and opt.frequency != null) {
        var freqs = try calculateHarmonicFrequencies(T, hess[0], sys.*, gpa);
        defer freqs.deinit(gpa);

        const method_str = if (opt.dft != null) "DFT" else "HARTREE-FOCK";

        try printHarmonicFrequencies(T, io, freqs, method_str);

        const thermo = try calculateThermochemistry(T, opt.frequency.?, sys.*, freqs, opt.multiplicity, gpa);

        try printThermochemistry(T, io, opt.frequency.?, thermo, energy, method_str);
    }

    return hess;
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

        if (dft) |pot| {
            try getFock(T, ws.F, ints, ws.P.*, opt, pot, gpa);

            e_new = getEnergy(T, ints, ws.F.*, ws.P.*, pot) + VN;
        }

        if (dft == null) {
            try getFock(T, ws.F, ints, ws.P.*, opt, null, gpa);

            e_new = getEnergy(T, ints, ws.F.*, ws.P.*, null) + VN;
        }

        if (opt.diis != null and opt.diis.? > 0) {
            try fck_hist.ensureUnusedCapacity(gpa, 1);
            try err_hist.ensureUnusedCapacity(gpa, 1);

            var f_diis = try Matrix(T).init(nbf, nbf, gpa);
            errdefer f_diis.deinit(gpa);

            var e_diis = try Matrix(T).init(nbf, nbf, gpa);
            errdefer e_diis.deinit(gpa);

            for (0..nbf) |j| for (0..nbf) |k| {
                f_diis.ptr(j, k).* = ws.F.at(j, k);
            };

            try getError(T, &e_diis, ws.F.*, ws.P.*, ints.S.?, gpa);

            if (fck_hist.items.len >= opt.diis.?) {
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

                for (0..nbf) |j| for (0..nbf) |k| {
                    f_retry.ptr(j, k).* = ws.F.at(j, k);
                };

                try getError(T, &e_retry, ws.F.*, ws.P.*, ints.S.?, gpa);

                fck_hist.appendAssumeCapacity(f_retry);
                err_hist.appendAssumeCapacity(e_retry);
            };
        }

        @memcpy(B.data, ints.S.?.data);

        try geigh(T, ws.e, ws.C, ws.F.*, &B);

        const delta_energy, const p_rm = .{ @abs(e_new - e_old), getDensity(T, ws.P, ws.C.*, nocc, opt.generalized) };

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
