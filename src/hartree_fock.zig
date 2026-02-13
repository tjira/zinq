//! Hartree-Fock self-consistent field method implementation.

const std = @import("std");

const basis_set = @import("basis_set.zig");
const classical_particle = @import("classical_particle.zig");
const contracted_gaussian = @import("contracted_gaussian.zig");
const device_write = @import("device_write.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const energy_derivative = @import("energy_derivative.zig");
const errror_handling = @import("error_handling.zig");
const frequency_analysis = @import("frequency_analysis.zig");
const global_variables = @import("global_variables.zig");
const integral_transform = @import("integral_transform.zig");
const linear_solve = @import("linear_solve.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const molecular_integrals = @import("molecular_integrals.zig");
const object_array = @import("object_array.zig");
const particle_optimization = @import("particle_optimization.zig");
const real_matrix = @import("real_matrix.zig");
const real_tensor_four = @import("real_tensor_four.zig");
const real_vector = @import("real_vector.zig");

const BasisSet = basis_set.BasisSet;
const ClassicalParticle = classical_particle.ClassicalParticle;
const ContractedGaussian = contracted_gaussian.ContractedGaussian;
const RealMatrix = real_matrix.RealMatrix;
const RealMatrixArray = object_array.RealMatrixArray;
const RealTensor4 = real_tensor_four.RealTensor4;
const RealVector = real_vector.RealVector;

const coulomb = molecular_integrals.coulomb;
const eigensystemHermitian = eigenproblem_solver.eigensystemHermitian;
const eigensystemHermitianAlloc = eigenproblem_solver.eigensystemHermitianAlloc;
const exportRealMatrix = device_write.exportRealMatrix;
const exportRealTensorFour = device_write.exportRealTensorFour;
const kinetic = molecular_integrals.kinetic;
const linearSolveSymmetric = linear_solve.linearSolveSymmetric;
const mm = matrix_multiplication.mm;
const mmAlloc = matrix_multiplication.mmAlloc;
const nuclear = molecular_integrals.nuclear;
const nuclearGradient = energy_derivative.nuclearGradient;
const nuclearHessian = energy_derivative.nuclearHessian;
const oneAO2AS = integral_transform.oneAO2AS;
const overlap = molecular_integrals.overlap;
const particleHarmonicFrequencies = frequency_analysis.particleHarmonicFrequencies;
const particleSteepestDescent = particle_optimization.particleSteepestDescent;
const print = device_write.print;
const printClassicalParticleAsMolecule = device_write.printClassicalParticleAsMolecule;
const printJson = device_write.printJson;
const printRealMatrix = device_write.printRealMatrix;
const throw = errror_handling.throw;
const twoAO2AS = integral_transform.twoAO2AS;

const MAX_POOL_SIZE = global_variables.MAX_POOL_SIZE;
const SINGULARITY_TOLERANCE = global_variables.SINGULARITY_TOLERANCE;
const TEST_TOLERANCE = global_variables.TEST_TOLERANCE;

/// Hartree-Fock target options.
pub fn Options(comptime T: type) type {
    return struct {
        const Diis = struct {
            size: u32 = 5,
            start: u32 = 1,
        };
        const Gradient = union(enum) {
            numeric: struct {
                step: T = 1e-3,
                nthread: u32 = 1
            },
            analytic: struct {}
        };
        const Hessian = union(enum) {
            numeric: struct {
                step: T = 1e-3,
                nthread: u32 = 1
            },
            analytic: struct {}
        };
        const Optimize = struct {
            maxiter: u32 = 100,
            step: T = 1,
            threshold: T = 1e-6,
        };
        const Write = struct {
            coefficient: ?[]const u8 = null,
            density: ?[]const u8 = null,
            fock: ?[]const u8 = null,
        };

        system: []const u8,
        basis: []const u8,

        charge: i32 = 0,
        generalized: bool = false,
        maxiter: u32 = 100,
        nthread: u32 = 1,
        threshold: T = 1e-12,

        diis: ?Diis = .{},
        gradient: ?Gradient = null,
        hessian: ?Hessian = null,
        optimize: ?Optimize = null,
        write: Write = .{}
    };
}

/// Hartree-Fock target output.
pub fn Output(comptime T: type) type {
    return struct {
        C: RealMatrix(T), F: RealMatrix(T), J: RealTensor4(T), K: RealMatrix(T), P: RealMatrix(T),
        S: RealMatrix(T), V: RealMatrix(T), G: ?RealMatrix(T) = null, H: ?RealMatrix(T) = null,
        energy: T, epsilon: RealMatrix(T), frequencies: ?RealVector(T) = null,

        /// Deinitialize the output struct.
        pub fn deinit(self: @This(), allocator: std.mem.Allocator) void {
            self.C.deinit(allocator); self.F.deinit(allocator); self.K.deinit(allocator); self.P.deinit(allocator);
            self.S.deinit(allocator); self.V.deinit(allocator); self.epsilon.deinit(allocator); self.J.deinit(allocator);
            if (self.G != null) self.G.?.deinit(allocator);
            if (self.H != null) self.H.?.deinit(allocator);
            if (self.frequencies != null) self.frequencies.?.deinit(allocator);
        }
    };
}

/// Run the Hartree-Fock target.
pub fn run(comptime T: type, opt: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(opt);
    
    if (opt.gradient != null and opt.gradient.? == .analytic) return throw(Output(T), "ANALYTIC GRADIENT NOT IMPLEMENTED", .{});
    if (opt.hessian != null and opt.hessian.? == .analytic) return throw(Output(T), "ANALYTIC HESSIAN NOT IMPLEMENTED", .{});

    var system = try classical_particle.read(T, opt.system, opt.charge, allocator); defer system.deinit(allocator);

    if (enable_printing) {try print("\nINPUT GEOMETRY (Å):\n", .{}); try printClassicalParticleAsMolecule(T, system, null);}

    if (opt.optimize != null) {

        const optimized_system = try particleSteepestDescent(T, opt, system, scf, "HARTREE-FOCK", enable_printing, allocator);

        system.deinit(allocator); system = optimized_system;
    }

    if (enable_printing and opt.optimize != null) {try print("\nOPTIMIZED GEOMETRY (Å):\n", .{}); try printClassicalParticleAsMolecule(T, system, null);}

    var output = try scf(T, opt, system, enable_printing, allocator); errdefer output.deinit(allocator);

    output.G = if (opt.gradient != null) try nuclearGradient(T, opt, system, scf, "HARTREE-FOCK", enable_printing, allocator) else null;

    if (output.G) |G| {try print("\nHARTREE-FOCK NUCLEAR GRADIENT (Eh/Bohr):\n", .{}); try printRealMatrix(T, G);}

    output.H = if (opt.hessian != null) try nuclearHessian(T, opt, system, scf, "HARTREE-FOCK", enable_printing, allocator) else null;

    if (output.H) |H| output.frequencies = try particleHarmonicFrequencies(T, system, H, allocator);

    if (output.frequencies) |freqs| {try print("\nHARTREE-FOCK VIBRATIONAL FREQUENCIES (CM^-1):\n", .{}); try printRealMatrix(T, freqs.asMatrix());}

    if (opt.write.coefficient) |path| try exportRealMatrix(T, path, output.C);
    if (opt.write.density) |path| try exportRealMatrix(T, path, output.P);
    if (opt.write.fock) |path| try exportRealMatrix(T, path, output.F);

    return output;
}

/// Calculates the energy from the density matrix, Fock matrix and core Hamiltonian.
pub fn calculateEnergy(comptime T: type, K: RealMatrix(T), V: RealMatrix(T), F: RealMatrix(T), P: RealMatrix(T), generalized: bool) T {
    var energy: T = 0;

    const factor: T = if (generalized) 0.5 else 1;

    for (0..P.rows) |i| for (0..P.cols) |j| {
        energy += factor * P.at(i, j) * (K.at(i, j) + V.at(i, j) + F.at(i, j));
    };

    return energy;
}

/// Extrapolate the DIIS error to obtain a new Fock matrix. The error vector from the first iteration is ignored, since it is zero.
pub fn diisExtrapolate(comptime T: type, F: *RealMatrix(T), DIIS_F: RealMatrixArray(T), DIIS_E: RealMatrixArray(T), iter: usize, allocator: std.mem.Allocator) !void {
    const size = @min(DIIS_F.len, iter);

    var A = try RealMatrix(T).init(size + 1, size + 1, allocator); defer A.deinit(allocator);
    var b = try RealVector(T).init(size + 1,           allocator); defer b.deinit(allocator);
    var c = try RealVector(T).init(size + 1,           allocator); defer c.deinit(allocator);

    var temporary = try RealVector(T).init(c.len, allocator); defer temporary.deinit(allocator);

    A.fill(1); b.fill(0); A.ptr(A.rows - 1, A.cols - 1).* = 0; b.ptr(b.len - 1).* = 1;

    for (0..size) |i| for (i..size) |j| {

        A.ptr(i, j).* = 0;

        const ii = (iter - size + i + 1) % DIIS_E.len;
        const jj = (iter - size + j + 1) % DIIS_E.len;

        for (0..DIIS_E.at(0).rows) |k| for (0..DIIS_E.at(0).cols) |l| {
            A.ptr(i, j).* += DIIS_E.at(ii).at(k, l) * DIIS_E.at(jj).at(k, l);
        };

        A.ptr(j, i).* = A.at(i, j);
    };

    const AJC = try eigensystemHermitianAlloc(T, A, allocator); defer AJC.J.deinit(allocator); defer AJC.C.deinit(allocator);

    for (0..b.len) |i| if (@abs(AJC.J.at(i, i)) < SINGULARITY_TOLERANCE) return;

    try linearSolveSymmetric(T, &c, A, AJC.J, AJC.C, b, &temporary); F.zero();

    for (0..size) |i| {

        const ii = (iter - size + i + 1) % DIIS_E.len;

        for (0..F.rows) |j| for (0..F.cols) |k| {
            F.ptr(j, k).* += c.at(i) * DIIS_F.at(ii).at(j, k);
        };
    }

    try F.symmetrize();
}

/// Function to calculate the error vector.
pub fn errorVector(comptime T: type, S: RealMatrix(T), F: RealMatrix(T), P: RealMatrix(T), allocator: std.mem.Allocator) !RealMatrix(T) {
    const SP = try mmAlloc(T, S, false, P, false, allocator); defer SP.deinit(allocator);
    const FP = try mmAlloc(T, F, false, P, false, allocator); defer FP.deinit(allocator);

    const SPF = try mmAlloc(T, SP, false, F, false, allocator); defer SPF.deinit(allocator);
    const FPS = try mmAlloc(T, FP, false, S, false, allocator); defer FPS.deinit(allocator);

    var e = try RealMatrix(T).initZero(P.rows, P.cols, allocator);

    for (0..e.rows) |i| for (i + 1..e.cols) |j| {
        e.ptr(i, j).* = SPF.at(i, j) - FPS.at(i, j); e.ptr(j, i).* = -e.at(i, j);
    };

    return e;
}

/// Calculate the density matrix.
pub fn getDensityMatrix(comptime T: type, P: *RealMatrix(T), C: RealMatrix(T), nocc: usize) void {
    for (0..P.rows) |i| for (i..P.cols) |j| {

        P.ptr(i, j).* = 0;

        for (0..nocc) |m| {
            P.ptr(i, j).* += C.at(i, m) * C.at(j, m);
        }

        P.ptr(j, i).* = P.at(i, j);
    };
}

/// Obtain the Fock matrix form core Hamiltonian and density matrix.
pub fn getFockMatrix(comptime T: type, F: *RealMatrix(T), K: RealMatrix(T), V: RealMatrix(T), P: RealMatrix(T), J: RealTensor4(T), basis: BasisSet(T)) !void {
    for (0..F.rows) |i| for (0..F.cols) |j| {
        F.ptr(i, j).* = V.at(i, j) + K.at(i, j);
    };

    const factor: T = if (P.rows == 2 * basis.nbf()) 1 else 2;

    for (0..J.shape[0]) |i| for (0..J.shape[1]) |j| for (0..J.shape[2]) |k| for (0..J.shape[3]) |l| {
        F.ptr(k, l).* += P.at(i, j) * (factor * J.at(i, j, k, l) - J.at(i, l, k, j));
    };

    try F.symmetrize();
}

/// Function to get the X matrix.
pub fn getXMatrix(comptime T: type, S: RealMatrix(T), allocator: std.mem.Allocator) !RealMatrix(T) {
    var X = try RealMatrix(T).init(S.rows, S.cols, allocator);

    var S_C = try RealMatrix(T).init(S.rows, S.cols, allocator); defer S_C.deinit(allocator);

    var mm_temp = try RealMatrix(T).init(S.rows, S.cols, allocator); defer mm_temp.deinit(allocator);

    try eigensystemHermitian(T, &X, &S_C, S);

    for (0..S.rows) |i| X.ptr(i, i).* = 1.0 / std.math.sqrt(X.at(i, i));

    try mm(T, &mm_temp, X, false, S_C, true);
    try mm(T, &X, S_C, false, mm_temp, false);

    return X;
}

/// Perform the SCF procedure and return the output.
pub fn scf(comptime T: type, opt: Options(T), system: ClassicalParticle(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    var basis = try BasisSet(T).init(system, opt.basis, allocator); defer basis.deinit(allocator);

    const nbf = if (opt.generalized) 2 * basis.nbf() else basis.nbf();
    const nocc = if (opt.generalized) try system.noccSpin() else try system.noccSpatial();

    if (enable_printing) try print("\nNUMBER OF BASIS FUNCTIONS: {d}\n", .{nbf});

    if (enable_printing) try print("\nONE-ELECTRON INTEGRALS: ", .{});

    var timer = try std.time.Timer.start();

    var S = try overlap(T, basis, opt.nthread, allocator);
    var K = try kinetic(T, basis, opt.nthread, allocator);
    var V = try nuclear(T, system, basis, opt.nthread, allocator);

    if (opt.generalized) {
        try oneAO2AS(T, &S, allocator);
        try oneAO2AS(T, &K, allocator);
        try oneAO2AS(T, &V, allocator);
    }

    if (enable_printing) try print("{D}\n", .{timer.read()}); timer.reset();

    if (enable_printing) try print("TWO-ELECTRON INTEGRALS: ", .{});

    var J = try coulomb(T, basis, opt.nthread, allocator);

    if (opt.generalized) {
        try twoAO2AS(T, &J, allocator);
    }

    if (enable_printing) try print("{D}\n", .{timer.read()});

    var C = try RealMatrix(T).initZero(nbf, nbf, allocator);
    var F = try RealMatrix(T).initZero(nbf, nbf, allocator);
    var P = try RealMatrix(T).initZero(nbf, nbf, allocator);

    var epsilon = try RealMatrix(T).initZero(nbf, nbf, allocator);

    var X = try getXMatrix(T, S, allocator); defer X.deinit(allocator);

    var DIIS_F = try RealMatrixArray(T).initZero(if (opt.diis != null) opt.diis.?.size else 0, .{.rows = nbf, .cols = nbf}, allocator); defer DIIS_F.deinit(allocator);
    var DIIS_E = try RealMatrixArray(T).initZero(if (opt.diis != null) opt.diis.?.size else 0, .{.rows = nbf, .cols = nbf}, allocator); defer DIIS_E.deinit(allocator);

    const VNN = system.nuclearRepulsionEnergy();

    var energy: T = 0; var energy_prev: T = 1; var iter: usize = 0;

    if (enable_printing) try print("\nSELF CONSISTENT FIELD:\n{s:4} {s:20} {s:8} {s:4}\n", .{"ITER", "ENERGY", "|DELTA E|", "TIME"});

    while (@abs(energy - energy_prev) > opt.threshold) : (iter += 1) {

        if (iter >= opt.maxiter) return throw(Output(T), "HARTREE-FOCK DID NOT CONVERGE IN {d} ITERATIONS", .{opt.maxiter});

        timer.reset();

        try getFockMatrix(T, &F, K, V, P, J, basis);

        if (opt.diis != null) {

            const e = try errorVector(T, S, F, P, allocator); defer e.deinit(allocator);

            try F.copyTo(DIIS_F.ptr(iter % DIIS_F.len));
            try e.copyTo(DIIS_E.ptr(iter % DIIS_E.len));

            if (iter >= opt.diis.?.start) try diisExtrapolate(T, &F, DIIS_F, DIIS_E, iter, allocator);
        }

        try solveRoothaan(T, &epsilon, &C, F, X, allocator);

        getDensityMatrix(T, &P, C, nocc);

        energy_prev = energy; energy = calculateEnergy(T, K, V, F, P, opt.generalized);

        if (enable_printing) try print("{d:4} {d:20.14} {e:9.3} {D}\n", .{iter + 1, energy + VNN, @abs(energy - energy_prev), timer.read()});
    }

    if (enable_printing) try print("\nHF ENERGY: {d:.14}\n", .{energy + VNN});

    return .{.C = C, .F = F, .J = J, .K = K, .P = P, .S = S, .V = V, .energy = energy + VNN, .epsilon = epsilon};
}

/// Solver for the Roothaan equations.
pub fn solveRoothaan(comptime T: type, E: *RealMatrix(T), C: *RealMatrix(T), F: RealMatrix(T), X: RealMatrix(T), allocator: std.mem.Allocator) !void {
    const FX = try mmAlloc(T, F, false, X, false, allocator); defer FX.deinit(allocator);
    var XFX = try mmAlloc(T, X, false, FX, false, allocator); defer XFX.deinit(allocator);

    try XFX.symmetrize();

    const XFXJC = try eigensystemHermitianAlloc(T, XFX, allocator); defer XFXJC.J.deinit(allocator); defer XFXJC.C.deinit(allocator);

    try mm(T, C, X, false, XFXJC.C, false); try XFXJC.J.copyTo(E);
}

test "Hartree-Fock Calculation for a Water Molecule with STO-3G Basis Set" {
    const opt = Options(f64){
        .system = "example/molecule/water.xyz",
        .basis = "sto-3g"
    };

    var output = try run(f64, opt, false, std.testing.allocator); defer output.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(output.energy, -74.96590121728507, TEST_TOLERANCE);
}

test "Hartree-Fock Calculation for a Methane Molecule with 6-31G* Basis Set" {
    const opt = Options(f64){
        .system = "example/molecule/methane.xyz",
        .basis = "6-31g*"
    };

    var output = try run(f64, opt, false, std.testing.allocator); defer output.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(output.energy, -40.19517074914403, TEST_TOLERANCE);
}

test "Generalized Hartree-Fock Calculation for a Water Molecule with STO-3G Basis Set" {
    const opt = Options(f64){
        .system = "example/molecule/water.xyz",
        .basis = "sto-3g",
        .generalized = true
    };

    var output = try run(f64, opt, false, std.testing.allocator); defer output.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(output.energy, -74.96590121728507, TEST_TOLERANCE);
}

test "Generalized Hartree-Fock Calculation for a Methane Molecule with 6-31G* Basis Set" {
    const opt = Options(f64){
        .system = "example/molecule/methane.xyz",
        .basis = "6-31g*",
        .generalized = true
    };

    var output = try run(f64, opt, false, std.testing.allocator); defer output.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(output.energy, -40.19517074914403, TEST_TOLERANCE);
}
