//! Module to perform the configuration interaction (CI) method for electronic structure calculations.

const std = @import("std");

const classical_particle = @import("classical_particle.zig");
const complete_active_space = @import("complete_active_space.zig");
const device_write = @import("device_write.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const energy_derivative = @import("energy_derivative.zig");
const hartree_fock = @import("hartree_fock.zig");
const frequency_analysis = @import("frequency_analysis.zig");
const global_variables = @import("global_variables.zig");
const integral_transform = @import("integral_transform.zig");
const math_functions = @import("math_functions.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const molecular_integrals = @import("molecular_integrals.zig");
const particle_optimization = @import("particle_optimization.zig");
const real_matrix = @import("real_matrix.zig");
const real_tensor_four = @import("real_tensor_four.zig");
const real_vector = @import("real_vector.zig");

const ClassicalParticle = classical_particle.ClassicalParticle;
const RealMatrix = real_matrix.RealMatrix;
const RealTensor4 = real_tensor_four.RealTensor4;
const RealVector = real_vector.RealVector;

const eigensystemHermitian = eigenproblem_solver.eigensystemHermitian;
const generateCasDeterminants = complete_active_space.generateCasDeterminants;
const kinetic = molecular_integrals.kinetic;
const scf = hartree_fock.scf;
const nuclear = molecular_integrals.nuclear;
const nuclearGradient = energy_derivative.nuclearGradient;
const nuclearHessian = energy_derivative.nuclearHessian;
const oneAO2MO = integral_transform.oneAO2MO;
const oneAO2MS = integral_transform.oneAO2MS;
const overlap = molecular_integrals.overlap;
const particleHarmonicFrequencies = frequency_analysis.particleHarmonicFrequencies;
const particleSteepestDescent = particle_optimization.particleSteepestDescent;
const print = device_write.print;
const printClassicalParticleAsMolecule = device_write.printClassicalParticleAsMolecule;
const printJson = device_write.printJson;
const printRealMatrix = device_write.printRealMatrix;
const twoAO2MO = integral_transform.twoAO2MO;
const twoAO2MS = integral_transform.twoAO2MS;

const AU2EV = global_variables.AU2EV;
const TEST_TOLERANCE = global_variables.TEST_TOLERANCE;

/// The Configuration Interaction (CI) method for electronic structure calculations.
pub fn Options(comptime T: type) type {
    return struct {
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

        active_space: ?[2]u32 = null,
        states: u32 = 5,

        hartree_fock: hartree_fock.Options(T),

        gradient: ?Gradient = null,
        hessian: ?Hessian = null,
        optimize: ?Optimize = null
    };
}

/// Configuration interaction output struct containing the energy, gradient, hessian, frequencies, and the underlying Hartree-Fock output.
pub fn Output(comptime T: type) type {
    return struct {
        G: ?RealMatrix(T) = null, H: ?RealMatrix(T) = null,
        energy: T, frequencies: ?RealVector(T) = null,
        hf_output: hartree_fock.Output(T),

        /// Deinitialize the output struct.
        pub fn deinit(self: @This(), allocator: std.mem.Allocator) void {
            if (self.G != null) self.G.?.deinit(allocator);
            if (self.H != null) self.H.?.deinit(allocator);
            if (self.frequencies != null) self.frequencies.?.deinit(allocator);
            self.hf_output.deinit(allocator);
        }
    };
}

/// Run the CI calculation with the given options for the provided system, returning the energy, gradient, hessian, and frequencies if requested.
pub fn run(comptime T: type, opt: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (opt.gradient != null and opt.gradient.? == .analytic) {

        std.log.err("ANALYTIC GRADIENT NOT IMPLEMENTED FOR CONFIGURATION INTERACTION", .{});

        return error.InvalidInput;
    }

    if (opt.hessian != null and opt.hessian.? == .analytic) {

        std.log.err("ANALYTIC HESSIAN NOT IMPLEMENTED FOR CONFIGURATION INTERACTION METHOD", .{});

        return error.InvalidInput;
    }

    var system = try classical_particle.read(T, opt.hartree_fock.system, opt.hartree_fock.charge, 0, allocator); defer system.deinit(allocator);

    if (enable_printing) {try print("\nINPUT GEOMETRY (A):\n", .{}); try printClassicalParticleAsMolecule(T, system, null);}

    if (opt.optimize != null) {

        const optimized_system = try particleSteepestDescent(T, opt, system, ci, "CONFIGURATION INTERACTION", enable_printing, allocator);

        system.deinit(allocator); system = optimized_system;
    }

    if (enable_printing and opt.optimize != null) {try print("\nOPTIMIZED GEOMETRY (A):\n", .{}); try printClassicalParticleAsMolecule(T, system, null);}

    var output = try ci(T, opt, system, enable_printing, allocator); errdefer output.deinit(allocator);

    output.G = if (opt.gradient != null) try nuclearGradient(T, opt, system, ci, "CONFIGURATION INTERACTION", enable_printing, allocator) else null;

    if (output.G) |G| {try print("\nCONFIGURATION INTERACTION NUCLEAR GRADIENT (Eh/Bohr):\n", .{}); try printRealMatrix(T, G);}

    output.H = if (opt.hessian != null) try nuclearHessian(T, opt, system, ci, "CONFIGURATION INTERACTION", enable_printing, allocator) else null;

    if (output.H) |H| output.frequencies = try particleHarmonicFrequencies(T, system, H, allocator);

    if (output.frequencies) |freqs| {try print("\nCONFIGURATION INTERACTION VIBRATIONAL FREQUENCIES (CM^-1):\n", .{}); try printRealMatrix(T, freqs.asMatrix());}

    return output;
}

/// Primary function to run the CI calculation, which first performs a Hartree-Fock calculation.
pub fn ci(comptime T: type, opt: Options(T), system: ClassicalParticle(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {

    const hf_output = try scf(T, opt.hartree_fock, system, enable_printing, allocator);

    const nbf = if (opt.hartree_fock.generalized) hf_output.S.rows else 2 * hf_output.S.rows; const nocc = try system.noccSpin();

    var AS = [2]usize{nocc, nbf}; if (opt.active_space) |as| {AS[0] = as[0]; AS[1] = as[1];}

    if (AS[0] > nocc) {

        std.log.err("ACTIVE SPACE CANNOT HAVE MORE ELECTRONS THAN OCCUPIED SPINORBITALS", .{});

        return error.InvalidInput;
    }

    if (AS[1] > nbf) {

        std.log.err("ACTIVE SPACE CANNOT HAVE MORE SPINORBITALS THAN BASIS FUNCTIONS", .{});

        return error.InvalidInput;
    }

    if (AS[1] - AS[0] > nbf - nocc) {

        std.log.err("ACTIVE SPACE CANNOT HAVE MORE VIRTUAL SPINORBITALS THAN AVAILABLE", .{});

        return error.InvalidInput;
    }

    if (enable_printing) try print("\nCI ACTIVE SPACE: {d} ELECTRONS IN {d} SPINORBITALS\n", .{AS[0], AS[1]});

    var D = try generateCasDeterminants(AS[0], AS[1], nocc, opt.hartree_fock.generalized, allocator); defer D.deinit(allocator);

    var H_AO = try hf_output.K.clone(allocator); defer H_AO.deinit(allocator); try H_AO.add(hf_output.V);

    var H_MS = try RealMatrix (T).init(nbf, nbf,                     allocator); defer H_MS.deinit(allocator);
    var J_MS = try RealTensor4(T).init([_]usize{nbf, nbf, nbf, nbf}, allocator); defer J_MS.deinit(allocator);

    try transform(T, &H_MS, &J_MS, H_AO, hf_output.J, hf_output.C);

    if (enable_printing) try print("\nNUMBER OF CI DETERMINANTS: {d}\n", .{D.rows});

    var A = try RealVector(usize).init(nocc,           allocator); defer A.deinit(allocator);
    var H = try RealMatrix(T    ).init(D.rows, D.rows, allocator); defer H.deinit(allocator);
    var E = try RealMatrix(T    ).init(D.rows, D.rows, allocator); defer E.deinit(allocator);
    var C = try RealMatrix(T    ).init(D.rows, D.rows, allocator); defer C.deinit(allocator);

    H.zero();

    for (0..D.rows) |i| for (i..D.rows) |j| {

        var so: [4]usize = undefined; var diff: u32 = 0; var k: usize = 0;

        const sign = try alignDeterminant(&A, D.row(i), D.row(j));

        for (0..nocc) |l| if (A.at(l) != D.at(i, l)) {diff += 1;};

        if (diff > 2) continue;

        if (diff == 1) for (0..nocc) |l| if (A.at(l) != D.at(i, l)) {
            so[0] = A.at(l); so[1] = D.at(i, l);
        };

        if (diff == 2) for (0..nocc) |l| if (A.at(l) != D.at(i, l)) {
            so[if (k == 0) 0 else 1] = A.at(l); so[if (k == 0) 2 else 3] = D.at(i, l); k += 1;
        };

        H.ptr(i, j).* = @as(T, @floatFromInt(sign)) * try slater(T, A, so[0..diff * 2], H_MS, J_MS); H.ptr(j, i).* = H.at(i, j);
    };
    
    try eigensystemHermitian(T, &E, &C, H); const energy = E.at(0, 0) + system.nuclearRepulsionEnergy();

    if (enable_printing) try print("\n", .{});

    if (enable_printing) {

        for (0..@min(opt.states, E.rows) + 1) |i| {

            const deltaE = E.at(i, i) - E.at(0, 0);

            try print("CASCI STATE {d:2}: {d:20.14} Eh (DELTA E = {d:7.4} Eh = {d:7.4} eV)\n", .{i, E.at(i, i) + system.nuclearRepulsionEnergy(), deltaE, deltaE * AU2EV});
        }
    }

    return .{
        .hf_output = hf_output, .energy = energy, .G = null, .H = null, .frequencies = null,
    };
}

/// Aligns the vector C to the vector B. The result is stored in the vector A and the sign of the permutation is returned.
fn alignDeterminant(A: *RealVector(usize), B: RealVector(usize), C: RealVector(usize)) !i32 {
    try C.copyTo(A); var k: i32 = 0; var sign: i32 = 1;

    while (k < A.len) : (k += 1) {
        if (A.at(@intCast(k)) != B.at(@intCast(k))) for (@as(usize, @intCast(k)) + 1..A.len) |l| if (A.at(@intCast(k)) == B.at(l) or A.at(l) == B.at(@intCast(k))) {
            std.mem.swap(usize, A.ptr(@intCast(k)), A.ptr(l)); sign *= -1; k -= 1; break;
        };
    }

    return sign;
}

/// Slater-Condon rules for the CI calculations.
fn slater(comptime T: type, A: RealVector(usize), so: []const usize, H_MS: RealMatrix(T), J_MS: RealTensor4(T)) !T {
    var hij: T = 0;

    if (so.len / 2 == 0) {

        for (0..A.len) |l| {
            hij += H_MS.at(A.at(l), A.at(l));
        }

        for (0..A.len) |l| for (0..A.len) |m| {

            const J_MS_A = J_MS.at(A.at(l), A.at(l), A.at(m), A.at(m)) - J_MS.at(A.at(l), A.at(m), A.at(m), A.at(l));

            hij += 0.5 * J_MS_A;
        };
    }

    else if (so.len / 2 == 1) {

        hij += H_MS.at(so[0], so[1]);

        for (0..A.len) |m| if (A.at(m) != so[0]) {

            const J_MS_A = J_MS.at(so[0], so[1], A.at(m), A.at(m)) - J_MS.at(so[0], A.at(m), A.at(m), so[1]);

            hij += J_MS_A;
        };
    }

    else if (so.len / 2 == 2) {

        const J_MS_A = J_MS.at(so[0], so[2], so[1], so[3]) - J_MS.at(so[0], so[3], so[1], so[2]);

        hij = J_MS_A;
    }
    
    return hij;
}

/// Function to perform all integrals transformations used in the CI method.
pub fn transform(comptime T: type, H_MS: *RealMatrix(T), J_MS: *RealTensor4(T), H: RealMatrix(T), J: RealTensor4(T), C: RealMatrix(T)) !void {
    if (H.rows != H_MS.rows) {oneAO2MS(T, H_MS, H, C);} else {oneAO2MO(T, H_MS, H, C);}
    if (H.rows != H_MS.rows) {twoAO2MS(T, J_MS, J, C);} else {twoAO2MO(T, J_MS, J, C);}
}
