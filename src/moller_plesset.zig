//! Module to perform the Moller-Plesset calculations.

const std = @import("std");

const basis_set = @import("basis_set.zig");
const classical_particle = @import("classical_particle.zig");
const contracted_gaussian = @import("contracted_gaussian.zig");
const device_write = @import("device_write.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const energy_derivative = @import("energy_derivative.zig");
const errror_handling = @import("error_handling.zig");
const hartree_fock = @import("hartree_fock.zig");
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
const eigensystemSymmetric = eigenproblem_solver.eigensystemSymmetric;
const eigensystemSymmetricAlloc = eigenproblem_solver.eigensystemSymmetricAlloc;
const exportRealMatrix = device_write.exportRealMatrix;
const exportRealTensorFour = device_write.exportRealTensorFour;
const kinetic = molecular_integrals.kinetic;
const linearSolveSymmetric = linear_solve.linearSolveSymmetric;
const mm = matrix_multiplication.mm;
const mmAlloc = matrix_multiplication.mmAlloc;
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
const throw = errror_handling.throw;
const twoAO2MO = integral_transform.twoAO2MO;
const twoAO2MS = integral_transform.twoAO2MS;

const MAX_POOL_SIZE = global_variables.MAX_POOL_SIZE;
const TEST_TOLERANCE = global_variables.TEST_TOLERANCE;

/// The Moller-Plesset options
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

        order: u32 = 2,

        hartree_fock: hartree_fock.Options(T),

        gradient: ?Gradient = null,
        hessian: ?Hessian = null,
        optimize: ?Optimize = null
    };
}

/// MollerPlesset target output.
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

/// Run the Moller-Plesset calculation with the given options.
pub fn run(comptime T: type, opt: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (opt.gradient != null and opt.gradient.? == .analytic) return throw(Output(T), "ANALYTIC GRADIENT NOT IMPLEMENTED", .{});
    if (opt.hessian != null and opt.hessian.? == .analytic) return throw(Output(T), "ANALYTIC HESSIAN NOT IMPLEMENTED", .{});

    var system = try classical_particle.read(T, opt.hartree_fock.system, opt.hartree_fock.charge, allocator); defer system.deinit(allocator);

    if (enable_printing) {try print("\nINPUT GEOMETRY (Å):\n", .{}); try printClassicalParticleAsMolecule(T, system, null);}

    if (opt.optimize != null) {

        const optimized_system = try particleSteepestDescent(T, opt, system, mp, "MOLLER PLESSET", enable_printing, allocator);

        system.deinit(allocator); system = optimized_system;
    }

    if (enable_printing and opt.optimize != null) {try print("\nOPTIMIZED GEOMETRY (Å):\n", .{}); try printClassicalParticleAsMolecule(T, system, null);}

    var output = try mp(T, opt, system, enable_printing, allocator); errdefer output.deinit(allocator);

    output.G = if (opt.gradient != null) try nuclearGradient(T, opt, system, mp, "MOLLER-PLESSET", enable_printing, allocator) else null;

    if (output.G) |G| {try print("\nMOLLER-PLESSET NUCLEAR GRADIENT (Eh/Bohr):\n", .{}); try printRealMatrix(T, G);}

    output.H = if (opt.hessian != null) try nuclearHessian(T, opt, system, mp, "MOLLER-PLESSET", enable_printing, allocator) else null;

    if (output.H) |H| output.frequencies = try particleHarmonicFrequencies(T, system, H, allocator);

    if (output.frequencies) |freqs| {try print("\nMOLLER-PLESSET VIBRATIONAL FREQUENCIES (CM^-1):\n", .{}); try printRealMatrix(T, freqs.asMatrix());}

    return output;
}

/// Primary function to run the Moller-Plesset calculation with the given options for the provided system.
pub fn mp(comptime T: type, opt: Options(T), system: ClassicalParticle(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    const hf_output = try scf(T, opt.hartree_fock, system, enable_printing, allocator);

    const nbf = if (opt.hartree_fock.generalized) hf_output.S.rows else 2 * hf_output.S.rows; const nocc = try system.noccSpin();

    var F_MS = try RealMatrix (T).init(nbf, nbf,                     allocator); defer F_MS.deinit(allocator);
    var J_MS = try RealTensor4(T).init([_]usize{nbf, nbf, nbf, nbf}, allocator); defer J_MS.deinit(allocator);

    try transform(T, &F_MS, &J_MS, hf_output.F, hf_output.J, hf_output.C);

    var energy: T = 0;

    if (opt.order < 2) return throw(Output(T), "MOLLER-PLESSET ORDER MUST BE >= 2", .{});

    if (opt.order >= 2) energy += mp2(T, F_MS, J_MS, nocc);

    if (opt.order > 2) return throw(Output(T), "MOLLER-PLESSET OF ORDER HIGHER THAN 2 NOT IMPLEMENTED", .{});

    if (enable_printing) try print("\nMP{d} ENERGY: {d:.14}\n", .{opt.order, hf_output.energy + energy});

    return .{
        .hf_output = hf_output, .energy = hf_output.energy + energy, .G = null, .H = null, .frequencies = null,
    };
}

/// Returns the second-order Moller-Plesset energy.
pub fn mp2(comptime T: type, F_MS: RealMatrix(T), J_MS: RealTensor4(T), nocc: usize) T {
    var energy: T = 0;

    for (0..nocc) |i| for (0..nocc) |j| for (nocc..J_MS.shape[0]) |a| for (nocc..J_MS.shape[0]) |b| {

        const J_MS_A_ijab = J_MS.at(i, a, j, b) - J_MS.at(i, b, j, a);

        energy += 0.25 * J_MS_A_ijab * J_MS_A_ijab / (F_MS.at(i, i) + F_MS.at(j, j) - F_MS.at(a, a) - F_MS.at(b, b));
    };

    return energy;
}

/// Function to perform all integrals transformations used in the Moller-Plesset calculations.
pub fn transform(comptime T: type, F_MS: *RealMatrix(T), J_MS: *RealTensor4(T), F: RealMatrix(T), J: RealTensor4(T), C: RealMatrix(T)) !void {
    if (F.rows != F_MS.rows) {oneAO2MS(T, F_MS, F, C);} else {oneAO2MO(T, F_MS, F, C);}
    if (F.rows != F_MS.rows) {twoAO2MS(T, J_MS, J, C);} else {twoAO2MO(T, J_MS, J, C);}
}
