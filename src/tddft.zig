//! Module to perform the Time-Dependent Density Functional Theory (TDDFT) calculations.

const std = @import("std");

const classical_particle = @import("classical_particle.zig");
const device_write = @import("device_write.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const energy_derivative = @import("energy_derivative.zig");
const hartree_fock = @import("hartree_fock.zig");
const frequency_analysis = @import("frequency_analysis.zig");
const global_variables = @import("global_variables.zig");
const integral_transform = @import("integral_transform.zig");
const molecular_integrals = @import("molecular_integrals.zig");
const particle_optimization = @import("particle_optimization.zig");
const real_matrix = @import("real_matrix.zig");
const real_tensor_four = @import("real_tensor_four.zig");
const real_vector = @import("real_vector.zig");

const ClassicalParticle = classical_particle.ClassicalParticle;
const RealMatrix = real_matrix.RealMatrix;
const RealTensor4 = real_tensor_four.RealTensor4;
const RealVector = real_vector.RealVector;

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

const TEST_TOLERANCE = global_variables.TEST_TOLERANCE;

/// The TDDFT target options struct.
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

        hartree_fock: hartree_fock.Options(T),

        gradient: ?Gradient = null,
        hessian: ?Hessian = null,
        optimize: ?Optimize = null
    };
}

/// The TDDFT target output struct.
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

/// Run the TDDFT calculation with the given options for the provided system, returning the output struct with energy, gradient, hessian, and frequencies if requested.
pub fn run(comptime T: type, opt: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (opt.gradient != null and opt.gradient.? == .analytic) {

        std.log.err("ANALYTIC GRADIENT NOT IMPLEMENTED FOR TDDFT METHOD", .{});

        return error.InvalidInput;
    }

    if (opt.hessian != null and opt.hessian.? == .analytic) {

        std.log.err("ANALYTIC HESSIAN NOT IMPLEMENTED FOR TDDFT METHOD", .{});

        return error.InvalidInput;
    }

    var system = try classical_particle.read(T, opt.hartree_fock.system, opt.hartree_fock.charge, 0, allocator); defer system.deinit(allocator);

    if (enable_printing) {try print("\nINPUT GEOMETRY (Å):\n", .{}); try printClassicalParticleAsMolecule(T, system, null);}

    if (opt.optimize != null) {

        const optimized_system = try particleSteepestDescent(T, opt, system, tddft, "TDDFT", enable_printing, allocator);

        system.deinit(allocator); system = optimized_system;
    }

    if (enable_printing and opt.optimize != null) {try print("\nOPTIMIZED GEOMETRY (Å):\n", .{}); try printClassicalParticleAsMolecule(T, system, null);}

    var output = try tddft(T, opt, system, enable_printing, allocator); errdefer output.deinit(allocator);

    output.G = if (opt.gradient != null) try nuclearGradient(T, opt, system, tddft, "TDDFT", enable_printing, allocator) else null;

    if (output.G) |G| {try print("\nTDDFT NUCLEAR GRADIENT (Eh/Bohr):\n", .{}); try printRealMatrix(T, G);}

    output.H = if (opt.hessian != null) try nuclearHessian(T, opt, system, tddft, "TDDFT", enable_printing, allocator) else null;

    if (output.H) |H| output.frequencies = try particleHarmonicFrequencies(T, system, H, allocator);

    if (output.frequencies) |freqs| {try print("\nTDDFT VIBRATIONAL FREQUENCIES (CM^-1):\n", .{}); try printRealMatrix(T, freqs.asMatrix());}

    return output;
}

/// Primary function to run the tddft calculation, returning the energy and the output of the underlying hartree-fock calculation.
pub fn tddft(comptime T: type, opt: Options(T), system: ClassicalParticle(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (opt.hartree_fock.dft == null) {

        std.log.err("TDDFT CALCULATION REQUIRES A DFT FUNCTIONAL TO BE SPECIFIED IN THE HARTREE-FOCK OPTIONS", .{});

        return error.InvalidInput;
    }

    const hf_output = try scf(T, opt.hartree_fock, system, enable_printing, allocator);

    try print("\nTDDFT CALCULATION NOT FULLY IMPLEMENTED, RETURNING HARTREE-FOCK ENERGY\n", .{});

    const energy: T = 0;

    return .{
        .hf_output = hf_output, .energy = hf_output.energy + energy, .G = null, .H = null, .frequencies = null,
    };
}

/// Function to perform all integrals transformations used in the Moller-Plesset calculations.
pub fn transform(comptime T: type, F_MS: *RealMatrix(T), J_MS: *RealTensor4(T), F: RealMatrix(T), J: RealTensor4(T), C: RealMatrix(T)) !void {
    if (F.rows != F_MS.rows) {oneAO2MS(T, F_MS, F, C);} else {oneAO2MO(T, F_MS, F, C);}
    if (F.rows != F_MS.rows) {twoAO2MS(T, J_MS, J, C);} else {twoAO2MO(T, J_MS, J, C);}
}
