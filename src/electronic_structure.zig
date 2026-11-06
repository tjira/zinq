//! Orchestrates electronic structure calculations, geometry optimizations, and nuclear derivative evaluations.

const std = @import("std");

const cblas = @import("cimport.zig").cblas;

const Allocator = std.mem.Allocator;

const casscf = @import("casscf.zig");
const configuration_interaction = @import("configuration_interaction.zig");
const hartree_fock = @import("hartree_fock.zig");
const moller_plesset = @import("moller_plesset.zig");
const validation = @import("validation.zig");

const DerivativeContext = @import("nuclear_derivative.zig").DerivativeContext;
const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;

const bfgs = @import("molecular_optimization.zig").bfgs;
const calculateHarmonicFrequencies = @import("frequency_analysis.zig").calculateHarmonicFrequencies;
const calculateNumericalGradient = @import("nuclear_derivative.zig").calculateNumericalGradient;
const calculateNumericalHessian = @import("nuclear_derivative.zig").calculateNumericalHessian;
const calculateNumericalNacv = @import("nuclear_derivative.zig").calculateNumericalNacv;
const calculateThermochemistry = @import("frequency_analysis.zig").calculateThermochemistry;
const exportIfBuiltin = @import("molecular_integrals.zig").exportIfBuiltin;
const generateCasDets = @import("determinant.zig").generateCasDets;
const generateDets = @import("determinant.zig").generateDets;
const printGeometry = @import("read_write.zig").printGeometry;
const printHarmonicFrequencies = @import("frequency_analysis.zig").printHarmonicFrequencies;
const printThermochemistry = @import("frequency_analysis.zig").printThermochemistry;
const printf = @import("read_write.zig").printf;
const steepestDescent = @import("molecular_optimization.zig").steepestDescent;
const writeMatrix = @import("read_write.zig").writeMatrix;
const writeXyzFile = @import("read_write.zig").writeXyzFile;

/// Categories of electronic structure calculation methods supported by the computation engine.
pub const MethodKind = enum {
    casscf,
    configuration_interaction,
    hartree_fock,
    moller_plesset,
};

/// Parameters defining the nuclear geometry, basis set, and spin state of a molecular system.
pub const SystemParams = struct {
    basis: []const u8,
    charge: i32,
    multiplicity: u32,
    system: []const u8,
};

/// Generates electronic Slater determinants for complete active space or truncated CI excitation spaces.
pub fn generateDeterminants(nel: usize, nbf: usize, opt: anytype, gpa: Allocator) !std.ArrayList([]const usize) {
    if (comptime @TypeOf(opt) == casscf.Options) {
        const ncore = (nel - opt.active_electrons) / 2;

        return generateCasDets(ncore, opt.active_orbitals, opt.active_electrons, opt.multiplicity, gpa);
    }

    if (comptime @TypeOf(opt) == configuration_interaction.Options) {
        return switch (opt.space) {
            .cas => |c| blk: {
                const ncore = (nel - c.active_electrons) / 2;

                break :blk generateCasDets(ncore, c.active_orbitals, c.active_electrons, c.multiplicity, gpa);
            },
            .excitations => |ex| generateDets(nel, 2 * nbf, ex, gpa),
        };
    }

    @compileError("DETERMINANT GENERATION NOT IMPLEMENTED FOR OPTIONS TYPE: " ++ @typeName(@TypeOf(opt)));
}

/// Extracts the one-particle reduced density matrix from an electronic structure result.
pub fn getDensity(comptime T: type, result: anytype) Matrix(T) {
    const ResType = switch (@typeInfo(@TypeOf(result))) {
        .pointer => |p| p.child,
        else => @TypeOf(result),
    };

    if (comptime @hasField(ResType, "P_ao")) {
        return result.P_ao;
    }

    return result.hartree_fock.P_ao;
}

/// Identifies the active electronic structure method type from configuration options.
pub fn getMethodKind(comptime OptType: type) MethodKind {
    return switch (OptType) {
        casscf.Options => .casscf,
        configuration_interaction.Options => .configuration_interaction,
        hartree_fock.Options => .hartree_fock,
        moller_plesset.Options => .moller_plesset,
        else => @compileError("UNSUPPORTED OPTIONS TYPE: " ++ @typeName(OptType)),
    };
}

/// Derives state-pair output paths for nonadiabatic coupling vectors by formatting filename indices.
pub fn getNacvPath(path: []const u8, p0: u32, p1: u32, gpa: Allocator) ![]const u8 {
    const ext = std.fs.path.extension(path);

    return try std.fmt.allocPrint(gpa, "{s}_{d}_{d}{s}", .{ path[0 .. path.len - ext.len], p0, p1, ext });
}

/// Extracts molecular orbital coefficient matrices from converged electronic structure results.
pub fn getOrbitals(comptime T: type, result: anytype) Matrix(T) {
    const ResType = switch (@typeInfo(@TypeOf(result))) {
        .pointer => |p| p.child,
        else => @TypeOf(result),
    };

    if (comptime @hasField(ResType, "C_mo")) {
        return result.C_mo;
    }

    return result.hartree_fock.C_mo;
}

/// Derives output paths for state nuclear derivatives by appending the state index before extension.
pub fn getStatePath(path: []const u8, state: u32, gpa: Allocator) ![]const u8 {
    const ext = std.fs.path.extension(path);

    return try std.fmt.allocPrint(gpa, "{s}_{d}{s}", .{ path[0 .. path.len - ext.len], state, ext });
}

/// Extracts the molecular system geometry and basis specification from calculation results.
pub fn getSystem(comptime T: type, result: anytype) MolecularSystem(T) {
    const ResType = switch (@typeInfo(@TypeOf(result))) {
        .pointer => |p| p.child,
        else => @TypeOf(result),
    };

    if (comptime @hasField(ResType, "hartree_fock")) {
        return result.hartree_fock.ints.sys;
    }

    return result.ints.sys;
}

/// Extracts the nuclear geometry, basis set, and electronic spin parameters from method options.
pub fn getSystemParams(opt: anytype) SystemParams {
    if (comptime isPostScf(@TypeOf(opt))) {
        return .{
            .basis = opt.hartree_fock.basis,
            .charge = opt.hartree_fock.charge,
            .multiplicity = opt.hartree_fock.multiplicity,
            .system = opt.hartree_fock.system,
        };
    }

    if (!comptime isPostScf(@TypeOf(opt))) {
        return .{
            .basis = opt.basis,
            .charge = opt.charge,
            .multiplicity = opt.multiplicity,
            .system = opt.system,
        };
    }
}

/// Determines whether spin-unrestricted or generalized orbital treatment is enabled.
pub fn isGeneralized(opt: anytype) bool {
    return if (comptime isPostScf(@TypeOf(opt))) opt.hartree_fock.generalized else opt.generalized;
}

/// Identifies whether the electronic structure method resolves multiple electronic states.
pub fn isMultiState(comptime Opt: type) bool {
    return switch (getMethodKind(Opt)) {
        .casscf, .configuration_interaction => true,
        .hartree_fock, .moller_plesset => false,
    };
}

/// Indicates whether the specified electronic structure method builds upon an SCF reference.
pub fn isPostScf(comptime Opt: type) bool {
    return switch (Opt) {
        casscf.Options, configuration_interaction.Options, moller_plesset.Options => true,
        hartree_fock.Options => false,
        else => @compileError("UNSUPPORTED OPTIONS TYPE: " ++ @typeName(Opt)),
    };
}

/// Runs an electronic structure method on a molecular system specified by file paths and basis set.
pub fn run(comptime T: type, io: std.Io, opt: anytype, log: bool, gpa: Allocator) !getMethod(@TypeOf(opt)).Result(T) {
    try validation.checkInvalidInput(opt);

    cblas.openblas_set_num_threads(@intCast(opt.nthreads));

    const params = getSystemParams(opt);

    const basis_path = try exportIfBuiltin(io, params.basis, gpa);

    defer if (std.mem.startsWith(u8, params.basis, "builtin:")) {
        std.Io.Dir.cwd().deleteFile(io, basis_path) catch {};
    };

    var sys = try MolecularSystem(T).init(io, params.system, basis_path, params.charge, params.multiplicity, gpa);
    defer sys.deinit(gpa);

    return try runFromSystem(T, io, opt, &sys, null, null, log, gpa);
}

/// Executes an electronic structure calculation and nuclear derivatives on an initialized system.
pub fn runFromSystem(comptime T: type, io: std.Io, opt: anytype, sys: *MolecularSystem(T), Pg: ?Matrix(T), C_mo: ?Matrix(T), log: bool, gpa: Allocator) !getMethod(@TypeOf(opt)).Result(T) {
    const Method = getMethod(@TypeOf(opt));

    try validation.checkInvalidInput(opt);

    cblas.openblas_set_num_threads(@intCast(opt.nthreads));

    if (log and comptime !isPostScf(@TypeOf(opt))) {
        try printGeometry(T, io, sys.atoms, sys.coors, "MOLECULAR GEOMETRY (Å)");
    }

    var final_Pg: ?Matrix(T) = null;
    defer if (final_Pg) |*p| p.deinit(gpa);

    if (opt.optimize) |o| {
        switch (o) {
            .steepest_descent => final_Pg = try steepestDescent(T, io, runFromSystem, opt, sys, Pg, log, gpa),
            .bfgs => final_Pg = try bfgs(T, io, runFromSystem, opt, sys, Pg, log, gpa),
        }

        if (log) {
            try printGeometry(T, io, sys.atoms, sys.coors, "OPTIMIZED MOLECULAR GEOMETRY (Å)");
        }
    }

    var result = if (comptime isPostScf(@TypeOf(opt))) blk: {
        const ref = try computeHfReference(T, io, opt, sys, final_Pg orelse Pg, log, gpa);

        if (comptime getMethodKind(@TypeOf(opt)) == .casscf) {
            break :blk try Method.solve(T, io, opt, sys, ref, C_mo, log, gpa);
        }

        break :blk try Method.solve(T, io, opt, sys, ref, log, gpa);
    } else try hartree_fock.solve(T, io, opt, sys, final_Pg orelse Pg, log, gpa);

    errdefer result.deinit(gpa);

    try evaluateGradients(T, io, opt, sys, &result, log, gpa);

    try evaluateHessians(T, io, opt, sys, &result, log, gpa);

    if (@hasField(@TypeOf(opt), "nacv") and @hasField(@TypeOf(result), "nacv")) {
        try evaluateNacvs(T, io, opt, sys, &result, log, gpa);
    }

    try exportOutputs(T, io, opt, sys, result, gpa);

    return result;
}

/// Evaluates baseline Hartree-Fock molecular orbitals and density matrices for correlated methods.
fn computeHfReference(comptime T: type, io: std.Io, opt: anytype, sys: *MolecularSystem(T), Pg: ?Matrix(T), log: bool, gpa: Allocator) !hartree_fock.Result(T) {
    var hf_opt = opt.hartree_fock;

    if (opt.gradient) |g| if (g == .analytic) {
        hf_opt.gradient = .{ .analytic = .{} };

        if (comptime getMethodKind(@TypeOf(opt)) != .casscf) {
            if (hf_opt.response == null) {
                hf_opt.response = .{};
            }
        }
    };

    const hfres = try runFromSystem(T, io, hf_opt, sys, Pg, null, log, gpa);

    cblas.openblas_set_num_threads(@intCast(opt.nthreads));

    return hfres;
}

/// Evaluates analytic or finite-difference nuclear gradients and prints coordinate components.
fn evaluateGradients(comptime T: type, io: std.Io, opt: anytype, sys: *MolecularSystem(T), result: anytype, log: bool, gpa: Allocator) !void {
    if (opt.gradient == null) return;

    const grad_states = getGradientStates(opt);

    switch (opt.gradient.?) {
        .analytic => {
            if (log) {
                try logGradients(T, io, opt, result.grad, grad_states, true);
            }
        },
        .numeric => {
            const Pg, const C_mo = .{ getDensity(T, result), getOrbitals(T, result) };

            const ctx: DerivativeContext(T) = .{ .gpa = gpa, .io = io, .log = log, .sys = sys, .Pg = Pg, .C_mo = C_mo };

            const num_grads = try calculateNumericalGradient(T, ctx, runFromSystem, opt, grad_states);

            errdefer {
                for (0..num_grads.len) |k| num_grads[k].deinit(gpa);

                gpa.free(num_grads);
            }

            if (log) {
                try logGradients(T, io, opt, num_grads, grad_states, false);
            }

            for (0..result.grad.len) |k| result.grad[k].deinit(gpa);

            if (result.grad.len > 0) gpa.free(result.grad);

            result.grad = num_grads;
        },
    }
}

/// Computes numerical nuclear Hessians, harmonic vibrational frequencies, and thermochemistry.
fn evaluateHessians(comptime T: type, io: std.Io, opt: anytype, sys: *MolecularSystem(T), result: anytype, log: bool, gpa: Allocator) !void {
    if (opt.hessian == null) return;

    const hess_states = getHessianStates(opt);

    const Pg, const C_mo = .{ getDensity(T, result), getOrbitals(T, result) };

    const ctx: DerivativeContext(T) = .{ .gpa = gpa, .io = io, .log = log, .sys = sys, .Pg = Pg, .C_mo = C_mo };

    const hess = try calculateNumericalHessian(T, ctx, runFromSystem, opt, hess_states);

    errdefer {
        for (0..hess.len) |k| hess[k].deinit(gpa);

        gpa.free(hess);
    }

    if (log and opt.frequency != null) {
        const mult = getSystemParams(opt).multiplicity;

        for (hess_states, 0..) |state, idx| {
            var freqs = try calculateHarmonicFrequencies(T, hess[idx], sys.*, gpa);
            defer freqs.deinit(gpa);

            const method_str = try getHessianMethodLabel(opt, state, gpa);
            defer gpa.free(method_str);

            try printHarmonicFrequencies(T, io, freqs, method_str);

            const energy_val = if (state < result.energy.len) result.energy[state] else result.energy[0];

            const thermo = try calculateThermochemistry(T, opt.frequency.?, sys.*, freqs, mult, gpa);

            try printThermochemistry(T, io, opt.frequency.?, thermo, energy_val, method_str);
        }
    }

    for (0..result.hess.len) |k| result.hess[k].deinit(gpa);

    if (result.hess.len > 0) gpa.free(result.hess);

    result.hess = hess;
}

/// Evaluates nonadiabatic coupling vectors between electronic states using numerical finite differences.
fn evaluateNacvs(comptime T: type, io: std.Io, opt: anytype, sys: *MolecularSystem(T), result: anytype, log: bool, gpa: Allocator) !void {
    if (opt.nacv == null) return;

    const nacv_pairs = getNacvPairs(opt);

    const Pg, const C_mo = .{ getDensity(T, result), getOrbitals(T, result) };

    const ctx: DerivativeContext(T) = .{ .gpa = gpa, .io = io, .log = log, .sys = sys, .Pg = Pg, .C_mo = C_mo };

    const nacv = try calculateNumericalNacv(T, ctx, runFromSystem, opt, nacv_pairs);

    errdefer {
        for (0..nacv.len) |k| nacv[k].deinit(gpa);

        gpa.free(nacv);
    }

    if (log) {
        const method_label = if (comptime getMethodKind(@TypeOf(opt)) == .casscf) "CASSCF" else "CI";

        for (nacv_pairs, 0..) |pair, idx| {
            const fmt = "\n{s} NONADIABATIC COUPLING VECTOR BETWEEN STATES {d} AND {d} (1/a0)\n";

            try printf(io, fmt, .{ method_label, pair[0], pair[1] });

            for (0..nacv[idx].nrow()) |j| for (0..nacv[idx].ncol()) |k| {
                try printf(io, "{d:20.14}{s}", .{ nacv[idx].at(j, k), if (k == 2) "\n" else " " });
            };
        }
    }

    for (0..result.nacv.len) |k| result.nacv[k].deinit(gpa);

    if (result.nacv.len > 0) gpa.free(result.nacv);

    result.nacv = nacv;
}

/// Exports computed molecular geometries, energy gradients, Hessians, and nonadiabatic couplings to disk.
fn exportOutputs(comptime T: type, io: std.Io, opt: anytype, sys: *MolecularSystem(T), result: anytype, gpa: Allocator) !void {
    if (opt.write.geometry) |fname| {
        try writeXyzFile(T, io, fname, sys.atoms, sys.coors);
    }

    if (opt.write.gradient) |fname| {
        if (result.grad.len == 1) {
            try writeMatrix(T, io, fname, result.grad[0]);
        }

        if (result.grad.len > 1) {
            const grad_states = getGradientStates(opt);

            for (result.grad, 0..) |g, idx| {
                const path = try getStatePath(fname, grad_states[idx], gpa);
                defer gpa.free(path);

                try writeMatrix(T, io, path, g);
            }
        }
    }

    if (opt.write.hessian) |fname| {
        if (result.hess.len == 1) {
            try writeMatrix(T, io, fname, result.hess[0]);
        }

        if (result.hess.len > 1) {
            const hess_states = getHessianStates(opt);

            for (result.hess, 0..) |h, idx| {
                const path = try getStatePath(fname, hess_states[idx], gpa);
                defer gpa.free(path);

                try writeMatrix(T, io, path, h);
            }
        }
    }

    if (@hasField(@TypeOf(opt.write), "nacv") and @hasField(@TypeOf(result), "nacv")) {
        if (opt.write.nacv) |fname| {
            const nacv_pairs = getNacvPairs(opt);

            if (result.nacv.len == 1) {
                try writeMatrix(T, io, fname, result.nacv[0]);
            }

            if (result.nacv.len > 1) {
                for (result.nacv, 0..) |n, idx| {
                    const path = try getNacvPath(fname, nacv_pairs[idx][0], nacv_pairs[idx][1], gpa);
                    defer gpa.free(path);

                    try writeMatrix(T, io, path, n);
                }
            }
        }
    }

    if (@hasField(@TypeOf(opt.write), "orbitals") and @hasField(@TypeOf(result), "C_mo")) {
        if (opt.write.orbitals) |fname| {
            try writeMatrix(T, io, fname, result.C_mo);
        }
    }
}

/// Extracts the electronic state indices targeted for nuclear energy gradient evaluation.
fn getGradientStates(opt: anytype) []const u32 {
    if (!@hasField(@TypeOf(opt), "gradient") or opt.gradient == null) return &.{};

    return switch (opt.gradient.?) {
        inline else => |g| if (@hasField(@TypeOf(g), "states")) g.states orelse &.{0} else &.{0},
    };
}

/// Generates descriptive method labels for printing harmonic vibrational frequency tables.
fn getHessianMethodLabel(opt: anytype, state: u32, gpa: Allocator) ![]const u8 {
    return switch (comptime getMethodKind(@TypeOf(opt))) {
        .casscf => try std.fmt.allocPrint(gpa, "CASSCF STATE {d}", .{state}),
        .configuration_interaction => try std.fmt.allocPrint(gpa, "CI STATE {d}", .{state}),
        .hartree_fock => try gpa.dupe(u8, if (opt.dft != null) "DFT" else "HARTREE-FOCK"),
        .moller_plesset => try std.fmt.allocPrint(gpa, "MP{d}", .{opt.order}),
    };
}

/// Extracts the electronic state indices targeted for nuclear Hessian evaluation.
fn getHessianStates(opt: anytype) []const u32 {
    if (!@hasField(@TypeOf(opt), "hessian") or opt.hessian == null) return &.{};

    return switch (opt.hessian.?) {
        inline else => |h| if (@hasField(@TypeOf(h), "states")) h.states orelse &.{0} else &.{0},
    };
}

/// Determines the electronic structure method module corresponding to an options configuration type.
fn getMethod(comptime OptType: type) type {
    return switch (OptType) {
        casscf.Options => casscf,
        configuration_interaction.Options => configuration_interaction,
        hartree_fock.Options => hartree_fock,
        moller_plesset.Options => moller_plesset,
        else => @compileError("UNSUPPORTED OPTIONS TYPE: " ++ @typeName(OptType)),
    };
}

/// Extracts electronic state index pairs targeted for nonadiabatic coupling evaluation.
fn getNacvPairs(opt: anytype) []const [2]u32 {
    if (!@hasField(@TypeOf(opt), "nacv") or opt.nacv == null) return &.{};

    return switch (opt.nacv.?) {
        inline else => |n| n.pairs orelse &[_][2]u32{.{ 0, 1 }},
    };
}

/// Formats and prints nuclear energy gradient matrices to standard output.
fn logGradients(comptime T: type, io: std.Io, opt: anytype, grad: []Matrix(T), grad_states: []const u32, is_analytic: bool) !void {
    const grad_type_str = if (is_analytic) "ANALYTIC" else "NUMERIC";

    for (0..grad.len) |i| {
        switch (comptime getMethodKind(@TypeOf(opt))) {
            .casscf => {
                const s = if (i < grad_states.len) grad_states[i] else @as(u32, @intCast(i));

                try printf(io, "\nCASSCF STATE {d} {s} NUCLEAR ENERGY GRADIENT (Eh/a0)\n", .{ s, grad_type_str });
            },
            .configuration_interaction => {
                const s = if (i < grad_states.len) grad_states[i] else @as(u32, @intCast(i));

                try printf(io, "\nCI STATE {d} {s} NUCLEAR ENERGY GRADIENT (Eh/a0)\n", .{ s, grad_type_str });
            },
            .hartree_fock => {
                const method_str = if (opt.dft != null) "DFT" else "HARTREE-FOCK";

                const hf_type_str = if (is_analytic) "ANALYTICAL" else "NUMERICAL";

                try printf(io, "\n{s} {s} NUCLEAR ENERGY GRADIENT (Eh/a0)\n", .{ method_str, hf_type_str });
            },
            .moller_plesset => {
                const mp_type_str = if (is_analytic) "ANALYTICAL" else "NUMERIC";

                try printf(io, "\nMP{d:02} {s} NUCLEAR ENERGY GRADIENT (Eh/a0)\n", .{ opt.order, mp_type_str });
            },
        }

        for (0..grad[i].nrow()) |j| for (0..grad[i].ncol()) |k| {
            try printf(io, "{d:20.14}{s}", .{ grad[i].at(j, k), if (k == 2) "\n" else " " });
        };
    }
}
