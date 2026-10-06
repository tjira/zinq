//! Calculates and exports molecular integrals including overlap, kinetic energy, nuclear attraction, and Coulomb repulsion.

const embed = @import("embed");
const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Tensor = @import("tensor.zig").Tensor;

const printf = @import("read_write.zig").printf;
const writeMatrix = @import("read_write.zig").writeMatrix;

/// Configuration parameters specifying molecular geometry, basis sets, spin properties, and task selectors for integral calculation.
pub const Options = struct {
    basis: []const u8,
    calculate: Calculate = .{},
    charge: i32 = 0,
    multiplicity: u32 = 1,
    nthreads: u32 = 1,
    spin: bool = false,
    system: []const u8,
    write: Write = .{},
};

/// Specifier of boolean flags to select which one-electron and two-electron molecular integrals and derivatives to compute.
const Calculate = struct {
    coulomb: bool = true,
    coulomb_d1: bool = false,
    hamiltonian: bool = true,
    hamiltonian_d1: bool = false,
    kinetic: bool = true,
    kinetic_d1: bool = false,
    nuclear: bool = true,
    nuclear_d1: bool = false,
    overlap: bool = true,
    overlap_d1: bool = false,
};

/// Options specifying the output file paths for writing calculated molecular integrals and their derivatives.
pub const Write = struct {
    coulomb: ?[]const u8 = null,
    coulomb_d1: ?[]const u8 = null,
    hamiltonian: ?[]const u8 = null,
    hamiltonian_d1: ?[]const u8 = null,
    kinetic: ?[]const u8 = null,
    kinetic_d1: ?[]const u8 = null,
    nuclear: ?[]const u8 = null,
    nuclear_d1: ?[]const u8 = null,
    overlap: ?[]const u8 = null,
    overlap_d1: ?[]const u8 = null,
};

/// Returns a generic type representing the calculated molecular integrals and their first-order nuclear derivatives.
pub fn Result(comptime T: type) type {
    return struct {
        sys: MolecularSystem(T),

        S: ?Matrix(T) = null,
        K: ?Matrix(T) = null,
        V: ?Matrix(T) = null,
        H: ?Matrix(T) = null,

        g: ?Tensor(T, 4) = null,

        dS: ?Tensor(T, 3) = null,
        dK: ?Tensor(T, 3) = null,
        dV: ?Tensor(T, 3) = null,
        dH: ?Tensor(T, 3) = null,

        dg: ?Tensor(T, 5) = null,

        /// Deallocates the computed one-electron and two-electron molecular integral tensors stored in the Result struct.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.sys.deinit(gpa);

            if (self.S) |*S| S.deinit(gpa);
            if (self.K) |*K| K.deinit(gpa);
            if (self.V) |*V| V.deinit(gpa);
            if (self.g) |*g| g.deinit(gpa);
            if (self.H) |*H| H.deinit(gpa);

            if (self.dS) |*dS| dS.deinit(gpa);
            if (self.dK) |*dK| dK.deinit(gpa);
            if (self.dV) |*dV| dV.deinit(gpa);
            if (self.dH) |*dH| dH.deinit(gpa);
            if (self.dg) |*dg| dg.deinit(gpa);
        }
    };
}

/// Automatically extracts and writes a builtin basis set to a temporary file if specified.
pub fn exportIfBuiltin(io: std.Io, basis: []const u8, gpa: Allocator) ![]const u8 {
    if (std.mem.startsWith(u8, basis, "builtin:")) {
        var result = try gpa.alloc(u8, basis["builtin:".len..].len);
        defer gpa.free(result);

        for (basis["builtin:".len..], 0..) |char, i| {
            if (char == '+') result[i] = 'p';
            if (char == '*') result[i] = 's';

            if (char != '+' and char != '*') {
                result[i] = std.ascii.toLower(char);
            }
        }

        const content = embed.bases.get(result) orelse {
            std.log.err("BUILTIN BASIS SET '{s}' NOT FOUND", .{result});

            return error.InvalidInput;
        };

        var file = try std.Io.Dir.cwd().createFile(io, "basis.g94", .{});
        defer file.close(io);

        try file.writeStreamingAll(io, content);

        return "basis.g94";
    }

    return basis;
}

/// Computes molecular integrals and nuclear derivatives from a geometry file and basis set path.
pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator) !Result(T) {
    try checkInvalidInput(opt);

    const basis_path = try exportIfBuiltin(io, opt.basis, gpa);

    defer if (std.mem.startsWith(u8, opt.basis, "builtin:")) {
        std.Io.Dir.cwd().deleteFile(io, basis_path) catch {};
    };

    var timer = std.Io.Timestamp.now(io, .real);

    var sys = try MolecularSystem(T).init(io, opt.system, basis_path, opt.charge, opt.multiplicity, gpa);
    defer sys.deinit(gpa);

    if (std.mem.startsWith(u8, opt.basis, "builtin:")) {
        try std.Io.Dir.cwd().deleteFile(io, basis_path);
    }

    if (log) try printf(io, "\nSYSTEM INITIALIZATION: {f}\n", .{timer.untilNow(io, .real)});

    return try runFromSystem(T, io, opt, sys, log, gpa);
}

/// Evaluates molecular integrals and gradients directly on an initialized molecular system.
pub fn runFromSystem(comptime T: type, io: std.Io, opt: Options, sys: MolecularSystem(T), log: bool, gpa: Allocator) !Result(T) {
    try checkInvalidInput(opt);

    var ints: Result(T) = .{ .sys = try sys.clone(gpa) };
    errdefer ints.deinit(gpa);

    const any_calc = blk: {
        inline for (@typeInfo(@TypeOf(opt.calculate)).@"struct".field_names) |f| {
            if (@field(opt.calculate, f)) break :blk true;
        }

        break :blk false;
    };

    if (log and any_calc) try std.Io.File.stdout().writeStreamingAll(io, "\n");

    var timer = std.Io.Timestamp.now(io, .real);

    if (opt.calculate.overlap) {
        ints.S = if (opt.spin) try sys.overlapSpin(opt.nthreads, gpa) else try sys.overlap(opt.nthreads, gpa);

        if (log) try printf(io, "OVERLAP INTEGRALS: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.calculate.kinetic or opt.calculate.hamiltonian) {
        ints.K = if (opt.spin) try sys.kineticSpin(opt.nthreads, gpa) else try sys.kinetic(opt.nthreads, gpa);

        if (log) try printf(io, "KINETIC INTEGRALS: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.calculate.nuclear or opt.calculate.hamiltonian) {
        ints.V = if (opt.spin) try sys.nuclearSpin(opt.nthreads, gpa) else try sys.nuclear(opt.nthreads, gpa);

        if (log) try printf(io, "NUCLEAR INTEGRALS: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.calculate.coulomb) {
        ints.g = if (opt.spin) try sys.coulombSpin(opt.nthreads, gpa) else try sys.coulomb(opt.nthreads, gpa);

        if (log) try printf(io, "COULOMB INTEGRALS: {f}\n", .{timer.untilNow(io, .real)});
    }

    if (opt.calculate.hamiltonian) {
        ints.H = try Matrix(T).init(ints.K.?.nrow(), ints.K.?.ncol(), gpa);

        for (0..ints.H.?.nrow()) |i| for (0..ints.H.?.ncol()) |j| {
            ints.H.?.ptr(i, j).* = ints.K.?.at(i, j) + ints.V.?.at(i, j);
        };
    }

    const any_deriv_calc = blk: {
        inline for (@typeInfo(@TypeOf(opt.calculate)).@"struct".field_names) |f| {
            if (comptime std.mem.endsWith(u8, f, "_d1")) {
                if (@field(opt.calculate, f)) break :blk true;
            }
        }

        break :blk false;
    };

    if (log and any_deriv_calc) try std.Io.File.stdout().writeStreamingAll(io, "\n");

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.calculate.overlap_d1) {
        ints.dS = if (opt.spin) try sys.overlapD1Spin(opt.nthreads, gpa) else try sys.overlapD1(opt.nthreads, gpa);

        if (log) try printf(io, "OVERLAP INTEGRALS DERIVATIVE: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.calculate.kinetic_d1 or opt.calculate.hamiltonian_d1) {
        ints.dK = if (opt.spin) try sys.kineticD1Spin(opt.nthreads, gpa) else try sys.kineticD1(opt.nthreads, gpa);

        if (log) try printf(io, "KINETIC INTEGRALS DERIVATIVE: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.calculate.nuclear_d1 or opt.calculate.hamiltonian_d1) {
        ints.dV = if (opt.spin) try sys.nuclearD1Spin(opt.nthreads, gpa) else try sys.nuclearD1(opt.nthreads, gpa);

        if (log) try printf(io, "NUCLEAR INTEGRALS DERIVATIVE: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.calculate.coulomb_d1) {
        ints.dg = if (opt.spin) try sys.coulombD1Spin(opt.nthreads, gpa) else try sys.coulombD1(opt.nthreads, gpa);

        if (log) try printf(io, "COULOMB INTEGRALS DERIVATIVE: {f}\n", .{timer.untilNow(io, .real)});
    }

    if (opt.calculate.hamiltonian_d1) {
        ints.dH = try Tensor(T, 3).init(ints.dK.?.shape, gpa);

        for (0..ints.dH.?.shape[0]) |k| for (0..ints.dH.?.shape[1]) |i| for (0..ints.dH.?.shape[2]) |j| {
            ints.dH.?.ptr(.{ k, i, j }).* = ints.dK.?.at(.{ k, i, j }) + ints.dV.?.at(.{ k, i, j });
        };
    }

    try writeIntegralsToFiles(T, io, opt, ints, log);

    return ints;
}

/// Validates that required options such as geometry and basis file paths are not empty.
fn checkInvalidInput(opt: Options) !void {
    if (opt.system.len == 0) {
        std.log.err("MOLECULAR SYSTEM XYZ PATH IS EMPTY", .{});

        return error.InvalidInput;
    }

    if (opt.basis.len == 0) {
        std.log.err("BASIS SET G94 PATH IS EMPTY", .{});

        return error.InvalidInput;
    }

    if (opt.nthreads == 0) {
        std.log.err("THREAD COUNT MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.write.coulomb != null and !opt.calculate.coulomb) {
        std.log.err("COULOMB INTEGRALS WRITE REQUESTED BUT COULOMB INTEGRALS ARE NOT CALCULATED", .{});

        return error.InvalidInput;
    }

    if (opt.write.coulomb_d1 != null and !opt.calculate.coulomb_d1) {
        std.log.err("COULOMB DERIVATIVE WRITE REQUESTED BUT COULOMB DERIVATIVE IS NOT CALCULATED", .{});

        return error.InvalidInput;
    }

    if (opt.write.hamiltonian != null and !opt.calculate.hamiltonian) {
        std.log.err("HAMILTONIAN INTEGRALS WRITE REQUESTED BUT HAMILTONIAN INTEGRALS ARE NOT CALCULATED", .{});

        return error.InvalidInput;
    }

    if (opt.write.hamiltonian_d1 != null and !opt.calculate.hamiltonian_d1) {
        std.log.err("HAMILTONIAN DERIVATIVE WRITE REQUESTED BUT HAMILTONIAN DERIVATIVE IS NOT CALCULATED", .{});

        return error.InvalidInput;
    }

    if (opt.write.kinetic != null and !opt.calculate.kinetic) {
        std.log.err("KINETIC INTEGRALS WRITE REQUESTED BUT KINETIC INTEGRALS ARE NOT CALCULATED", .{});

        return error.InvalidInput;
    }

    if (opt.write.kinetic_d1 != null and !opt.calculate.kinetic_d1) {
        std.log.err("KINETIC DERIVATIVE WRITE REQUESTED BUT KINETIC DERIVATIVE IS NOT CALCULATED", .{});

        return error.InvalidInput;
    }

    if (opt.write.nuclear != null and !opt.calculate.nuclear) {
        std.log.err("NUCLEAR INTEGRALS WRITE REQUESTED BUT NUCLEAR INTEGRALS ARE NOT CALCULATED", .{});

        return error.InvalidInput;
    }

    if (opt.write.nuclear_d1 != null and !opt.calculate.nuclear_d1) {
        std.log.err("NUCLEAR DERIVATIVE WRITE REQUESTED BUT NUCLEAR DERIVATIVE IS NOT CALCULATED", .{});

        return error.InvalidInput;
    }

    if (opt.write.overlap != null and !opt.calculate.overlap) {
        std.log.err("OVERLAP INTEGRALS WRITE REQUESTED BUT OVERLAP INTEGRALS ARE NOT CALCULATED", .{});

        return error.InvalidInput;
    }

    if (opt.write.overlap_d1 != null and !opt.calculate.overlap_d1) {
        std.log.err("OVERLAP DERIVATIVE WRITE REQUESTED BUT OVERLAP DERIVATIVE IS NOT CALCULATED", .{});

        return error.InvalidInput;
    }
}

/// Exports calculated molecular integral and derivative tensors to target files specified in options.
fn writeIntegralsToFiles(comptime T: type, io: std.Io, opt: Options, ints: Result(T), log: bool) !void {
    const any_write = blk: {
        inline for (@typeInfo(@TypeOf(opt.write)).@"struct".field_names) |f| {
            if (@field(opt.write, f) != null) break :blk true;
        }

        break :blk false;
    };

    if (log and any_write) try std.Io.File.stdout().writeStreamingAll(io, "\n");

    var timer = std.Io.Timestamp.now(io, .real);

    if (opt.write.overlap) |fname| {
        const S = ints.S orelse return error.OverlapMatrixNotCalculated;

        try writeMatrix(T, io, fname, S);

        if (log) try printf(io, "OVERLAP INTEGRALS WRITING: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.write.kinetic) |fname| {
        const K = ints.K orelse return error.KineticMatrixNotCalculated;

        try writeMatrix(T, io, fname, K);

        if (log) try printf(io, "KINETIC INTEGRALS WRITING: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.write.nuclear) |fname| {
        const V = ints.V orelse return error.NuclearMatrixNotCalculated;

        try writeMatrix(T, io, fname, V);

        if (log) try printf(io, "NUCLEAR INTEGRALS WRITING: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.write.coulomb) |fname| {
        const g = ints.g orelse return error.CoulombMatrixNotCalculated;

        try writeMatrix(T, io, fname, g.asMatrix());

        if (log) try printf(io, "COULOMB INTEGRALS WRITING: {f}\n", .{timer.untilNow(io, .real)});
    }

    if (opt.write.hamiltonian) |fname| {
        const H = ints.H orelse return error.HamiltonianNotCalculated;

        try writeMatrix(T, io, fname, H);

        if (log) try printf(io, "HAMILTONIAN INTEGRALS WRITING: {f}\n", .{timer.untilNow(io, .real)});
    }

    const any_deriv_write = blk: {
        inline for (@typeInfo(@TypeOf(opt.write)).@"struct".field_names) |f| {
            if (comptime std.mem.endsWith(u8, f, "_d1")) {
                if (@field(opt.write, f) != null) break :blk true;
            }
        }

        break :blk false;
    };

    if (log and any_deriv_write) try std.Io.File.stdout().writeStreamingAll(io, "\n");

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.write.overlap_d1) |fname| {
        const dS = ints.dS orelse return error.OverlapDerivativeMatrixNotCalculated;

        try writeMatrix(T, io, fname, dS.asMatrix());

        if (log) try printf(io, "OVERLAP DERIVATIVE INTEGRALS WRITING: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.write.kinetic_d1) |fname| {
        const dK = ints.dK orelse return error.KineticDerivativeMatrixNotCalculated;

        try writeMatrix(T, io, fname, dK.asMatrix());

        if (log) try printf(io, "KINETIC DERIVATIVE INTEGRALS WRITING: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.write.nuclear_d1) |fname| {
        const dV = ints.dV orelse return error.NuclearDerivativeMatrixNotCalculated;

        try writeMatrix(T, io, fname, dV.asMatrix());

        if (log) try printf(io, "NUCLEAR DERIVATIVE INTEGRALS WRITING: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.write.coulomb_d1) |fname| {
        const dg = ints.dg orelse return error.CoulombDerivativeMatrixNotCalculated;

        try writeMatrix(T, io, fname, dg.asMatrix());

        if (log) try printf(io, "COULOMB DERIVATIVE INTEGRALS WRITING: {f}\n", .{timer.untilNow(io, .real)});
    }

    if (opt.write.hamiltonian_d1) |fname| {
        const dH = ints.dH orelse return error.HamiltonianDerivativeNotCalculated;

        try writeMatrix(T, io, fname, dH.asMatrix());

        if (log) try printf(io, "HAMILTONIAN DERIVATIVE INTEGRALS WRITING: {f}\n", .{timer.untilNow(io, .real)});
    }
}
