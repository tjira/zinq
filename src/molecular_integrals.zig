const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Tensor = @import("tensor.zig").Tensor;

const printf = @import("read_write.zig").printf;
const writeMatrix = @import("read_write.zig").writeMatrix;

// OPTIONS =============================================================================================================

const Calculate = struct {
    kinetic: bool = true,
    overlap: bool = true,
    coulomb: bool = true,
    nuclear: bool = true,
    hmatrix: bool = true,

    kinetic_d1: bool = false,
    overlap_d1: bool = false,
    coulomb_d1: bool = false,
    nuclear_d1: bool = false,
    hmatrix_d1: bool = false,
};

const Write = struct {
    kinetic: ?[]const u8 = null,
    overlap: ?[]const u8 = null,
    coulomb: ?[]const u8 = null,
    nuclear: ?[]const u8 = null,
    hmatrix: ?[]const u8 = null,

    kinetic_d1: ?[]const u8 = null,
    overlap_d1: ?[]const u8 = null,
    coulomb_d1: ?[]const u8 = null,
    nuclear_d1: ?[]const u8 = null,
    hmatrix_d1: ?[]const u8 = null,
};

pub const Options = struct {
    system: []const u8,
    basis: []const u8,

    calculate: Calculate = .{},
    write: Write = .{},
    spin: bool = false,
};

// INTEGRALS ===========================================================================================================

pub fn Integrals(comptime T: type) type {
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

        pub fn deinit(self: *@This(), gpa: Allocator) void {
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

            self.sys.deinit(gpa);
        }
    };
}

// RUN =================================================================================================================

pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator) !Integrals(T) {
    var timer = std.Io.Timestamp.now(io, .real);

    var ints: Integrals(T) = .{ .sys = try MolecularSystem(T).init(opt.system, opt.basis, gpa) };
    errdefer ints.deinit(gpa);

    if (log) try printf(io, "\nSYSTEM INITIALIZATION: {f}\n", .{timer.untilNow(io, .real)});

    const any_calc = blk: {
        inline for (std.meta.fields(@TypeOf(opt.calculate))) |f| {
            if (@field(opt.calculate, f.name)) break :blk true;
        }

        break :blk false;
    };

    if (log and any_calc) try std.Io.File.stdout().writeStreamingAll(io, "\n");

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.calculate.overlap) {
        ints.S = if (opt.spin) try ints.sys.overlapSpin(gpa) else try ints.sys.overlap(gpa);

        if (log) try printf(io, "OVERLAP INTEGRALS: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.calculate.kinetic or opt.calculate.hmatrix) {
        ints.K = if (opt.spin) try ints.sys.kineticSpin(gpa) else try ints.sys.kinetic(gpa);

        if (log) try printf(io, "KINETIC INTEGRALS: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.calculate.nuclear or opt.calculate.hmatrix) {
        ints.V = if (opt.spin) try ints.sys.nuclearSpin(gpa) else try ints.sys.nuclear(gpa);

        if (log) try printf(io, "NUCLEAR INTEGRALS: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.calculate.coulomb) {
        ints.g = if (opt.spin) try ints.sys.coulombSpin(gpa) else try ints.sys.coulomb(gpa);

        if (log) try printf(io, "COULOMB INTEGRALS: {f}\n", .{timer.untilNow(io, .real)});
    }

    if (opt.calculate.hmatrix) {
        ints.H = try Matrix(T).init(ints.K.?.nrow(), ints.K.?.ncol(), gpa);

        for (0..ints.H.?.nrow()) |i| for (0..ints.H.?.ncol()) |j| {
            ints.H.?.ptr(i, j).* = ints.K.?.at(i, j) + ints.V.?.at(i, j);
        };
    }

    const any_deriv_calc = blk: {
        inline for (std.meta.fields(@TypeOf(opt.calculate))) |f| {
            if (comptime std.mem.endsWith(u8, f.name, "_d1")) {
                if (@field(opt.calculate, f.name)) break :blk true;
            }
        }

        break :blk false;
    };

    if (log and any_deriv_calc) try std.Io.File.stdout().writeStreamingAll(io, "\n");

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.calculate.overlap_d1) {
        ints.dS = if (opt.spin) try ints.sys.overlapD1Spin(gpa) else try ints.sys.overlapD1(gpa);

        if (log) try printf(io, "OVERLAP INTEGRALS DERIVATIVE: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.calculate.kinetic_d1 or opt.calculate.hmatrix_d1) {
        ints.dK = if (opt.spin) try ints.sys.kineticD1Spin(gpa) else try ints.sys.kineticD1(gpa);

        if (log) try printf(io, "KINETIC INTEGRALS DERIVATIVE: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.calculate.nuclear_d1 or opt.calculate.hmatrix_d1) {
        ints.dV = if (opt.spin) try ints.sys.nuclearD1Spin(gpa) else try ints.sys.nuclearD1(gpa);

        if (log) try printf(io, "NUCLEAR INTEGRALS DERIVATIVE: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.calculate.coulomb_d1) {
        ints.dg = if (opt.spin) try ints.sys.coulombD1Spin(gpa) else try ints.sys.coulombD1(gpa);

        if (log) try printf(io, "COULOMB INTEGRALS DERIVATIVE: {f}\n", .{timer.untilNow(io, .real)});
    }

    if (opt.calculate.hmatrix_d1) {
        ints.dH = try Tensor(T, 3).init(ints.dK.?.shape, gpa);

        for (0..ints.dH.?.shape[0]) |k| for (0..ints.dH.?.shape[1]) |i| for (0..ints.dH.?.shape[2]) |j| {
            ints.dH.?.ptr(.{ k, i, j }).* = ints.dK.?.at(.{ k, i, j }) + ints.dV.?.at(.{ k, i, j });
        };
    }

    const any_write = blk: {
        inline for (std.meta.fields(@TypeOf(opt.write))) |f| {
            if (@field(opt.write, f.name) != null) break :blk true;
        }

        break :blk false;
    };

    if (log and any_write) try std.Io.File.stdout().writeStreamingAll(io, "\n");

    timer = std.Io.Timestamp.now(io, .real);

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

    if (opt.write.hmatrix) |fname| {
        const H = ints.H orelse return error.hmatrixMatrixNotCalculated;

        try writeMatrix(T, io, fname, H);

        if (log) try printf(io, "HMATRIX INTEGRALS WRITING: {f}\n", .{timer.untilNow(io, .real)});
    }

    const any_deriv_write = blk: {
        inline for (std.meta.fields(@TypeOf(opt.write))) |f| {
            if (comptime std.mem.endsWith(u8, f.name, "_d1")) {
                if (@field(opt.write, f.name) != null) break :blk true;
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

    if (opt.write.hmatrix_d1) |fname| {
        const dH = ints.dH orelse return error.hmatrixDerivativeMatrixNotCalculated;

        try writeMatrix(T, io, fname, dH.asMatrix());

        if (log) try printf(io, "HMATRIX DERIVATIVE INTEGRALS WRITING: {f}\n", .{timer.untilNow(io, .real)});
    }

    return ints;
}
