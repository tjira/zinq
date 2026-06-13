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
};

const Write = struct {
    kinetic: ?[]const u8 = null,
    overlap: ?[]const u8 = null,
    coulomb: ?[]const u8 = null,
    nuclear: ?[]const u8 = null,
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

        J: ?Tensor(T, 4) = null,

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            if (self.S) |*S| S.deinit(gpa);
            if (self.K) |*K| K.deinit(gpa);
            if (self.V) |*V| V.deinit(gpa);
            if (self.J) |*J| J.deinit(gpa);

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

    if (opt.calculate.kinetic) {
        ints.K = if (opt.spin) try ints.sys.kineticSpin(gpa) else try ints.sys.kinetic(gpa);

        if (log) try printf(io, "KINETIC INTEGRALS: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.calculate.nuclear) {
        ints.V = if (opt.spin) try ints.sys.nuclearSpin(gpa) else try ints.sys.nuclear(gpa);

        if (log) try printf(io, "NUCLEAR INTEGRALS: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.calculate.coulomb) {
        ints.J = if (opt.spin) try ints.sys.coulombSpin(gpa) else try ints.sys.coulomb(gpa);

        if (log) try printf(io, "COULOMB INTEGRALS: {f}\n", .{timer.untilNow(io, .real)});
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
        const S = ints.S orelse @panic("OVERLAP MATRIX NOT CALCULATED");

        try writeMatrix(T, io, fname, S);

        if (log) try printf(io, "OVERLAP INTEGRALS WRITING: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.write.kinetic) |fname| {
        const K = ints.K orelse @panic("KINETIC MATRIX NOT CALCULATED");

        try writeMatrix(T, io, fname, K);

        if (log) try printf(io, "KINETIC INTEGRALS WRITING: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.write.nuclear) |fname| {
        const V = ints.V orelse @panic("NUCLEAR MATRIX NOT CALCULATED");

        try writeMatrix(T, io, fname, V);

        if (log) try printf(io, "NUCLEAR INTEGRALS WRITING: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.write.coulomb) |fname| {
        const J = ints.J orelse @panic("COULOMB MATRIX NOT CALCULATED");

        try writeMatrix(T, io, fname, J.asMatrix());

        if (log) try printf(io, "COULOMB INTEGRALS WRITING: {f}\n", .{timer.untilNow(io, .real)});
    }

    return ints;
}
