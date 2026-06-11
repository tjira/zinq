const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("libint.zig").MolecularSystem;

const printf = @import("read_write.zig").printf;
const writeMatrix = @import("read_write.zig").writeMatrix;

// OPTIONS =============================================================================================================

const Write = struct {
    kinetic: ?[]const u8 = null,
    overlap: ?[]const u8 = null,
    coulomb: ?[]const u8 = null,
    nuclear: ?[]const u8 = null,
};

pub const Options = struct {
    system: [:0]const u8,
    basis: [:0]const u8,

    write: Write = .{},
};

// RUN =================================================================================================================

pub fn run(comptime T: type, io: std.Io, opt: Options, _: bool, gpa: Allocator, _: Allocator) !void {
    var timer = std.Io.Timestamp.now(io, .real);

    var sys = try MolecularSystem(T).init(opt.system, opt.basis);
    defer sys.deinit();

    if (opt.write.overlap) |fname| {
        var S = try sys.overlap(gpa);
        defer S.deinit(gpa);

        try writeMatrix(T, io, fname, S);
    }

    if (opt.write.kinetic) |fname| {
        var K = try sys.kinetic(gpa);
        defer K.deinit(gpa);

        try writeMatrix(T, io, fname, K);
    }

    if (opt.write.nuclear) |fname| {
        var V = try sys.nuclear(gpa);
        defer V.deinit(gpa);

        try writeMatrix(T, io, fname, V);
    }

    if (opt.write.coulomb) |fname| {
        var J = try sys.coulomb(gpa);
        defer J.deinit(gpa);

        try writeMatrix(T, io, fname, J);
    }

    try printf(io, "\n{f}\n", .{timer.untilNow(io, .real)});
}
