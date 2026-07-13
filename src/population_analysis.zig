const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Vector = @import("tensor.zig").Vector;

const dot = @import("linear_algebra.zig").dot;
const getSymbol = @import("constant.zig").getSymbol;
const printf = @import("read_write.zig").printf;

pub fn mulliken(comptime T: type, sys: MolecularSystem(T), P: Matrix(T), S: Matrix(T), gpa: Allocator) !Vector(T) {
    var net_populations = try gpa.alloc(T, sys.atoms.len);
    defer gpa.free(net_populations);

    @memset(net_populations, 0);

    for (0..P.shape[0]) |u| {
        net_populations[@intCast(sys.bf2at[u % sys.nbf])] += dot(T, P.row(u), S.row(u));
    }

    var charges = try Vector(T).init(sys.atoms.len, gpa);
    errdefer charges.deinit(gpa);

    for (0..sys.atoms.len) |i| {
        charges.data[i] = @as(T, @floatFromInt(sys.atoms[i])) - net_populations[i];
    }

    return charges;
}

pub fn printMullikenCharges(comptime T: type, io: std.Io, sys: MolecularSystem(T), charges: Vector(T), method_str: []const u8) !void {
    try printf(io, "\n{s} MULLIKEN POPULATION ANALYSIS\n", .{method_str});

    for (0..sys.atoms.len) |i| {
        const sym = try getSymbol(sys.atoms[i]);

        try printf(io, "{s:4} {d:20.14}\n", .{ sym, charges.data[i] });
    }
}
