//! Implements Mulliken and Löwdin population analysis to calculate partial atomic charges from the density and overlap matrices.

const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Vector = @import("tensor.zig").Vector;

const dot = @import("linear_algebra.zig").dot;
const eigh = @import("linear_algebra.zig").eigh;
const getSymbol = @import("constant.zig").getSymbol;
const mm = @import("linear_algebra.zig").mm;
const printf = @import("read_write.zig").printf;

/// Computes Löwdin net atomic charges using symmetric orthogonalization of the overlap matrix.
pub fn lowdin(comptime T: type, sys: MolecularSystem(T), P: Matrix(T), S: Matrix(T), gpa: Allocator) !Vector(T) {
    const n = S.nrow();

    var eigvals = try Vector(T).init(n, gpa);
    defer eigvals.deinit(gpa);

    var U = try Matrix(T).init(n, n, gpa);
    defer U.deinit(gpa);

    try eigh(T, &eigvals, &U, S);

    var U_scaled = try Matrix(T).init(n, n, gpa);
    defer U_scaled.deinit(gpa);

    for (0..n) |i| for (0..n) |j| {
        U_scaled.ptr(i, j).* = U.at(i, j) * @sqrt(@max(0, eigvals.at(j)));
    };

    var Shalf = try Matrix(T).init(n, n, gpa);
    defer Shalf.deinit(gpa);

    mm(T, &Shalf, U_scaled, U, 1.0, 0.0, false, true);

    var SP = try Matrix(T).init(n, n, gpa);
    defer SP.deinit(gpa);

    mm(T, &SP, Shalf, P, 1.0, 0.0, false, false);

    var net_populations = try gpa.alloc(T, sys.atoms.len);
    defer gpa.free(net_populations);

    @memset(net_populations, 0);

    for (0..n) |u| {
        net_populations[@intCast(sys.bf2at[u % sys.nbf])] += dot(T, SP.row(u), Shalf.row(u));
    }

    var charges = try Vector(T).init(sys.atoms.len, gpa);
    errdefer charges.deinit(gpa);

    for (0..sys.atoms.len) |i| {
        charges.data[i] = @as(T, @floatFromInt(sys.atoms[i])) - net_populations[i];
    }

    return charges;
}

/// Computes Mulliken net atomic charges by partitioning the electronic density matrix using the overlap matrix.
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

/// Formats and prints the calculated Löwdin atomic charges to the output.
pub fn printLowdinCharges(comptime T: type, io: std.Io, sys: MolecularSystem(T), charges: Vector(T), method_str: []const u8) !void {
    try printf(io, "\n{s} LÖWDIN POPULATION ANALYSIS\n", .{method_str});

    for (0..sys.atoms.len) |i| {
        const sym = try getSymbol(sys.atoms[i]);

        try printf(io, "{s:4} {d:20.14}\n", .{ sym, charges.data[i] });
    }
}

/// Formats and prints the calculated Mulliken atomic charges to the output.
pub fn printMullikenCharges(comptime T: type, io: std.Io, sys: MolecularSystem(T), charges: Vector(T), method_str: []const u8) !void {
    try printf(io, "\n{s} MULLIKEN POPULATION ANALYSIS\n", .{method_str});

    for (0..sys.atoms.len) |i| {
        const sym = try getSymbol(sys.atoms[i]);

        try printf(io, "{s:4} {d:20.14}\n", .{ sym, charges.data[i] });
    }
}
