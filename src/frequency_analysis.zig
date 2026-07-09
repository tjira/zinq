const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const Vector = @import("tensor.zig").Vector;

const eighSlice = @import("linear_algebra.zig").eighSlice;
const getMass = @import("constant.zig").getMass;
const printf = @import("read_write.zig").printf;

const AN2M = @import("constant.zig").AN2M;
const AN2SM = @import("constant.zig").AN2SM;
const AU2CM = @import("constant.zig").AU2CM;

pub fn calculateHarmonicFrequencies(comptime T: type, hessian: Matrix(T), atoms: []const i32, gpa: Allocator) !Vector(T) {
    std.debug.assert(hessian.nrow() == atoms.len * 3);

    var HM = try Matrix(T).init(hessian.nrow(), hessian.nrow(), gpa);
    defer HM.deinit(gpa);

    for (0..atoms.len) |i| for (0..3) |xyz_i| for (0..atoms.len) |j| for (0..3) |xyz_j| {
        const row = i * 3 + xyz_i;
        const col = j * 3 + xyz_j;

        const mass_i = try getMass(T, atoms[i]);
        const mass_j = try getMass(T, atoms[j]);

        HM.ptr(row, col).* = hessian.at(row, col) / std.math.sqrt(mass_i * mass_j);
    };

    var w = try Vector(T).init(hessian.nrow(), gpa);
    errdefer w.deinit(gpa);

    var u = try Matrix(T).init(hessian.nrow(), hessian.nrow(), gpa);
    defer u.deinit(gpa);

    try eighSlice(T, w.data, u.data, HM.data);

    for (0..hessian.nrow()) |i| {
        w.ptr(i).* = std.math.sign(w.at(i)) * @as(T, @floatCast(std.math.sqrt(@abs(w.at(i)))));
    }

    return w;
}

pub fn printHarmonicFrequencies(comptime T: type, io: std.Io, freqs: Vector(T), method_str: []const u8) !void {
    try printf(io, "\n{s} HARMONIC VIBRATIONAL FREQUENCIES (cm^-1)\n", .{method_str});

    for (0..freqs.length()) |i| {
        try printf(io, "MODE {d:3}: {d:12.4}\n", .{ i + 1, AU2CM * freqs.at(i) });
    }
}
