const std = @import("std");

const Matrix = @import("tensor.zig").Matrix;

pub fn printf(io: std.Io, comptime format: []const u8, args: anytype) !void {
    var buffer: [4096]u8 = undefined;
    var writer = std.Io.File.stdout().writer(io, &buffer);

    try writer.interface.print(format, args);
    try writer.interface.flush();
}

pub fn writeMatrix(comptime T: type, io: std.Io, fname: []const u8, A: Matrix(T)) !void {
    if (@typeInfo(T) == .@"struct") {
        return try writeMatrixComplex(T, io, fname, A);
    }

    try writeMatrixReal(T, io, fname, A);
}

pub fn writeMatrixReal(comptime T: type, io: std.Io, fname: []const u8, A: Matrix(T)) !void {
    var file = try std.Io.Dir.cwd().createFile(io, fname, .{});
    defer file.close(io);

    var buffer: [65536]u8 = undefined;
    var writer = file.writer(io, &buffer);

    try writer.interface.print("{d} {d}\n", .{ A.nrow(), A.ncol() });

    for (0..A.nrow()) |i| for (0..A.ncol()) |j| {
        const sep = if (j == A.ncol() - 1) "\n" else " ";

        try writer.interface.print("{d:20.14}{s}", .{ A.at(i, j), sep });
    };

    try writer.interface.flush();
}

pub fn writeMatrixComplex(comptime T: type, io: std.Io, fname: []const u8, A: Matrix(T)) !void {
    var file = try std.Io.Dir.cwd().createFile(io, fname, .{});
    defer file.close(io);

    var buffer: [65536]u8 = undefined;
    var writer = file.writer(io, &buffer);

    try writer.interface.print("{d} {d}\n", .{ A.nrow(), 2 * A.ncol() });

    for (0..A.nrow()) |i| for (0..A.ncol()) |j| {
        const sep = if (j == A.ncol() - 1) "\n" else " ";

        try writer.interface.print("{d:20.14} {d:20.14}{s}", .{ A.at(i, j).re, A.at(i, j).im, sep });
    };

    try writer.interface.flush();
}
