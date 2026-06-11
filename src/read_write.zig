const std = @import("std");

const Matrix = @import("tensor.zig").Matrix;

pub fn printf(io: std.Io, comptime format: []const u8, args: anytype) !void {
    var buffer: [4096]u8 = undefined;
    var writer = std.Io.File.stdout().writer(io, &buffer);

    try writeAndFlush(&writer, format, args);
}

pub fn writeMatrixHjoin(comptime T: type, io: std.Io, fname: []const u8, A: Matrix(T), B: Matrix(T)) !void {
    std.debug.assert(A.nrow() == B.nrow());

    if (@typeInfo(T) == .@"struct") {
        return try writeMatrixHjoinComplex(T, io, fname, A, B);
    }

    try writeMatrixHjoinReal(T, io, fname, A, B);
}

pub fn writeMatrixHjoinReal(comptime T: type, io: std.Io, fname: []const u8, A: Matrix(T), B: Matrix(T)) !void {
    var file = try std.Io.Dir.cwd().createFile(io, fname, .{});
    defer file.close(io);

    var buffer: [65536]u8 = undefined;
    var writer = file.writer(io, &buffer);

    try writer.interface.print("{d} {d}\n", .{ A.nrow(), A.ncol() + B.ncol() });

    for (0..A.nrow()) |i| {
        for (0..A.ncol()) |j| {
            try writer.interface.print("{d:20.14} ", .{A.at(i, j)});
        }

        for (0..B.ncol()) |j| {
            const sep = if (j == B.ncol() - 1) "\n" else " ";

            try writer.interface.print("{d:20.14}{s}", .{ B.at(i, j), sep });
        }
    }

    try writer.interface.flush();
}

pub fn writeMatrixHjoinComplex(comptime T: type, io: std.Io, fname: []const u8, A: Matrix(T), B: Matrix(T)) !void {
    var file = try std.Io.Dir.cwd().createFile(io, fname, .{});
    defer file.close(io);

    var buffer: [65536]u8 = undefined;
    var writer = file.writer(io, &buffer);

    try writer.interface.print("{d} {d}\n", .{ A.nrow(), 2 * A.ncol() + 2 * B.ncol() });

    for (0..A.nrow()) |i| {
        for (0..A.ncol()) |j| {
            try writer.interface.print("{d:20.14} {d:20.14} ", .{ A.at(i, j).re, A.at(i, j).im });
        }

        for (0..B.ncol()) |j| {
            const sep = if (j == B.ncol() - 1) "\n" else " ";

            try writer.interface.print(" {d:20.14} {d:20.14}{s}", .{ B.at(i, j).re, B.at(i, j).im, sep });
        }
    }

    try writer.interface.flush();
}

pub fn writeMatrixLspace(comptime T: type, io: std.Io, fname: []const u8, A: Matrix(T), start: T, end: T) !void {
    if (@typeInfo(T) == .@"struct") {
        return try writeMatrixLspaceComplex(T, io, fname, A, start, end);
    }

    try writeMatrixLspaceReal(T, io, fname, A, start, end);
}

pub fn writeMatrixLspaceReal(comptime T: type, io: std.Io, fname: []const u8, A: Matrix(T), start: T, end: T) !void {
    var file = try std.Io.Dir.cwd().createFile(io, fname, .{});
    defer file.close(io);

    var buffer: [65536]u8 = undefined;
    var writer = file.writer(io, &buffer);

    try writer.interface.print("{d} {d}\n", .{ A.nrow(), A.ncol() + 1 });

    const dt: T = if (A.nrow() > 1) (end - start) / @as(T, @floatFromInt(A.nrow() - 1)) else 0;

    for (0..A.nrow()) |i| {
        try writer.interface.print("{d:20.14} ", .{start + dt * @as(T, @floatFromInt(i))});

        for (0..A.ncol()) |j| {
            const sep = if (j == A.ncol() - 1) "\n" else " ";

            try writer.interface.print("{d:20.14}{s}", .{ A.at(i, j), sep });
        }
    }

    try writer.interface.flush();
}

pub fn writeMatrixLspaceComplex(comptime T: type, io: std.Io, fname: []const u8, A: Matrix(T), start: T, end: T) !void {
    var file = try std.Io.Dir.cwd().createFile(io, fname, .{});
    defer file.close(io);

    var buffer: [65536]u8 = undefined;
    var writer = file.writer(io, &buffer);

    try writer.interface.print("{d} {d}\n", .{ A.nrow(), 2 * A.ncol() + 1 });

    const dt: T = if (A.nrow() > 1) (end - start) / @as(T, @floatFromInt(A.nrow() - 1)) else std.mem.zeroes(T);

    for (0..A.nrow()) |i| {
        try writer.interface.print("{d:20.14} ", .{start + dt * @as(T, @floatFromInt(i))});

        for (0..A.ncol()) |j| {
            const sep = if (j == A.ncol() - 1) "\n" else " ";

            try writer.interface.print(" {d:20.14} {d:20.14}{s}", .{ A.at(i, j).re, A.at(i, j).im, sep });
        }
    }

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

fn writeAndFlush(device: *std.Io.File.Writer, comptime format: []const u8, args: anytype) !void {
    try device.interface.print(format, args);

    try device.interface.flush();
}
