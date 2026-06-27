const std = @import("std");

const Matrix = @import("tensor.zig").Matrix;

const isComplex = @import("value.zig").isComplex;
const primType = @import("value.zig").primType;

pub fn printf(io: std.Io, comptime format: []const u8, args: anytype) !void {
    var buffer: [4096]u8 = undefined;
    var writer = std.Io.File.stdout().writer(io, &buffer);

    try writer.interface.print(format, args);

    try writer.interface.flush();
}

pub fn writeMatrix(comptime T: type, io: std.Io, fname: []const u8, A: Matrix(T)) !void {
    var file = try std.Io.Dir.cwd().createFile(io, fname, .{});
    defer file.close(io);

    var buffer: [65536]u8 = undefined;
    var writer = file.writer(io, &buffer);

    const ncol = if (comptime isComplex(T)) 2 * A.ncol() else A.ncol();

    try writer.interface.print("{d} {d}\n", .{ A.nrow(), ncol });

    for (0..A.nrow()) |i| for (0..A.ncol()) |j| {
        try writeElement(&writer, A.at(i, j));

        try writer.interface.print("{s}", .{if (j == A.ncol() - 1) "\n" else " "});
    };

    try writer.interface.flush();
}

pub fn writeMatrixHjoin(comptime T: type, io: std.Io, fname: []const u8, A: Matrix(T), nca: ?usize, B: Matrix(T), ncb: ?usize) !void {
    std.debug.assert(A.nrow() == B.nrow());

    var file = try std.Io.Dir.cwd().createFile(io, fname, .{});
    defer file.close(io);

    var buffer: [65536]u8 = undefined;
    var writer = file.writer(io, &buffer);

    const cols_A = if (nca) |limit| limit else (if (comptime isComplex(T)) 2 * A.ncol() else A.ncol());
    const cols_B = if (ncb) |limit| limit else (if (comptime isComplex(T)) 2 * B.ncol() else B.ncol());

    try writer.interface.print("{d} {d}\n", .{ A.nrow(), cols_A + cols_B });

    const a_cols = if (comptime isComplex(T)) cols_A / 2 else cols_A;
    const b_cols = if (comptime isComplex(T)) cols_B / 2 else cols_B;

    for (0..A.nrow()) |i| {
        for (0..a_cols) |j| {
            try writeElement(&writer, A.at(i, j));

            try writer.interface.print(" ", .{});
        }

        for (0..b_cols) |j| {
            try writeElement(&writer, B.at(i, j));

            try writer.interface.print("{s}", .{if (j == b_cols - 1) "\n" else " "});
        }
    }

    try writer.interface.flush();
}

pub fn writeMatrixLspace(comptime T: type, io: std.Io, fname: []const u8, A: Matrix(T), start: primType(T), end: primType(T)) !void {
    var file = try std.Io.Dir.cwd().createFile(io, fname, .{});
    defer file.close(io);

    var buffer: [65536]u8 = undefined;
    var writer = file.writer(io, &buffer);

    const ncol = if (comptime isComplex(T)) 2 * A.ncol() else A.ncol();

    try writer.interface.print("{d} {d}\n", .{ A.nrow(), ncol + 1 });

    const dt = if (A.nrow() > 1) (end - start) / @as(primType(T), @floatFromInt(A.nrow() - 1)) else 0;

    for (0..A.nrow()) |i| {
        try writer.interface.print("{d:20.14} ", .{start + dt * @as(primType(T), @floatFromInt(i))});

        for (0..A.ncol()) |j| {
            try writeElement(&writer, A.at(i, j));

            try writer.interface.print("{s}", .{if (j == A.ncol() - 1) "\n" else " "});
        }
    }

    try writer.interface.flush();
}

fn writeElement(writer: anytype, val: anytype) !void {
    if (comptime isComplex(@TypeOf(val))) {
        try writer.interface.print("{d:20.14} {d:20.14}", .{ val.re, val.im });
    }

    if (comptime !isComplex(@TypeOf(val))) {
        try writer.interface.print("{d:20.14}", .{val});
    }
}
