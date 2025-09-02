//! Functions for writing stuff to files.

const std = @import("std");

const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

/// Exports the real matrix to a file.
pub fn exportRealMatrix(comptime T: type, path: []const u8, A: RealMatrix(T)) !void {
    var file = try std.fs.cwd().createFile(path, .{}); defer file.close();

    try writeRealMatrix(T, file, A);
}

/// Exports the real matrix to a file with the leftmost column being linspaced values from start to end with points number of points.
pub fn exportRealMatrixWithLinspacedLeftColumn(comptime T: type, path: []const u8, A: RealMatrix(T), start: T, end: T, points: usize) !void {
    var file = try std.fs.cwd().createFile(path, .{}); defer file.close();

    try writeRealMatrixWithLinspacedLeftColumn(T, file, A, start, end, points);
}

/// Print the formatted line to the standard output.
pub fn print(comptime format: []const u8, args: anytype) !void {
    try write(std.fs.File.stdout(), format, args);
}

/// Print the formatted real matrix to the standard output.
pub fn printRealMatrix(comptime T: type, A: RealMatrix(T)) !void {
    try writeRealMatrix(T, std.fs.File.stdout(), A);
}

/// Print the formatted real vector to the standard output.
pub fn printRealVector(comptime T: type, v: RealVector(T)) !void {
    try writeRealVector(T, std.fs.File.stdout(), v);
}

/// Print a formatted line into the specified device.
pub fn write(device: std.fs.File, comptime format: []const u8, args: anytype) !void {
    var buffer: [32768]u8 = undefined;

    var writer = device.writer(&buffer); var writer_interface = &writer.interface;

    try writer_interface.print(format, args);

    try writer_interface.flush();
}

/// Print the formatted real matrix to the specified device.
pub fn writeRealMatrix(comptime T: type, device: std.fs.File, A: RealMatrix(T)) !void {
    var buffer: [32768]u8 = undefined;

    var writer = device.writer(&buffer); var writer_interface = &writer.interface;

    try writer_interface.print("{d} {d}\n", .{A.rows, A.cols});

    for (0..A.rows) |i| for (0..A.cols) |j| {
        try writer_interface.print("{d:20.14}{s}", .{A.at(i, j), if(j == A.cols - 1) "\n" else " "});
    };

    try writer_interface.flush();
}

/// Write the real matrix to the specified device with the leftmost column being linspaced values from start to end with points number of points.
pub fn writeRealMatrixWithLinspacedLeftColumn(comptime T: type, device: std.fs.File, A: RealMatrix(T), start: T, end: T, points: usize) !void {
    var buffer: [32768]u8 = undefined;

    var writer = device.writer(&buffer); var writer_interface = &writer.interface;

    try writer_interface.print("{d} {d}\n", .{A.rows, A.cols + 1});

    for (0..A.rows) |i| {

        const x = start + (end - start) * @as(T, @floatFromInt(i)) / @as(T, @floatFromInt(points - 1));

        try writer_interface.print("{d:20.14}", .{x});

        for (0..A.cols) |j| {
            try writer_interface.print(" {d:20.14}{s}", .{A.at(i, j), if(j == A.cols - 1) "\n" else ""});
        }
    }

    try writer_interface.flush();
}

/// Print the formatted real vector to the specified device.
pub fn writeRealVector(comptime T: type, device: std.fs.File, v: RealVector(T)) !void {
    var buffer: [32768]u8 = undefined;

    var writer = device.writer(&buffer); var writer_interface = &writer.interface;

    try writer_interface.print("{d}\n", .{v.len});

    for (0..v.len) |i| {
        try writer_interface.print("{d:20.14}\n", .{v.at(i)});
    }

    try writer_interface.flush();
}
