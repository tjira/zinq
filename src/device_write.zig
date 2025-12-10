//! Functions for writing stuff to files.

const std = @import("std");

const classical_particle = @import("classical_particle.zig");
const complex_matrix = @import("complex_matrix.zig");
const complex_vector = @import("complex_vector.zig");
const global_variables = @import("global_variables.zig");
const real_matrix = @import("real_matrix.zig");
const real_tensor_four = @import("real_tensor_four.zig");
const real_tensor_three = @import("real_tensor_three.zig");
const real_vector = @import("real_vector.zig");

const ClassicalParticle = classical_particle.ClassicalParticle;
const ComplexMatrix = complex_matrix.ComplexMatrix;
const RealMatrix = real_matrix.RealMatrix;
const RealTensor3 = real_tensor_three.RealTensor3;
const RealTensor4 = real_tensor_four.RealTensor4;

const AN2SM = global_variables.AN2SM;
const A2AU = global_variables.A2AU;

const WRITE_BUFFER_SIZE = global_variables.WRITE_BUFFER_SIZE;

/// Exports the complex matrix to a file with the leftmost column being linspaced values from start to end.
pub fn exportComplexMatrixWithLinspacedLeftColumn(comptime T: type, path: []const u8, A: ComplexMatrix(T), start: T, end: T) !void {
    var file = try std.fs.cwd().createFile(path, .{}); defer file.close();

    try writeComplexMatrixWithLinspacedLeftColumn(T, file, A, start, end);
}

/// Exports the real matrix to a file.
pub fn exportRealMatrix(comptime T: type, path: []const u8, A: RealMatrix(T)) !void {
    var file = try std.fs.cwd().createFile(path, .{}); defer file.close();

    try writeRealMatrix(T, file, A);
}

/// Exports the real matrix to a file with the leftmost column being linspaced values from start to end.
pub fn exportRealMatrixWithLinspacedLeftColumn(comptime T: type, path: []const u8, A: RealMatrix(T), start: T, end: T) !void {
    var file = try std.fs.cwd().createFile(path, .{}); defer file.close();

    try writeRealMatrixWithLinspacedLeftColumn(T, file, A, start, end);
}

/// Exports the real 4th order tensor to a file.
pub fn exportRealTensorFour(comptime T: type, path: []const u8, A: RealTensor4(T)) !void {
    var file = try std.fs.cwd().createFile(path, .{}); defer file.close();

    try writeRealTensorFour(T, file, A);
}

/// Exports the real 3rd order tensor to a file as a PPM image.
pub fn exportRealTensorThreeAsPPM(comptime T: type, path: []const u8, A: RealTensor3(T)) !void {
    var file = try std.fs.cwd().createFile(path, .{}); defer file.close();

    try writeRealTensorThreeAsPPM(T, file, A);
}

/// Print the formatted line to the standard output.
pub fn print(comptime format: []const u8, args: anytype) !void {
    try write(std.fs.File.stdout(), format, args);
}

/// Print a classical particle into a terminal.
pub fn printClassicalParticleAsMolecule(comptime T: type, object: ClassicalParticle(T), comment: ?[]const u8) !void {
    try writeClassicalParticleAsMolecule(T, std.fs.File.stdout(), object, comment);
}

/// Print the formatted complex matrix to the standard output.
pub fn printComplexMatrix(comptime T: type, A: ComplexMatrix(T)) !void {
    try writeComplexMatrix(T, std.fs.File.stdout(), A);
}

/// Print the formatted json to the standard output.
pub fn printJson(object: anytype) !void {
    try writeJson(std.fs.File.stdout(), object);
}

/// Print the formatted real matrix to the standard output.
pub fn printRealMatrix(comptime T: type, A: RealMatrix(T)) !void {
    try writeRealMatrix(T, std.fs.File.stdout(), A);
}

/// Print a formatted line into the specified device.
pub fn write(device: std.fs.File, comptime format: []const u8, args: anytype) !void {
    var buffer: [WRITE_BUFFER_SIZE]u8 = undefined;

    var writer = device.writer(&buffer); var writer_interface = &writer.interface;

    try writer_interface.print(format, args);

    try writer_interface.flush();
}

/// Write the classical particle as a .xyz molecule file into the specified device with an optional comment line.
pub fn writeClassicalParticleAsMolecule(comptime T: type, device: std.fs.File, object: ClassicalParticle(T), comment: ?[]const u8) !void {
    var buffer: [WRITE_BUFFER_SIZE]u8 = undefined;

    var writer = device.writer(&buffer); var writer_interface = &writer.interface;

    if (comment != null) {
        try writer_interface.print("{d}\n{s}\n", .{object.atoms.?.len, comment.?});
    }

    for (0..object.atoms.?.len) |i| {

        const x = object.position.at(3 * i + 0) / A2AU;
        const y = object.position.at(3 * i + 1) / A2AU;
        const z = object.position.at(3 * i + 2) / A2AU;

        try writer_interface.print("{s} {d:20.14} {d:20.14} {d:20.14}\n", .{try AN2SM(object.atoms.?[i]), x, y, z});
    }

    try writer_interface.flush();
}

/// Print the formatted complex matrix to the specified device.
pub fn writeComplexMatrix(comptime T: type, device: std.fs.File, A: ComplexMatrix(T)) !void {
    var buffer: [WRITE_BUFFER_SIZE]u8 = undefined;

    var writer = device.writer(&buffer); var writer_interface = &writer.interface;

    try writer_interface.print("{d} {d}\n", .{A.rows, 2 * A.cols});

    for (0..A.rows) |i| for (0..A.cols) |j| {
        try writer_interface.print("{d:20.14} {d:20.14}{s}", .{A.at(i, j).re, A.at(i, j).im, if (j == A.cols - 1) "\n" else " "});
    };

    try writer_interface.flush();
}

/// Print the formatted complex matrix to the specified device.
pub fn writeComplexMatrixWithLinspacedLeftColumn(comptime T: type, device: std.fs.File, A: ComplexMatrix(T), start: T, end: T) !void {
    var buffer: [WRITE_BUFFER_SIZE]u8 = undefined;

    var writer = device.writer(&buffer); var writer_interface = &writer.interface;

    try writer_interface.print("{d} {d}\n", .{A.rows, 2 * A.cols});

    for (0..A.rows) |i| {

        const x = start + (end - start) * @as(T, @floatFromInt(i)) / @as(T, @floatFromInt(A.rows - 1));

        try writer_interface.print("{d:20.14}", .{x});

        for (0..A.cols) |j| {
            try writer_interface.print(" {d:20.14} {d:20.14}{s}", .{A.at(i, j).re, A.at(i, j).im, if (j == A.cols - 1) "\n" else " "});
        }
    }

    try writer_interface.flush();
}

/// Print a formatted json into the specified device.
pub fn writeJson(device: std.fs.File, object: anytype) !void {
    var buffer: [WRITE_BUFFER_SIZE]u8 = undefined;

    var writer = device.writer(&buffer); var writer_interface = &writer.interface;

    try std.json.Stringify.value(object, .{.whitespace = .indent_2}, writer_interface);

    try writer_interface.print("\n", .{});

    try writer_interface.flush();
}

/// Print the formatted real matrix to the specified device.
pub fn writeRealMatrix(comptime T: type, device: std.fs.File, A: RealMatrix(T)) !void {
    var buffer: [WRITE_BUFFER_SIZE]u8 = undefined;

    var writer = device.writer(&buffer); var writer_interface = &writer.interface;

    try writer_interface.print("{d} {d}\n", .{A.rows, A.cols});

    for (0..A.rows) |i| for (0..A.cols) |j| {
        try writer_interface.print("{d:20.14}{s}", .{A.at(i, j), if (j == A.cols - 1) "\n" else " "});
    };

    try writer_interface.flush();
}

/// Write the real matrix to the specified device with the leftmost column being linspaced values from start to end.
pub fn writeRealMatrixWithLinspacedLeftColumn(comptime T: type, device: std.fs.File, A: RealMatrix(T), start: T, end: T) !void {
    var buffer: [WRITE_BUFFER_SIZE]u8 = undefined;

    var writer = device.writer(&buffer); var writer_interface = &writer.interface;

    try writer_interface.print("{d} {d}\n", .{A.rows, A.cols + 1});

    for (0..A.rows) |i| {

        const x = start + (end - start) * @as(T, @floatFromInt(i)) / @as(T, @floatFromInt(A.rows - 1));

        try writer_interface.print("{d:20.14}", .{x});

        for (0..A.cols) |j| {
            try writer_interface.print(" {d:20.14}{s}", .{A.at(i, j), if (j == A.cols - 1) "\n" else ""});
        }
    }

    try writer_interface.flush();
}

/// Print the formatted real 4th order tensor to the specified device.
pub fn writeRealTensorFour(comptime T: type, device: std.fs.File, A: RealTensor4(T)) !void {
    var buffer: [WRITE_BUFFER_SIZE]u8 = undefined;

    var writer = device.writer(&buffer); var writer_interface = &writer.interface;

    try writer_interface.print("{d} {d} {d} {d}\n", .{A.shape[0], A.shape[1], A.shape[2], A.shape[3]});

    for (0..A.shape[0]) |i| for (0..A.shape[1]) |j| for (0..A.shape[2]) |k| for (0..A.shape[3]) |l| {
        try writer_interface.print("{d:20.14}{s}", .{A.at(i, j, k, l), if (k == A.shape[2] - 1 and l == A.shape[3] - 1) "\n" else " "});
    };

    try writer_interface.flush();
}

/// Print the formatted real 3rd order tensor to the specified device as a PPM image.
pub fn writeRealTensorThreeAsPPM(comptime T: type, device: std.fs.File, A: RealTensor3(T)) !void {
    var buffer: [WRITE_BUFFER_SIZE]u8 = undefined;

    var writer = device.writer(&buffer); var writer_interface = &writer.interface;

    try writer_interface.print("P3 {d} {d} 255\n{d}", .{A.shape[0], A.shape[1], A.data[0]});

    for (A.data[1..A.data.len]) |value| try writer_interface.print(" {d}", .{value});

    try writer_interface.print("\n", .{});

    try writer_interface.flush();
}
