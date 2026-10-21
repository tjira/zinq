//! Utility functions for reading and writing matrix and tensor data, formatted printing, and molecular coordinates.

const std = @import("std");

const Complex = std.math.Complex;

const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;

const getSymbol = @import("constant.zig").getSymbol;
const isComplex = @import("value.zig").isComplex;
const primType = @import("value.zig").primType;

const A2BOHR = @import("constant.zig").A2BOHR;
const AMU2AU = @import("constant.zig").AMU2AU;
const AU2FS = @import("constant.zig").AU2FS;

/// Resolves the trajectory output path with automatic numerical suffixing when ntraj > 1.
pub fn getTrajectoryPath(path: []const u8, i: usize, ntraj: usize, gpa: std.mem.Allocator) ![]const u8 {
    if (ntraj <= 1) return path;

    const ext = std.fs.path.extension(path);

    const width = std.math.log10_int(ntraj) + 1;

    var num_buf: [32]u8 = undefined;

    const len = std.fmt.printInt(&num_buf, i + 1, 10, .lower, .{ .width = width, .fill = '0' });

    return try std.fmt.allocPrint(gpa, "{s}_{s}{s}", .{ path[0 .. path.len - ext.len], num_buf[0..len], ext });
}

/// Formats and prints a string to standard output, flushing the writer buffer.
pub fn printf(io: std.Io, comptime format: []const u8, args: anytype) !void {
    var buffer: [4096]u8 = undefined;
    var writer = std.Io.File.stdout().writer(io, &buffer);

    try writer.interface.print(format, args);

    try writer.interface.flush();
}

/// Prints complex matrix dimensions and space-separated elements to standard output.
pub fn printComplexMatrix(comptime T: type, io: std.Io, A: Matrix(Complex(primType(T)))) !void {
    try printComplexTensor(T, 2, io, A.asTensor());
}

/// Prints multi-dimensional complex tensor dimensions and space-separated elements to standard output.
pub fn printComplexTensor(comptime T: type, comptime N: usize, io: std.Io, A: Tensor(Complex(primType(T)), N)) !void {
    var buffer: [65536]u8 = undefined;
    var writer = std.Io.File.stdout().writer(io, &buffer);

    try writeComplexTensorWriter(T, N, &writer, A);
}

/// Formats and prints atomic symbols and 3D nuclear coordinates in angstroms to standard output.
pub fn printGeometry(comptime T: type, io: std.Io, atoms: []const i32, coors: []const T, title: ?[]const u8) !void {
    const header = title orelse "MOLECULAR GEOMETRY (Å)";

    try printf(io, "\n{s}\n", .{header});

    for (0..atoms.len) |i| {
        const sym = getSymbol(atoms[i]) catch "X";

        const x = coors[3 * i + 0] / A2BOHR;
        const y = coors[3 * i + 1] / A2BOHR;
        const z = coors[3 * i + 2] / A2BOHR;

        try printf(io, "{s:2} {d:20.14} {d:20.14} {d:20.14}\n", .{ sym, x, y, z });
    }
}

/// Prints real or complex matrix dimensions and elements to standard output.
pub fn printMatrix(comptime T: type, io: std.Io, A: Matrix(T)) !void {
    if (comptime isComplex(T)) {
        try printComplexMatrix(primType(T), io, A);
    }

    if (comptime !isComplex(T)) {
        try printRealMatrix(T, io, A);
    }
}

/// Prints real matrix dimensions and space-separated elements to standard output.
pub fn printRealMatrix(comptime T: type, io: std.Io, A: Matrix(T)) !void {
    try printRealTensor(T, 2, io, A.asTensor());
}

/// Prints multi-dimensional real tensor dimensions and space-separated elements to standard output.
pub fn printRealTensor(comptime T: type, comptime N: usize, io: std.Io, A: Tensor(T, N)) !void {
    var buffer: [65536]u8 = undefined;
    var writer = std.Io.File.stdout().writer(io, &buffer);

    try writeRealTensorWriter(T, N, &writer, A);
}

/// Prints multi-dimensional real or complex tensor dimensions and elements to standard output.
pub fn printTensor(comptime T: type, comptime N: usize, io: std.Io, A: Tensor(T, N)) !void {
    if (comptime isComplex(T)) {
        try printComplexTensor(primType(T), N, io, A);
    }

    if (comptime !isComplex(T)) {
        try printRealTensor(T, N, io, A);
    }
}

/// Reads a complex matrix from a text file, parsing dimensions and complex elements.
pub fn readComplexMatrix(comptime T: type, io: std.Io, path: []const u8, allocator: std.mem.Allocator) !Matrix(Complex(primType(T))) {
    const ten = try readComplexTensor(T, 2, io, path, allocator);

    return ten.asMatrix();
}

/// Reads a complex multi-dimensional tensor from a file, parsing header shape and elements.
pub fn readComplexTensor(comptime T: type, comptime N: usize, io: std.Io, path: []const u8, allocator: std.mem.Allocator) !Tensor(Complex(primType(T)), N) {
    var file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);

    var buffer: [65536]u8 = undefined;
    var reader = file.reader(io, &buffer);

    var hdr_it = std.mem.tokenizeAny(u8, try reader.interface.takeDelimiterInclusive('\n'), " \t\r\n");

    var shape: [N]usize = undefined;

    inline for (0..N) |d| {
        shape[d] = try std.fmt.parseInt(usize, hdr_it.next() orelse return error.InvalidRank, 10);
    }

    if (hdr_it.next() != null) return error.InvalidRank;

    var A = try Tensor(Complex(primType(T)), N).init(shape, allocator);
    errdefer A.deinit(allocator);

    var i: usize, var pending_re: ?primType(T) = .{ 0, null };

    while (true) {
        const line = reader.interface.takeDelimiterExclusive('\n') catch break;

        reader.interface.toss(1);

        var line_it = std.mem.tokenizeAny(u8, std.mem.trim(u8, line, "\r"), " \t");

        while (line_it.next()) |elem| {
            const val = try std.fmt.parseFloat(primType(T), elem);

            if (pending_re) |re| {
                if (i >= A.data.len) return error.InvalidDataLength;

                A.data[i] = Complex(primType(T)).init(re, val);

                i, pending_re = .{ i + 1, null };
            } else {
                pending_re = val;
            }
        }
    }

    if (pending_re != null or i != A.data.len) return error.IncompleteData;

    return A;
}

/// Reads a real or complex matrix from a text file, parsing dimensions and elements.
pub fn readMatrix(comptime T: type, io: std.Io, path: []const u8, allocator: std.mem.Allocator) !Matrix(T) {
    if (comptime isComplex(T)) {
        return try readComplexMatrix(primType(T), io, path, allocator);
    }

    if (comptime !isComplex(T)) {
        return try readRealMatrix(T, io, path, allocator);
    }
}

/// Reads a real matrix from a text file, parsing dimensions and floating-point elements.
pub fn readRealMatrix(comptime T: type, io: std.Io, path: []const u8, allocator: std.mem.Allocator) !Matrix(T) {
    const ten = try readRealTensor(T, 2, io, path, allocator);

    return ten.asMatrix();
}

/// Reads a multi-dimensional real tensor from a text file, parsing shape and elements.
pub fn readRealTensor(comptime T: type, comptime N: usize, io: std.Io, path: []const u8, allocator: std.mem.Allocator) !Tensor(T, N) {
    var file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);

    var buffer: [65536]u8 = undefined;
    var reader = file.reader(io, &buffer);

    var hdr_it = std.mem.tokenizeAny(u8, try reader.interface.takeDelimiterInclusive('\n'), " \t\r\n");

    var shape: [N]usize = undefined;

    inline for (0..N) |d| {
        shape[d] = try std.fmt.parseInt(usize, hdr_it.next() orelse return error.InvalidRank, 10);
    }

    if (hdr_it.next() != null) return error.InvalidRank;

    var A = try Tensor(T, N).init(shape, allocator);
    errdefer A.deinit(allocator);

    var i: usize = 0;

    while (true) {
        const line = reader.interface.takeDelimiterExclusive('\n') catch break;

        reader.interface.toss(1);

        var line_it = std.mem.tokenizeAny(u8, std.mem.trim(u8, line, "\r"), " \t");

        while (line_it.next()) |elem| : (i += 1) {
            if (i >= A.data.len) return error.InvalidDataLength;

            A.data[i] = try std.fmt.parseFloat(T, elem);
        }
    }

    if (i != A.data.len) return error.IncompleteData;

    return A;
}

/// Reads a multi-dimensional real or complex tensor from a text file.
pub fn readTensor(comptime T: type, comptime N: usize, io: std.Io, path: []const u8, allocator: std.mem.Allocator) !Tensor(T, N) {
    if (comptime isComplex(T)) {
        return try readComplexTensor(primType(T), N, io, path, allocator);
    }

    if (comptime !isComplex(T)) {
        return try readRealTensor(T, N, io, path, allocator);
    }
}

/// Reads the header line of a tensor data file and determines the dimensional rank.
pub fn readTensorRank(io: std.Io, path: []const u8) !usize {
    var file = try std.Io.Dir.cwd().openFile(io, path, .{});
    defer file.close(io);

    var buffer: [65536]u8 = undefined;
    var reader = file.reader(io, &buffer);

    const hdr_line = try reader.interface.takeDelimiterInclusive('\n');
    var hdr_it = std.mem.tokenizeAny(u8, hdr_line, " \t\r\n");

    var rank: usize = 0;

    while (hdr_it.next()) |_| rank += 1;

    if (rank == 0) return error.InvalidFormat;

    return rank;
}

/// Writes complex matrix dimensions and elements to a file in standardized representation.
pub fn writeComplexMatrix(comptime T: type, io: std.Io, fname: []const u8, A: Matrix(Complex(primType(T)))) !void {
    try writeComplexTensor(T, 2, io, fname, A.asTensor());
}

/// Writes a multi-dimensional complex tensor to a file, serializing shape and elements.
pub fn writeComplexTensor(comptime T: type, comptime N: usize, io: std.Io, fname: []const u8, A: Tensor(Complex(primType(T)), N)) !void {
    var file = try std.Io.Dir.cwd().createFile(io, fname, .{});
    defer file.close(io);

    var buffer: [65536]u8 = undefined;
    var writer = file.writer(io, &buffer);

    try writeComplexTensorWriter(T, N, &writer, A);
}

/// Writes matrix dimensions and space-separated elements to a file, handling real or complex.
pub fn writeMatrix(comptime T: type, io: std.Io, fname: []const u8, A: Matrix(T)) !void {
    if (comptime isComplex(T)) {
        try writeComplexMatrix(primType(T), io, fname, A);
    }

    if (comptime !isComplex(T)) {
        try writeRealMatrix(T, io, fname, A);
    }
}

/// Horizontally concatenates two matrices and writes the merged matrix to a file.
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

/// Writes a matrix to a file prepended with an evenly spaced coordinate/time grid column.
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

/// Writes real matrix dimensions and space-separated elements to a file.
pub fn writeRealMatrix(comptime T: type, io: std.Io, fname: []const u8, A: Matrix(T)) !void {
    try writeRealTensor(T, 2, io, fname, A.asTensor());
}

/// Writes a multi-dimensional real tensor to a file, serializing the shape and elements.
pub fn writeRealTensor(comptime T: type, comptime N: usize, io: std.Io, fname: []const u8, A: Tensor(T, N)) !void {
    var file = try std.Io.Dir.cwd().createFile(io, fname, .{});
    defer file.close(io);

    var buffer: [65536]u8 = undefined;
    var writer = file.writer(io, &buffer);

    try writeRealTensorWriter(T, N, &writer, A);
}

/// Writes a multi-dimensional real or complex tensor to a file.
pub fn writeTensor(comptime T: type, comptime N: usize, io: std.Io, fname: []const u8, A: Tensor(T, N)) !void {
    if (comptime isComplex(T)) {
        try writeComplexTensor(primType(T), N, io, fname, A);
    }

    if (comptime !isComplex(T)) {
        try writeRealTensor(T, N, io, fname, A);
    }
}

/// Writes atomic symbols and 3D nuclear coordinates to a file in standard XYZ format.
pub fn writeXyzFile(comptime T: type, io: std.Io, fname: []const u8, atoms: []const i32, coors: []const T) !void {
    var file = try std.Io.Dir.cwd().createFile(io, fname, .{});
    defer file.close(io);

    var buffer: [65536]u8 = undefined;
    var writer = file.writer(io, &buffer);

    try writer.interface.print("{d}\n\n", .{atoms.len});

    for (0..atoms.len) |i| {
        const sym = getSymbol(atoms[i]) catch "X";

        const x = coors[3 * i + 0] / A2BOHR;
        const y = coors[3 * i + 1] / A2BOHR;
        const z = coors[3 * i + 2] / A2BOHR;

        try writer.interface.print("{s:2} {d:20.14} {d:20.14} {d:20.14}\n", .{ sym, x, y, z });
    }

    try writer.interface.flush();
}

/// Appends a single XYZ frame to a writer stream with atomic coordinates, timestamp, and optional electronic state.
pub fn writeXyzFrame(comptime T: type, writer: anytype, atoms: []const i32, coors: []const T, time: T, state: ?usize) !void {
    if (state) |s| {
        try writer.interface.print("{d}\nTIME = {d:.4} a.u. = {d:.4} fs, STATE = {d}\n", .{ atoms.len, time, time * AU2FS, s });
    }

    if (state == null) {
        try writer.interface.print("{d}\nTIME = {d:.4} a.u. = {d:.4} fs\n", .{ atoms.len, time, time * AU2FS });
    }

    for (0..atoms.len) |i| {
        const sym = getSymbol(atoms[i]) catch "X";

        const x = coors[3 * i + 0] / A2BOHR;
        const y = coors[3 * i + 1] / A2BOHR;
        const z = coors[3 * i + 2] / A2BOHR;

        try writer.interface.print("{s:2} {d:20.14} {d:20.14} {d:20.14}\n", .{ sym, x, y, z });
    }
}

/// Appends XYZ frame to a writer stream with atomic momentum vectors, timestamp, and optional electronic state.
pub fn writeXyzMomentumFrame(comptime T: type, writer: anytype, atoms: []const i32, mom: []const T, time: T, state: ?usize) !void {
    if (state) |s| {
        try writer.interface.print("{d}\nTIME = {d:.4} a.u. = {d:.4} fs, STATE = {d}\n", .{ atoms.len, time, time * AU2FS, s });
    }

    if (state == null) {
        try writer.interface.print("{d}\nTIME = {d:.4} a.u. = {d:.4} fs\n", .{ atoms.len, time, time * AU2FS });
    }

    const conv = 1.0 / (AMU2AU * A2BOHR * AU2FS);

    for (0..atoms.len) |i| {
        const sym = getSymbol(atoms[i]) catch "X";

        const px = mom[3 * i + 0] * conv;
        const py = mom[3 * i + 1] * conv;
        const pz = mom[3 * i + 2] * conv;

        try writer.interface.print("{s:2} {d:20.14} {d:20.14} {d:20.14}\n", .{ sym, px, py, pz });
    }
}

/// Writes trajectory momentum vectors across all simulation time steps to a file in standard multi-frame XYZ format.
pub fn writeXyzMomentumTrajectory(comptime T: type, io: std.Io, fname: []const u8, atoms: []const i32, P: Matrix(T), dt: primType(T)) !void {
    var file = try std.Io.Dir.cwd().createFile(io, fname, .{});
    defer file.close(io);

    var buffer: [65536]u8 = undefined;
    var writer = file.writer(io, &buffer);

    for (0..P.nrow()) |i| {
        const time = dt * @as(primType(T), @floatFromInt(i));

        try writeXyzMomentumFrame(T, &writer, atoms, P.rowSlice(i), time, null);
    }

    try writer.interface.flush();
}

/// Writes trajectory coordinates across all simulation time steps to a file in standard multi-frame XYZ format.
pub fn writeXyzTrajectory(comptime T: type, io: std.Io, fname: []const u8, atoms: []const i32, R: Matrix(T), dt: primType(T)) !void {
    var file = try std.Io.Dir.cwd().createFile(io, fname, .{});
    defer file.close(io);

    var buffer: [65536]u8 = undefined;
    var writer = file.writer(io, &buffer);

    for (0..R.nrow()) |i| {
        const time = dt * @as(primType(T), @floatFromInt(i));

        try writeXyzFrame(T, &writer, atoms, R.rowSlice(i), time, null);
    }

    try writer.interface.flush();
}

/// Writes multi-dimensional complex tensor dimensions and elements to a buffered writer stream.
fn writeComplexTensorWriter(comptime T: type, comptime N: usize, writer: anytype, A: Tensor(Complex(primType(T)), N)) !void {
    inline for (0..N) |d| {
        try writer.interface.print("{d}{s}", .{ A.shape[d], if (d == N - 1) "\n" else " " });
    }

    const mat = A.asMatrix();

    for (0..mat.nrow()) |i| for (0..mat.ncol()) |j| {
        try writeElement(writer, mat.at(i, j));

        try writer.interface.print("{s}", .{if (j == mat.ncol() - 1) "\n" else " "});
    };

    try writer.interface.flush();
}

/// Writes a single real or complex floating-point element to the output writer stream.
fn writeElement(writer: anytype, val: anytype) !void {
    if (comptime isComplex(@TypeOf(val))) {
        try writer.interface.print("{d:20.14} {d:20.14}", .{ val.re, val.im });
    }

    if (comptime !isComplex(@TypeOf(val))) {
        try writer.interface.print("{d:20.14}", .{val});
    }
}

/// Writes multi-dimensional real tensor dimensions and elements to a buffered writer stream.
fn writeRealTensorWriter(comptime T: type, comptime N: usize, writer: anytype, A: Tensor(T, N)) !void {
    inline for (0..N) |d| {
        try writer.interface.print("{d}{s}", .{ A.shape[d], if (d == N - 1) "\n" else " " });
    }

    const mat = A.asMatrix();

    for (0..mat.nrow()) |i| for (0..mat.ncol()) |j| {
        try writeElement(writer, mat.at(i, j));

        try writer.interface.print("{s}", .{if (j == mat.ncol() - 1) "\n" else " "});
    };

    try writer.interface.flush();
}
