//! Functions for writing stuff to files.

const std = @import("std");

const complex_matrix = @import("complex_matrix.zig");
const global_variables = @import("global_variables.zig");
const real_matrix = @import("real_matrix.zig");
const string_manipulation = @import("string_manipulation.zig");

const ComplexMatrix = complex_matrix.ComplexMatrix;
const RealMatrix = real_matrix.RealMatrix;

const uncr = string_manipulation.uncr;

const WRITE_BUFFER_SIZE = global_variables.WRITE_BUFFER_SIZE;

/// Reads the complex matrix from a file.
pub fn readComplexMatrix(comptime T: type, path: []const u8, allocator: std.mem.Allocator) !ComplexMatrix(T) {
    var file = try std.fs.cwd().openFile(path, .{}); defer file.close();

    var buffer: [WRITE_BUFFER_SIZE]u8 = undefined; var reader = file.reader(&buffer); var reader_interface = &reader.interface;

    const rowstr = try reader_interface.takeDelimiterExclusive(' ' ); const rows = try std.fmt.parseInt(usize, uncr(rowstr), 10); reader_interface.toss(1);
    const colstr = try reader_interface.takeDelimiterExclusive('\n'); const cols = try std.fmt.parseInt(usize, uncr(colstr), 10); reader_interface.toss(1);

    const A = try ComplexMatrix(T).init(rows, cols, allocator); var i: usize = 0;

    while (true) {

        const line = reader_interface.takeDelimiterExclusive('\n') catch {break;}; reader_interface.toss(1);

        var line_iterator = std.mem.tokenizeAny(u8, line, " "); 

        while (line_iterator.next()) |element| : (i += 1) {
            if (i % 2 == 0) A.data[i / 2].re = try std.fmt.parseFloat(T, uncr(element));
            if (i % 2 == 1) A.data[i / 2].im = try std.fmt.parseFloat(T, uncr(element));
        }
    }

    return A;
}

/// Reads the real matrix from a file.
pub fn readRealMatrix(comptime T: type, path: []const u8, allocator: std.mem.Allocator) !RealMatrix(T) {
    var file = try std.fs.cwd().openFile(path, .{}); defer file.close();

    var buffer: [WRITE_BUFFER_SIZE]u8 = undefined; var reader = file.reader(&buffer); var reader_interface = &reader.interface;

    const rowstr = try reader_interface.takeDelimiterExclusive(' ' ); const rows = try std.fmt.parseInt(usize, uncr(rowstr), 10); reader_interface.toss(1);
    const colstr = try reader_interface.takeDelimiterExclusive('\n'); const cols = try std.fmt.parseInt(usize, uncr(colstr), 10); reader_interface.toss(1);

    const A = try RealMatrix(T).init(rows, cols, allocator); var i: usize = 0;

    while (true) {

        const line = reader_interface.takeDelimiterExclusive('\n') catch {break;}; reader_interface.toss(1);

        var line_iterator = std.mem.tokenizeAny(u8, line, " "); 

        while (line_iterator.next()) |element| : (i += 1) {
            A.data[i] = try std.fmt.parseFloat(T, uncr(element));
        }
    }

    return A;
}
