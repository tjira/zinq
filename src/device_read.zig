//! Functions for writing stuff to files.

const std = @import("std");

const real_matrix = @import("real_matrix.zig");
const string_manipulation = @import("string_manipulation.zig");

const RealMatrix = real_matrix.RealMatrix;

const uncr = string_manipulation.uncr;

/// Reads the real matrix from a file.
pub fn readRealMatrix(comptime T: type, path: []const u8, allocator: std.mem.Allocator) !RealMatrix(T) {
    var file = try std.fs.cwd().openFile(path, .{}); defer file.close();

    var buffer: [32768]u8 = undefined; var reader = file.reader(&buffer); var reader_interface = &reader.interface;

    const rowstr = try reader_interface.takeDelimiterExclusive(' ' ); const rows = try std.fmt.parseInt(usize, uncr(rowstr), 10);
    const colstr = try reader_interface.takeDelimiterExclusive('\n'); const cols = try std.fmt.parseInt(usize, uncr(colstr), 10);

    const A = try RealMatrix(T).init(rows, cols, allocator); var i: usize = 0;

    while (true) {

        const line = reader_interface.takeDelimiterExclusive('\n') catch {break;};

        var line_iterator = std.mem.tokenizeAny(u8, line, " "); 

        while (line_iterator.next()) |element| : (i += 1) {
            A.data[i] = try std.fmt.parseFloat(T, uncr(element));
        }
    }

    return A;
}
