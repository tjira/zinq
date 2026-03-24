//! Functions for aligning trajectories and systems.

const std = @import("std");

const device_write = @import("device_write.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");
const determinant = @import("determinant.zig");
const global_variables = @import("global_variables.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const singular_value_decomposition = @import("singular_value_decomposition.zig");
const string_manipulation = @import("string_manipulation.zig");

const determinant3x3 = determinant.determinant3x3;
const mmAlloc = matrix_multiplication.mmAlloc;
const mm = matrix_multiplication.mm;
const svdAlloc = singular_value_decomposition.svdAlloc;
const uncr = string_manipulation.uncr;

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const A2AU = global_variables.A2AU;
const SM2AN = global_variables.SM2AN;
const AN2M = global_variables.AN2M;
const U2AU = global_variables.U2AU;

/// Align trajectory.
pub fn alignTrajectory(comptime T: type, positions: RealMatrix(T), allocator: std.mem.Allocator) !RealMatrix(T) {
    var aligned = try positions.clone(allocator); errdefer aligned.deinit(allocator);

    for (0..aligned.rows) |i| {

        var means: [3]T = .{0} ** 3;

        for (0..aligned.cols / 3) |j| {
            means[0] += aligned.at(i, j * 3 + 0);
            means[1] += aligned.at(i, j * 3 + 1);
            means[2] += aligned.at(i, j * 3 + 2);
        }

        means[0] /= @as(T, @floatFromInt(aligned.cols / 3));
        means[1] /= @as(T, @floatFromInt(aligned.cols / 3));
        means[2] /= @as(T, @floatFromInt(aligned.cols / 3));

        for (0..aligned.cols / 3) |j| {
            aligned.ptr(i, j * 3 + 0).* -= means[0];
            aligned.ptr(i, j * 3 + 1).* -= means[1];
            aligned.ptr(i, j * 3 + 2).* -= means[2];
        }
    }

    var reference = aligned.row(0).asMatrix(); try reference.reshape(reference.rows / 3, 3);

    for (1..aligned.rows) |i| {

        var row = aligned.row(i).asMatrix(); try row.reshape(row.rows / 3, 3);

        const H = try mmAlloc(T, reference, true, row, false, allocator); defer H.deinit(allocator);

        const SVD = try svdAlloc(T, H, allocator); defer SVD.U.deinit(allocator); defer SVD.S.deinit(allocator); defer SVD.VT.deinit(allocator);

        var R = try mmAlloc(T, SVD.VT, true, SVD.U, true, allocator); defer R.deinit(allocator);

        var CORR = try RealMatrix(T).initZero(3, 3, allocator); defer CORR.deinit(allocator); CORR.identity();

        CORR.ptr(0, 0).* = std.math.sign(try determinant3x3(T, R));

        const CORRU = try mmAlloc(T, CORR, false, SVD.U, true, allocator); defer CORRU.deinit(allocator);

        try mm(T, &R, SVD.VT, true, CORRU, false);

        const row_aligned = try mmAlloc(T, row, false, R, false, allocator); defer row_aligned.deinit(allocator);

        try row_aligned.copyTo(&row);
    }

    return aligned;
}

/// Function to read a trajectory from an XYZ file and return it as a RealMatrix.
pub fn readTrajectoryFromXYZ(comptime T: type, path: []const u8, allocator: std.mem.Allocator) !struct{positions: RealMatrix(T), masses: RealVector(T)} {
    const file = std.fs.cwd().openFile(path, .{}) catch |err| {

        std.log.err("FILE '{s}' NOT FOUND", .{path});

        return err;
    };

    defer file.close();

    var buffer: [1024]u8 = undefined; var reader = file.reader(&buffer); var reader_interface = &reader.interface;

    const natom = try std.fmt.parseInt(u32, uncr(try reader_interface.peekDelimiterExclusive('\n')), 10);

    var atoms = try allocator.alloc(usize, natom); defer allocator.free(atoms);

    var positions = try RealMatrix(T).init(0, 3 * natom, allocator); errdefer positions.deinit(allocator);

    while (true) {

        const peek = reader_interface.takeDelimiterExclusive('\n') catch |err| if (err == error.EndOfStream) break else return err;

        if (natom != try std.fmt.parseInt(u32, uncr(peek), 10)) {

            std.log.err("INCORRECTLY FORMATTED TRAJECTORY FILE", .{});

            return error.InvalidInput;
        }

        reader_interface.toss(1); _ = try reader_interface.discardDelimiterInclusive('\n');

        try positions.addRow(allocator);

        for (0..natom) |i| {

            var it = std.mem.tokenizeAny(u8, try reader_interface.takeDelimiterExclusive('\n'), " "); 

            var next = it.next() orelse {

                std.log.err("INCORRECTLY FORMATTED TRAJECTORY FILE", .{});

                return error.InvalidInput;
            };

            if (positions.rows == 1) atoms[i] = SM2AN.get(next) orelse {

                std.log.err("UNKNOWN ATOMIC SYMBOL '{s}' IN TRAJECTORY FILE", .{next});

                return error.InvalidInput;
            };

            for (0..3) |j| {

                next = it.next() orelse {

                    std.log.err("INCORRECTLY FORMATTED TRAJECTORY FILE", .{});

                    return error.InvalidInput;
                };

                positions.ptr(positions.rows - 1, 3 * i + j).* = try std.fmt.parseFloat(T, uncr(next)) * A2AU;
            }

            reader_interface.toss(1);
        }
    }

    var masses = try RealVector(T).init(3 * natom, allocator); errdefer masses.deinit(allocator);

    for (0..3 * natom) |i| masses.ptr(i).* = AN2M[atoms[i / 3]] * U2AU;

    return .{.positions = positions, .masses = masses};
}
