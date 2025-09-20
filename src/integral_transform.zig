//! Molecular integral transform file.

const std = @import("std");

const real_matrix = @import("real_matrix.zig");
const real_tensor_four = @import("real_tensor_four.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealTensor4 = real_tensor_four.RealTensor4;

/// Extend one-electron integral matrices to spinorbital basis.
pub fn oneAO2AS(comptime T: type, A: *RealMatrix(T)) !void {
    var AS = try RealMatrix(T).initZero(2 * A.rows, 2 * A.cols, A.allocator);

    for (0..A.rows) |i| for (0..A.cols) |j| {
        AS.ptr(i         , j         ).* = A.at(i, j);
        AS.ptr(i + A.rows, j + A.cols).* = A.at(i, j);
    };

    A.deinit(); A.* = AS;
}

/// Extend two-electron integral tensor to spinorbital basis.
pub fn twoAO2AS(comptime T: type, A: *RealTensor4(T)) !void {
    var AS = try RealTensor4(T).initZero([4]usize{2 * A.shape[0], 2 * A.shape[1], 2 * A.shape[2], 2 * A.shape[3]}, A.allocator);

    for (0..A.shape[0]) |i| for (0..A.shape[1]) |j| for (0..A.shape[2]) |k| for (0..A.shape[3]) |l| {
        AS.ptr(i,              j,              k             , l             ).* = A.at(i, j, k, l);
        AS.ptr(i,              j             , k + A.shape[2], l + A.shape[3]).* = A.at(i, j, k, l);
        AS.ptr(i + A.shape[0], j + A.shape[1], k,              l             ).* = A.at(i, j, k, l);
        AS.ptr(i + A.shape[0], j + A.shape[1], k + A.shape[2], l + A.shape[3]).* = A.at(i, j, k, l);
    };

    A.deinit(); A.* = AS;
}
