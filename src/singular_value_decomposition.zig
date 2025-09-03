//! Singular Value Decomposition (SVD) implementation in Zig.

const std = @import("std");

const eigenproblem_solver = @import("eigenproblem_solver.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const real_matrix = @import("real_matrix.zig");

const RealMatrix = real_matrix.RealMatrix;

const eigensystemSymmetric = eigenproblem_solver.eigensystemSymmetric;
const mmRealRealTrans = matrix_multiplication.mmRealRealTrans;
const mmRealRealTransAlloc = matrix_multiplication.mmRealRealTransAlloc;
const mmRealTransRealAlloc = matrix_multiplication.mmRealTransRealAlloc;

/// Function to perform singular value decomposition.
pub fn svd(comptime T: type, U: *RealMatrix(T), S: *RealMatrix(T), Vt: *RealMatrix(T), A: RealMatrix(T)) !void {
    const m = A.rows; const n = A.cols;

    S.zero();

    var B = if (m > n) try mmRealTransRealAlloc(T, A, A) else try mmRealRealTransAlloc(T, A, A); defer B.deinit();
    var C = if (m < n) try mmRealTransRealAlloc(T, A, A) else try mmRealRealTransAlloc(T, A, A); defer C.deinit();

    var BJ = try RealMatrix(T).init(B.rows, B.cols, B.allocator); defer BJ.deinit();
    var CJ = try RealMatrix(T).init(C.rows, C.cols, B.allocator); defer CJ.deinit();

    try eigensystemSymmetric(T, &BJ, U, B); try eigensystemSymmetric(T, &CJ, Vt, C);

    for (0..BJ.rows) |i| {if (BJ.at(i, i) < 16 * std.math.floatEps(T)) return error.NotImplemented;}

    for (0..BJ.rows - 1) |_| for (0..BJ.rows - 1) |j| if (BJ.at(j, j) < BJ.at(j + 1, j + 1)) {
        std.mem.swap(T, BJ.ptr(j, j), BJ.ptr(j + 1, j + 1));
    };

    for (0..CJ.rows - 1) |_| for (0..CJ.rows - 1) |j| if (CJ.at(j, j) < CJ.at(j + 1, j + 1)) {

        std.mem.swap(T, CJ.ptr(j, j), CJ.ptr(j + 1, j + 1));

        for (0..CJ.rows) |k| {
            std.mem.swap(T, Vt.ptr(k, j), Vt.ptr(k, j + 1));
        }
    };

    for (0..S.rows) |i| S.ptr(i, i).* = std.math.sqrt(BJ.at(i, i));

    for (0..Vt.rows) |i| for (i + 1..Vt.cols) |j| {
        std.mem.swap(T, Vt.ptr(i, j), Vt.ptr(j, i));
    };

    for (0..U.rows) |i| {

        var row = U.row(i).asMatrix();

        mmRealRealTrans(T, &row, A, Vt.row(i).asMatrix());

        row.divs(S.at(i, i));
    }

    for (0..U.rows) |i| for (i + 1..U.cols) |j| {
        std.mem.swap(T, U.ptr(i, j), U.ptr(j, i));
    };
}
