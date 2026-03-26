//! File with determinant calculation functions.

const std = @import("std");

const eigenproblem_solver = @import("eigenproblem_solver.zig");
const real_matrix = @import("real_matrix.zig");
const singular_value_decomposition = @import("singular_value_decomposition.zig");

const RealMatrix = real_matrix.RealMatrix;

const eigensystemHermitianAlloc = eigenproblem_solver.eigensystemHermitianAlloc;
const svdAlloc = singular_value_decomposition.svdAlloc;

/// Calculate the determinant of a symmetric matrix using its eigenvalue decomposition. The matrix must be hermitian (symmetric in the real case).
pub fn determinantHermitian(comptime T: type, A: RealMatrix(T), allocator: std.mem.Allocator) !T {
    const AJC = try eigensystemHermitianAlloc(T, A, allocator); defer AJC.J.deinit(allocator); defer AJC.C.deinit(allocator);

    var det: T = 1;

    for (0..AJC.J.rows) |i| det *= AJC.J.at(i, i);

    return det;
}

/// Calculate the determinant of a matrix using its singular value decomposition.
pub fn determinantAbs(comptime T: type, A: RealMatrix(T), allocator: std.mem.Allocator) !T {
    const SVD = try svdAlloc(T, A, allocator); defer SVD.U.deinit(allocator); defer SVD.S.deinit(allocator); defer SVD.VT.deinit(allocator);

    var det: T = 1;

    for (0..SVD.S.rows) |i| det *= SVD.S.at(i, i);

    return det;
}

/// Calculate the natural logarithm of the determinant of a matrix using its singular value decomposition.
pub fn determinantLog(comptime T: type, A: RealMatrix(T), allocator: std.mem.Allocator) !T {
    const SVD = try svdAlloc(T, A, allocator); defer SVD.U.deinit(allocator); defer SVD.S.deinit(allocator); defer SVD.VT.deinit(allocator);

    var det: T = 0;

    for (0..SVD.S.rows) |i| det += std.math.log(T, std.math.e, SVD.S.at(i, i));

    return det;
}

/// Calculate thedeterminant of a 3x3 matrix.
pub fn determinant3x3(comptime T: type, A: RealMatrix(T)) !T {
    if (A.rows != 3 or A.cols != 3) {

        std.log.err("DETERMINANT3X3 IS ONLY DEFINED FOR 3x3 MATRICES, BUT GOT {d}x{d}", .{A.rows, A.cols});

        return error.InvalidInput;
    }

    const a = A.at(0, 0);
    const b = A.at(0, 1);
    const c = A.at(0, 2);
    const d = A.at(1, 0);
    const e = A.at(1, 1);
    const f = A.at(1, 2);
    const g = A.at(2, 0);
    const h = A.at(2, 1);
    const i = A.at(2, 2);

    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
}
