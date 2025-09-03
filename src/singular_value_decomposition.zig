//! Singular Value Decomposition (SVD) implementation in Zig.

const std = @import("std");

const eigenproblem_solver = @import("eigenproblem_solver.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const mmRealTransReal = matrix_multiplication.mmRealTransReal;
const mmRealRealTrans = matrix_multiplication.mmRealRealTrans;
const diagonalizeSymmetric = eigenproblem_solver.diagonalizeSymmetric;
const eigensystemSymmetric = eigenproblem_solver.eigensystemSymmetric;

pub fn svd(comptime T: type, U: *RealMatrix(T), S: *RealMatrix(T), Vt: *RealMatrix(T), A: RealMatrix(T)) !void {
    _ = U; _ = S; _ = Vt; _ = A;
}
