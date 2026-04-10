//! Wrapper for OpenBLAS, a high-performance implementation of BLAS and LAPACK.

const std = @import("std");

const config = @import("config");

const lapacke = if (config.use_openblas) @cImport(@cInclude("lapacke.h")) else struct {};

const real_matrix = @import("real_matrix.zig");

const RealMatrix = real_matrix.RealMatrix;

/// Diagonalize a symmetric matrix A, returning the eigenvalues in J and the eigenvectors in C. The provided matrix is not modified.
pub fn dsyevd(comptime T: type, J: *RealMatrix(T), C: *RealMatrix(T), A: RealMatrix(T)) !void {
    if (comptime !config.use_openblas) @compileError("OPENBLAS SUPPORT IS DISABLED IN CONFIG");

    const n: i32 = @intCast(A.rows); try A.copyTo(C); J.fill(0);

    const info = lapacke.LAPACKE_dsyevd(lapacke.LAPACK_ROW_MAJOR, 'V', 'U', n, &C.data[0], n, &J.data[0]);

    for (0..J.rows) |i| std.mem.swap(T, J.ptr(i, i), J.ptr(0, i));

    if (info != 0) return error.NumericalError;
}
