//! Wrapper for OpenBLAS, a high-performance implementation of BLAS and LAPACK.

const std = @import("std");

const config = @import("config");

const cblas = if (config.use_openblas) @cImport(@cInclude("cblas.h")) else struct {};
const lapacke = if (config.use_openblas) @cImport(@cInclude("lapacke.h")) else struct {};

const complex_matrix = @import("complex_matrix.zig");
const real_matrix = @import("real_matrix.zig");

const ComplexMatrix = complex_matrix.ComplexMatrix;
const RealMatrix = real_matrix.RealMatrix;

/// Calculate the matrix-matrix product C = A * B and store the result in C. The AT and BT flags control whether the matrices A and B are transposed or not.
pub fn dgemm(comptime T: type, C: *RealMatrix(T), A: RealMatrix(T), comptime at: bool, B: RealMatrix(T), comptime bt: bool) !void {
    if (comptime !config.use_openblas) @compileError("OPENBLAS SUPPORT IS DISABLED IN CONFIG");

    const m: i32 = @intCast(C.rows);
    const n: i32 = @intCast(C.cols);
    const k: i32 = @intCast(if (!at) A.cols else A.rows);

    const AT: c_uint = if (at) cblas.CblasTrans else cblas.CblasNoTrans;
    const BT: c_uint = if (bt) cblas.CblasTrans else cblas.CblasNoTrans;

    const LDA = if (!at) k else m;
    const LDB = if (!bt) n else k;

    cblas.cblas_dgemm(cblas.CblasRowMajor, AT, BT, m, n, k, 1.0, &A.data[0], LDA, &B.data[0], LDB, 0.0, &C.data[0], n);
}

/// Diagonalize a symmetric matrix A, returning the eigenvalues in J and the eigenvectors in C. The provided matrix is not modified.
pub fn dsyevd(comptime T: type, J: *RealMatrix(T), C: *RealMatrix(T), A: RealMatrix(T)) !void {
    if (comptime !config.use_openblas) @compileError("OPENBLAS SUPPORT IS DISABLED IN CONFIG");

    const n: i32 = @intCast(A.rows);
    try A.copyTo(C);
    J.fill(0);

    const info = lapacke.LAPACKE_dsyevd(lapacke.LAPACK_ROW_MAJOR, 'V', 'U', n, &C.data[0], n, &J.data[0]);

    for (0..J.rows) |i| std.mem.swap(T, J.ptr(i, i), J.ptr(0, i));

    if (info != 0) return error.NumericalError;
}

/// Calculate the matrix-matrix product C = A * B and store the result in C. The AH and BH flags control whether the matrices A and B are Hermitian transposed or not.
pub fn zgemm(comptime T: type, C: *ComplexMatrix(T), A: ComplexMatrix(T), ah: bool, B: ComplexMatrix(T), bh: bool) !void {
    if (comptime !config.use_openblas) @compileError("OPENBLAS SUPPORT IS DISABLED IN CONFIG");

    const m: i32 = @intCast(C.rows);
    const n: i32 = @intCast(C.cols);
    const k: i32 = @intCast(if (!ah) A.cols else A.rows);

    const AH: c_uint = if (ah) cblas.CblasConjTrans else cblas.CblasNoTrans;
    const BH: c_uint = if (bh) cblas.CblasConjTrans else cblas.CblasNoTrans;

    const LDA = if (!ah) k else m;
    const LDB = if (!bh) n else k;

    cblas.cblas_zgemm(cblas.CblasRowMajor, AH, BH, m, n, k, &std.math.Complex(T).init(1.0, 0.0), &A.data[0], LDA, &B.data[0], LDB, &std.math.Complex(T).init(0.0, 0.0), &C.data[0], n);
}
