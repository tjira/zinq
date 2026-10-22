//! BLAS and LAPACK wrapper for matrix operations, eigendecompositions, and linear system solvers.

const std = @import("std");

const cblas = @import("cimport.zig").cblas;
const lapacke = @import("cimport.zig").lapacke;

const Matrix = @import("tensor.zig").Matrix;
const Vector = @import("tensor.zig").Vector;

const isComplex = @import("value.zig").isComplex;
const primType = @import("value.zig").primType;

const ROW_MAJOR = lapacke.LAPACK_ROW_MAJOR;

/// Scales a vector and adds it to another: y = alpha * x + y (AXPY operation).
pub fn axpy(comptime T: type, y: *Vector(T), x: Vector(T), alpha: T) void {
    std.debug.assert(x.length() == y.length());

    axpySlice(T, y.data, x.data, alpha);
}

/// Computes the determinant of a square matrix from its LU factorization diagonal and pivot parity.
pub fn det(comptime T: type, A: *Matrix(T), ipiv: []i32) !T {
    std.debug.assert(A.nrow() == A.ncol());
    std.debug.assert(ipiv.len == A.nrow());

    if (comptime primType(T) != f64) @compileError("DET NOW ONLY SUPPORTS F64 NUMBERS");

    const n = A.nrow();

    if (n == 0) {
        return 1;
    }

    if (n == 1) {
        return A.at(0, 0);
    }

    if (n == 2) {
        return A.at(0, 0) * A.at(1, 1) - A.at(0, 1) * A.at(1, 0);
    }

    const info = lapacke.LAPACKE_dgetrf(ROW_MAJOR, @intCast(n), @intCast(n), A.data.ptr, @intCast(n), ipiv.ptr);

    if (info < 0) {
        return error.LapackError;
    }

    if (info > 0) return 0;

    var d: T = 1;

    for (0..n) |k| {
        if (ipiv[k] != @as(i32, @intCast(k + 1))) d = -d;

        d *= A.at(k, k);
    }

    return d;
}

/// Computes the dot product (inner product) of two vectors: x^T * y using CBLAS ddot.
pub fn dot(comptime T: type, x: Vector(T), y: Vector(T)) T {
    std.debug.assert(x.length() == y.length());

    return dotSlice(T, x.data, y.data);
}

/// Performs general matrix-matrix multiplication: C = alpha * op(A) * op(B) + beta * C (GEMM).
pub fn gemm(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T), alpha: T, beta: T, ta: bool, tb: bool) void {
    const m = if (ta) A.ncol() else A.nrow();
    const k = if (ta) A.nrow() else A.ncol();
    const n = if (tb) B.nrow() else B.ncol();

    const k2 = if (tb) B.ncol() else B.nrow();

    std.debug.assert(k == k2);

    std.debug.assert(C.nrow() == m);
    std.debug.assert(C.ncol() == n);

    gemmSlice(T, C.data, A.data, B.data, m, n, k, alpha, beta, ta, tb);
}

/// Performs general matrix-vector multiplication: y = alpha * op(A) * x + beta * y (GEMV).
pub fn gemv(comptime T: type, y: *Vector(T), A: Matrix(T), x: Vector(T), alpha: T, beta: T, ta: bool) void {
    if (ta) {
        std.debug.assert(A.nrow() == x.length());
        std.debug.assert(A.ncol() == y.length());
    }

    if (!ta) {
        std.debug.assert(A.ncol() == x.length());
        std.debug.assert(A.nrow() == y.length());
    }

    gemvSlice(T, y.data, A.data, x.data, A.nrow(), A.ncol(), alpha, beta, ta);
}

/// Computes the LU factorization of a square matrix with partial pivoting: P * A = L * U.
pub fn getrf(comptime T: type, A: *Matrix(T), ipiv: []i32) !void {
    std.debug.assert(A.nrow() == A.ncol());
    std.debug.assert(ipiv.len == A.nrow());

    try getrfSlice(T, A.data, ipiv);
}

/// Solves a system of linear equations A * X = B using precomputed LU factors from getrf.
pub fn getrs(comptime T: type, X: *Matrix(T), LU: Matrix(T), ipiv: []const i32, B: Matrix(T)) !void {
    std.debug.assert(LU.nrow() == LU.ncol());

    std.debug.assert(LU.nrow() == B.nrow());
    std.debug.assert(ipiv.len == LU.nrow());

    std.debug.assert(X.nrow() == B.nrow());
    std.debug.assert(X.ncol() == B.ncol());

    try getrsSlice(T, X.data, LU.data, ipiv, B.data);
}

/// Computes the Euclidean L2 norm of a vector using CBLAS dnrm2.
pub fn nrm2(comptime T: type, x: Vector(T)) T {
    return nrm2Slice(T, x.data);
}

/// Computes eigenvalues and eigenvectors of a real symmetric matrix via divide-and-conquer.
pub fn syevd(comptime T: type, W: *Vector(T), U: *Matrix(T), V: Matrix(T)) !void {
    std.debug.assert(V.ncol() == W.length());

    std.debug.assert(U.nrow() == V.nrow());
    std.debug.assert(U.ncol() == V.ncol());

    try syevdSlice(T, W.data, U.data, V.data);
}

/// Computes eigenvalues and eigenvectors for a batch of symmetric matrices.
pub fn syevdBatch(comptime T: type, W: *Matrix(T), U: *Matrix(T), V: Matrix(T)) !void {
    std.debug.assert(V.ncol() == W.ncol() * W.ncol());

    std.debug.assert(W.nrow() == V.nrow());
    std.debug.assert(U.nrow() == V.nrow());
    std.debug.assert(U.ncol() == V.ncol());

    for (0..V.nrow()) |i| {
        const Vi = V.rowSlice(i);
        const Wi = W.rowSlice(i);
        const Ui = U.rowSlice(i);

        try syevdSlice(T, Wi, Ui, Vi);
    }
}

/// Solves the generalized symmetric-definite eigenproblem: A * x = lambda * B * x.
pub fn sygvd(comptime T: type, W: *Vector(T), U: *Matrix(T), V: Matrix(T), B: *Matrix(T)) !void {
    std.debug.assert(V.ncol() == W.length());

    std.debug.assert(U.nrow() == V.nrow());
    std.debug.assert(U.ncol() == V.ncol());
    std.debug.assert(B.nrow() == V.nrow());
    std.debug.assert(B.ncol() == V.ncol());

    try sygvdSlice(T, W.data, U.data, V.data, B.data);
}

/// Performs the AXPY operation (y = alpha * x + y) on contiguous slices using CBLAS daxpy or zaxpy.
fn axpySlice(comptime T: type, y: []T, x: []const T, alpha: T) void {
    std.debug.assert(x.len == y.len);

    if (comptime primType(T) != f64) @compileError("AXPY ONLY SUPPORTS F64 PRIMITIVE TYPES");

    if (comptime T == f64) {
        return cblas.cblas_daxpy(@intCast(x.len), alpha, x.ptr, 1, y.ptr, 1);
    }

    if (comptime isComplex(T)) {
        return cblas.cblas_zaxpy(@intCast(x.len), &alpha, x.ptr, 1, y.ptr, 1);
    }
}

/// Computes the dot product of two contiguous slices using CBLAS ddot.
fn dotSlice(comptime T: type, x: []const T, y: []const T) T {
    std.debug.assert(x.len == y.len);

    if (comptime primType(T) != f64) @compileError("DOT ONLY SUPPORTS F64 PRIMITIVE TYPES");

    return cblas.cblas_ddot(@intCast(x.len), x.ptr, 1, y.ptr, 1);
}

/// Performs general matrix-matrix multiplication on raw contiguous slices using CBLAS dgemm.
fn gemmSlice(comptime T: type, C: []T, A: []const T, B: []const T, m: usize, n: usize, k: usize, alpha: T, beta: T, ta: bool, tb: bool) void {
    std.debug.assert(A.len == m * k);
    std.debug.assert(B.len == k * n);
    std.debug.assert(C.len == m * n);

    if (comptime primType(T) != f64) @compileError("GEMM ONLY SUPPORTS F64 NUMBERS");

    const at: c_uint = if (ta) @intCast(cblas.CblasTrans) else @intCast(cblas.CblasNoTrans);
    const bt: c_uint = if (tb) @intCast(cblas.CblasTrans) else @intCast(cblas.CblasNoTrans);

    const mi: i32 = @intCast(m);
    const ni: i32 = @intCast(n);
    const ki: i32 = @intCast(k);

    const Ap: [*]const f64 = @ptrCast(A.ptr);
    const Bp: [*]const f64 = @ptrCast(B.ptr);

    const Cp: [*]f64 = @ptrCast(C.ptr);

    const ldai: i32 = @intCast(if (ta) m else k);
    const ldbi: i32 = @intCast(if (tb) k else n);

    const ldci: i32 = @intCast(n);

    cblas.cblas_dgemm(cblas.CblasRowMajor, at, bt, mi, ni, ki, alpha, Ap, ldai, Bp, ldbi, beta, Cp, ldci);
}

/// Performs general matrix-vector multiplication on raw slices using CBLAS dgemv.
fn gemvSlice(comptime T: type, y: []T, A: []const T, x: []const T, m: usize, n: usize, alpha: T, beta: T, ta: bool) void {
    if (comptime primType(T) != f64) @compileError("GEMV ONLY SUPPORTS F64 NUMBERS");

    const trans: c_uint = if (ta) @intCast(cblas.CblasTrans) else @intCast(cblas.CblasNoTrans);

    const mi: i32 = @intCast(m);
    const ni: i32 = @intCast(n);

    cblas.cblas_dgemv(cblas.CblasRowMajor, trans, mi, ni, alpha, A.ptr, ni, x.ptr, 1, beta, y.ptr, 1);
}

/// Computes the LU factorization of a matrix slice using LAPACKE dgetrf.
fn getrfSlice(comptime T: type, A: []T, ipiv: []i32) !void {
    std.debug.assert(ipiv.len * ipiv.len == A.len);

    if (primType(T) != f64) @compileError("GETRF NOW ONLY SUPPORTS F64 NUMBERS");

    const n: i32 = @intCast(std.math.sqrt(A.len));

    const info = lapacke.LAPACKE_dgetrf(lapacke.LAPACK_ROW_MAJOR, n, n, A.ptr, n, ipiv.ptr);

    if (info != 0) return error.LapackError;
}

/// Solves a system of linear equations on slices using precomputed LU factors via LAPACKE dgetrs.
fn getrsSlice(comptime T: type, X: []T, LU: []const T, ipiv: []const i32, B: []const T) !void {
    std.debug.assert(ipiv.len * ipiv.len == LU.len);

    if (comptime primType(T) != f64) @compileError("GETRS NOW ONLY SUPPORTS F64 NUMBERS");

    const n: i32 = @intCast(std.math.sqrt(LU.len));
    const nofrhs: i32 = @intCast(B.len / ipiv.len);

    for (0..B.len) |i| {
        X[i] = B[i];
    }

    const info = lapacke.LAPACKE_dgetrs(ROW_MAJOR, 'N', n, nofrhs, LU.ptr, n, ipiv.ptr, X.ptr, nofrhs);

    if (info != 0) return error.LapackError;
}

/// Computes the Euclidean L2 norm of a slice using CBLAS dnrm2.
fn nrm2Slice(comptime T: type, x: []const T) primType(T) {
    if (comptime primType(T) != f64) @compileError("NRM2 ONLY SUPPORTS F64 PRIMITIVE TYPES");

    return cblas.cblas_dnrm2(@intCast(x.len), x.ptr, 1);
}

/// Solves the symmetric eigenvalue problem on contiguous slices using LAPACKE dsyevd.
fn syevdSlice(comptime T: type, W: []T, U: []T, V: []T) !void {
    std.debug.assert(V.len == W.len * W.len);
    std.debug.assert(U.len == W.len * W.len);

    if (primType(T) != f64) @compileError("SYEVD NOW ONLY SUPPORTS F64 NUMBERS");

    if (W.len == 1) {
        W[0], U[0] = .{ V[0], 1 };

        return;
    }

    if (W.len == 2) {
        const A, const B, const C = .{ V[0], 0.5 * (V[1] + V[2]), V[3] };

        const theta = 0.5 * std.math.atan2(2 * B, A - C);

        const c = @cos(theta);
        const s = @sin(theta);

        W[0] = A * s * s - 2 * B * c * s + C * c * c;
        W[1] = A * c * c + 2 * B * c * s + C * s * s;

        U[0], U[1], U[2], U[3] = .{ -s, c, c, s };

        return;
    }

    const n: i32 = @intCast(W.len);

    for (0..V.len) |i| {
        U[i] = V[i];
    }

    const info = lapacke.LAPACKE_dsyevd(ROW_MAJOR, 'V', 'U', n, U.ptr, n, W.ptr);

    if (info != 0) return error.LapackError;
}

/// Solves the generalized symmetric-definite eigenproblem on raw slices using LAPACKE dsygvd.
fn sygvdSlice(comptime T: type, W: []T, U: []T, V: []T, B: []T) !void {
    std.debug.assert(V.len == W.len * W.len);
    std.debug.assert(U.len == W.len * W.len);
    std.debug.assert(B.len == W.len * W.len);

    if (comptime primType(T) != f64) @compileError("SYGVD NOW ONLY SUPPORTS F64 NUMBERS");

    const n: i32 = @intCast(W.len);

    for (0..V.len) |i| {
        U[i] = V[i];
    }

    const info = lapacke.LAPACKE_dsygvd(lapacke.LAPACK_ROW_MAJOR, 1, 'V', 'U', n, U.ptr, n, B.ptr, n, W.ptr);

    if (info != 0) return error.LapackError;
}
