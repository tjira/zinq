const std = @import("std");

const lapacke = @import("cimport.zig").lapacke;
const cblas = @import("cimport.zig").cblas;

const Matrix = @import("tensor.zig").Matrix;
const Vector = @import("tensor.zig").Vector;

const primType = @import("value.zig").primType;

pub fn eighBatch(comptime T: type, W: *Matrix(T), U: *Matrix(T), V: Matrix(T)) !void {
    std.debug.assert(V.ncol() == W.ncol() * W.ncol());

    std.debug.assert(W.nrow() == V.nrow());
    std.debug.assert(U.nrow() == V.nrow());
    std.debug.assert(U.ncol() == V.ncol());

    for (0..V.nrow()) |i| {
        const Vi = V.rowSlice(i);
        const Wi = W.rowSlice(i);
        const Ui = U.rowSlice(i);

        try eighSlice(T, Wi, Ui, Vi);
    }
}

pub fn eighSlice(comptime T: type, W: []T, U: []T, V: []T) !void {
    std.debug.assert(V.len == W.len * W.len);
    std.debug.assert(U.len == W.len * W.len);

    if (primType(T) != f64) @compileError("EIGH NOW ONLY SUPPORTS F64 NUMBERS");

    const n: i32 = @intCast(W.len);

    for (0..V.len) |i| {
        U[i] = V[i];
    }

    const info = lapacke.LAPACKE_dsyevd(lapacke.LAPACK_ROW_MAJOR, 'V', 'U', n, U.ptr, n, W.ptr);

    if (info != 0) return error.LapackError;
}

pub fn geigh(comptime T: type, W: *Vector(T), U: *Matrix(T), V: Matrix(T), B: *Matrix(T)) !void {
    std.debug.assert(V.ncol() == W.length());

    std.debug.assert(U.nrow() == V.nrow());
    std.debug.assert(U.ncol() == V.ncol());
    std.debug.assert(B.nrow() == V.nrow());
    std.debug.assert(B.ncol() == V.ncol());

    try geighSlice(T, W.data, U.data, V.data, B.data);
}

pub fn luFactorize(comptime T: type, A: *Matrix(T), ipiv: []i32) !void {
    std.debug.assert(A.nrow() == A.ncol());
    std.debug.assert(ipiv.len == A.nrow());

    try luFactorizeSlice(T, A.data, ipiv);
}

pub fn luSolve(comptime T: type, X: *Matrix(T), LU: Matrix(T), ipiv: []const i32, B: Matrix(T)) !void {
    std.debug.assert(LU.nrow() == LU.ncol());

    std.debug.assert(LU.nrow() == B.nrow());
    std.debug.assert(ipiv.len == LU.nrow());

    std.debug.assert(X.nrow() == B.nrow());
    std.debug.assert(X.ncol() == B.ncol());

    try luSolveSlice(T, X.data, LU.data, ipiv, B.data);
}

pub fn mm(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T)) void {
    std.debug.assert(A.ncol() == B.nrow());
    std.debug.assert(C.nrow() == A.nrow());
    std.debug.assert(C.ncol() == B.ncol());

    mmSlice(T, C.data, A.data, B.data, A.nrow(), B.ncol(), A.ncol());
}

pub fn mmSlice(comptime T: type, C: []T, A: []const T, B: []const T, m: usize, n: usize, k: usize) void {
    std.debug.assert(A.len == m * k);
    std.debug.assert(B.len == k * n);
    std.debug.assert(C.len == m * n);

    if (primType(T) != f64) @compileError("MM ONLY SUPPORTS F64 NUMBERS");

    const at = cblas.CblasNoTrans;
    const bt = cblas.CblasNoTrans;

    const mi: i32 = @intCast(m);
    const ni: i32 = @intCast(n);
    const ki: i32 = @intCast(k);

    const Ap: [*]const f64 = @ptrCast(A.ptr);
    const Bp: [*]const f64 = @ptrCast(B.ptr);

    const Cp: [*]f64 = @ptrCast(C.ptr);

    cblas.cblas_dgemm(cblas.CblasRowMajor, at, bt, mi, ni, ki, 1, Ap, ki, Bp, ni, 0, Cp, ni);
}

fn geighSlice(comptime T: type, W: []T, U: []T, V: []T, B: []T) !void {
    std.debug.assert(V.len == W.len * W.len);
    std.debug.assert(U.len == W.len * W.len);
    std.debug.assert(B.len == W.len * W.len);

    if (primType(T) != f64) @compileError("GEIGH NOW ONLY SUPPORTS F64 NUMBERS");

    const n: i32 = @intCast(W.len);

    for (0..V.len) |i| {
        U[i] = V[i];
    }

    const info = lapacke.LAPACKE_dsygvd(lapacke.LAPACK_ROW_MAJOR, 1, 'V', 'U', n, U.ptr, n, B.ptr, n, W.ptr);

    if (info != 0) return error.LapackError;
}

fn luFactorizeSlice(comptime T: type, A: []T, ipiv: []i32) !void {
    std.debug.assert(ipiv.len * ipiv.len == A.len);

    if (primType(T) != f64) @compileError("LU FACTORIZE NOW ONLY SUPPORTS F64 NUMBERS");

    const n: i32 = @intCast(std.math.sqrt(A.len));

    const info = lapacke.LAPACKE_dgetrf(lapacke.LAPACK_ROW_MAJOR, n, n, A.ptr, n, ipiv.ptr);

    if (info != 0) return error.LapackError;
}

fn luSolveSlice(comptime T: type, X: []T, LU: []const T, ipiv: []const i32, B: []const T) !void {
    std.debug.assert(ipiv.len * ipiv.len == LU.len);

    if (primType(T) != f64) @compileError("LU SOLVE NOW ONLY SUPPORTS F64 NUMBERS");

    const n: i32 = @intCast(std.math.sqrt(LU.len));
    const nofrhs: i32 = @intCast(B.len / ipiv.len);

    for (0..B.len) |i| {
        X[i] = B[i];
    }

    const info = lapacke.LAPACKE_dgetrs(lapacke.LAPACK_ROW_MAJOR, 'N', n, nofrhs, LU.ptr, n, ipiv.ptr, X.ptr, nofrhs);

    if (info != 0) return error.LapackError;
}
