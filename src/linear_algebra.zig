const std = @import("std");

const lapacke = @cImport(@cInclude("lapacke.h"));

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

pub fn geighSlice(comptime T: type, W: []T, U: []T, V: []T, B: []T) !void {
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

pub fn luFactorize(comptime T: type, A: *Matrix(T), ipiv: []i32) !void {
    std.debug.assert(A.nrow() == A.ncol());
    std.debug.assert(ipiv.len == A.nrow());

    try luFactorizeSlice(T, A.data, ipiv);
}

pub fn luFactorizeSlice(comptime T: type, A: []T, ipiv: []i32) !void {
    std.debug.assert(ipiv.len * ipiv.len == A.len);

    if (primType(T) != f64) @compileError("LU FACTORIZE NOW ONLY SUPPORTS F64 NUMBERS");

    const n: i32 = @intCast(std.math.sqrt(A.len));

    const info = lapacke.LAPACKE_dgetrf(lapacke.LAPACK_ROW_MAJOR, n, n, A.ptr, n, ipiv.ptr);

    if (info != 0) return error.LapackError;
}

pub fn luSolve(comptime T: type, X: *Matrix(T), LU: Matrix(T), ipiv: []const i32, B: Matrix(T)) !void {
    std.debug.assert(LU.nrow() == LU.ncol());

    std.debug.assert(LU.nrow() == B.nrow());
    std.debug.assert(ipiv.len == LU.nrow());

    std.debug.assert(X.nrow() == B.nrow());
    std.debug.assert(X.ncol() == B.ncol());

    try luSolveSlice(T, X.data, LU.data, ipiv, B.data);
}

pub fn luSolveSlice(comptime T: type, X: []T, LU: []const T, ipiv: []const i32, B: []const T) !void {
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
