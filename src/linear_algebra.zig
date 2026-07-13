const std = @import("std");

const lapacke = @import("cimport.zig").lapacke;
const cblas = @import("cimport.zig").cblas;

const Matrix = @import("tensor.zig").Matrix;
const Vector = @import("tensor.zig").Vector;

const primType = @import("value.zig").primType;

pub fn addScaled(comptime T: type, alpha: T, x: Vector(T), y: *Vector(T)) void {
    std.debug.assert(x.length() == y.length());

    addScaledSlice(T, x.length(), alpha, x.data, y.data);
}

pub fn dot(comptime T: type, x: Vector(T), y: Vector(T)) T {
    std.debug.assert(x.length() == y.length());

    return dotSlice(T, x.length(), x.data, y.data);
}

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

pub fn eigh(comptime T: type, W: *Vector(T), U: *Matrix(T), V: Matrix(T)) !void {
    std.debug.assert(V.ncol() == W.length());

    std.debug.assert(U.nrow() == V.nrow());
    std.debug.assert(U.ncol() == V.ncol());

    try eighSlice(T, W.data, U.data, V.data);
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

pub fn mm(comptime T: type, C: *Matrix(T), A: Matrix(T), B: Matrix(T), alpha: T, beta: T, ta: bool, tb: bool) void {
    const m = if (ta) A.ncol() else A.nrow();
    const k = if (ta) A.nrow() else A.ncol();
    const n = if (tb) B.nrow() else B.ncol();

    const k2 = if (tb) B.ncol() else B.nrow();

    std.debug.assert(k == k2);

    std.debug.assert(C.nrow() == m);
    std.debug.assert(C.ncol() == n);

    mmSlice(T, C.data, A.data, B.data, m, n, k, alpha, beta, ta, tb);
}

pub fn mmSlice(comptime T: type, C: []T, A: []const T, B: []const T, m: usize, n: usize, k: usize, alpha: T, beta: T, ta: bool, tb: bool) void {
    std.debug.assert(A.len == m * k);
    std.debug.assert(B.len == k * n);
    std.debug.assert(C.len == m * n);

    if (comptime primType(T) != f64) @compileError("MM ONLY SUPPORTS F64 NUMBERS");

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

pub fn mmv(comptime T: type, y: *Vector(T), A: Matrix(T), x: Vector(T), a: T, b: T, tr: bool) void {
    if (tr) {
        std.debug.assert(A.nrow() == x.length());
        std.debug.assert(A.ncol() == y.length());
    }

    if (!tr) {
        std.debug.assert(A.ncol() == x.length());
        std.debug.assert(A.nrow() == y.length());
    }

    mmvSlice(T, y.data, A.data, x.data, A.nrow(), A.ncol(), a, b, tr, 1);
}

pub fn norm(comptime T: type, x: Vector(T)) T {
    return normSlice(T, x.length(), x.data);
}

fn addScaledSlice(comptime T: type, n: usize, alpha: T, x: []const T, y: []T) void {
    if (comptime primType(T) != f64) @compileError("ADDSCALED ONLY SUPPORTS F64 PRIMITIVE TYPES");

    cblas.cblas_daxpy(@intCast(n), alpha, x.ptr, 1, y.ptr, 1);
}

fn dotSlice(comptime T: type, n: usize, x: []const T, y: []const T) T {
    if (comptime primType(T) != f64) @compileError("DOT ONLY SUPPORTS F64 PRIMITIVE TYPES");

    return cblas.cblas_ddot(@intCast(n), x.ptr, 1, y.ptr, 1);
}

fn geighSlice(comptime T: type, W: []T, U: []T, V: []T, B: []T) !void {
    std.debug.assert(V.len == W.len * W.len);
    std.debug.assert(U.len == W.len * W.len);
    std.debug.assert(B.len == W.len * W.len);

    if (comptime primType(T) != f64) @compileError("GEIGH NOW ONLY SUPPORTS F64 NUMBERS");

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

    if (comptime primType(T) != f64) @compileError("LU SOLVE NOW ONLY SUPPORTS F64 NUMBERS");

    const n: i32 = @intCast(std.math.sqrt(LU.len));
    const nofrhs: i32 = @intCast(B.len / ipiv.len);

    for (0..B.len) |i| {
        X[i] = B[i];
    }

    const info = lapacke.LAPACKE_dgetrs(lapacke.LAPACK_ROW_MAJOR, 'N', n, nofrhs, LU.ptr, n, ipiv.ptr, X.ptr, nofrhs);

    if (info != 0) return error.LapackError;
}

fn mmvSlice(comptime T: type, y: []T, A: []const T, x: []const T, m: usize, n: usize, a: T, b: T, tr: bool, sx: i32) void {
    if (comptime primType(T) != f64) @compileError("MMV ONLY SUPPORTS F64 NUMBERS");

    const trans: c_uint = if (tr) @intCast(cblas.CblasTrans) else @intCast(cblas.CblasNoTrans);

    const mi: i32 = @intCast(m);
    const ni: i32 = @intCast(n);

    cblas.cblas_dgemv(cblas.CblasRowMajor, trans, mi, ni, a, A.ptr, ni, x.ptr, sx, b, y.ptr, 1);
}

fn normSlice(comptime T: type, n: usize, x: []const T) primType(T) {
    if (comptime primType(T) != f64) @compileError("NORM ONLY SUPPORTS F64 PRIMITIVE TYPES");

    return cblas.cblas_dnrm2(@intCast(n), x.ptr, 1);
}
