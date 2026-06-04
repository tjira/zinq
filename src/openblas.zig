const std = @import("std");

const lapacke = @cImport(@cInclude("lapacke.h"));

const Matrix = @import("tensor.zig").Matrix;

pub fn eighBatch(comptime T: type, W: *Matrix(T), U: *Matrix(T), V: Matrix(T)) !void {
    for (0..V.nrow()) |i| {
        const Vi = V.rowSlice(i);
        const Wi = W.rowSlice(i);
        const Ui = U.rowSlice(i);

        try eighSlice(T, Wi, Ui, Vi);
    }
}

pub fn eighSlice(comptime T: type, W: []T, U: []T, V: []T) !void {
    const n: i32 = @intCast(W.len);

    for (0..V.len) |i| {
        U[i] = V[i];
    }

    for (0..W.len) |i| {
        W[i] = V[i];
    }

    const info = lapacke.LAPACKE_dsyevd(lapacke.LAPACK_ROW_MAJOR, 'V', 'U', n, U.ptr, n, W.ptr);

    if (info != 0) return error.LapackError;
}
