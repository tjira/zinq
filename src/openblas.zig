const std = @import("std");

const Matrix = @import("tensor.zig").Matrix;

pub fn eigh2x2(comptime T: type, W: []T, U: []T, V: []T) !void {
    const a = V[0];
    const b = V[1];
    const c = V[3];

    const R = std.math.sqrt(4 * b * b + (a - c) * (a - c));

    W[0] = 0.5 * (a + c - R);
    W[1] = 0.5 * (a + c + R);

    const theta = 0.5 * std.math.atan2(2 * b, a - c);

    const cos_t = @cos(theta);
    const sin_t = @sin(theta);

    // zig fmt: off
    U[0] =  sin_t;
    U[1] = -cos_t;
    U[2] =  cos_t;
    U[3] =  sin_t;
    // zig fmt: on
}

pub fn eighBatch(comptime T: type, W: *Matrix(T), U: *Matrix(T), V: Matrix(T)) !void {
    if (V.nrow() == 1) return eighBatch1x1(T, W, U, V);
    if (V.nrow() == 4) return eighBatch2x2(T, W, U, V);
}

pub fn eighBatch1x1(comptime T: type, W: *Matrix(T), U: *Matrix(T), V: Matrix(T)) !void {
    for (0..V.ncol()) |j| {
        W.ptr(0, j).* = V.at(0, j);
    }

    U.fill(1);
}

pub fn eighBatch2x2(comptime T: type, W: *Matrix(T), U: *Matrix(T), V: Matrix(T)) !void {
    for (0..V.ncol()) |j| {
        const Vj = V.colSlice(j);
        const Wj = W.colSlice(j);
        const Uj = U.colSlice(j);

        try eigh2x2(T, Wj, Uj, Vj);
    }
}
