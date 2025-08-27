//! Collection of more advanced linear algebra operations.

const std = @import("std");

const real_matrix = @import("real_matrix.zig");

const RealMatrix = real_matrix.RealMatrix;

/// Multiply two matrices A and B and store the result in C. The function returns an error if the dimensions are incompatible.
pub fn mmReal(comptime T: type, C: *RealMatrix(T), A: RealMatrix(T), B: RealMatrix(T)) !void {
    if (A.cols != B.rows) return error.IncompatibleDimensions;

    for (0..A.rows) |i| for (0..B.cols) |j| {
        var sum: T = 0;

        for (0..A.cols) |k| {
            sum += A.at(i, k) * B.at(k, j);
        }

        C.ptr(i, j).* = sum;
    };
}

/// Multiply two matrices A and B and return the result. The function returns an error if the allocation fails or if the dimensions are incompatible.
pub fn mmRealAlloc(comptime T: type, A: RealMatrix(T), B: RealMatrix(T)) !RealMatrix(T) {
    var C = try RealMatrix(T).init(A.rows, B.cols, A.allocator);

    try mmReal(T, &C, A, B);

    return C;
}

test "mmReal, mmRealAlloc" {
    var A = try RealMatrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try RealMatrix(f64).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try RealMatrix(f64).init(2, 2, std.testing.allocator); defer C.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2; A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;
    B.ptr(0, 0).* = 5; B.ptr(0, 1).* = 6; B.ptr(1, 0).* = 7; B.ptr(1, 1).* = 8;

    C.ptr(0, 0).* = 19; C.ptr(0, 1).* = 22; C.ptr(1, 0).* = 43; C.ptr(1, 1).* = 50;

    var D = try mmRealAlloc(f64, A, B); defer D.deinit();

    try std.testing.expect(C.eq(D));
}
