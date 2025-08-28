//! Collection of more advanced linear algebra operations.

const std = @import("std");

const real_matrix = @import("real_matrix.zig");

const RealMatrix = real_matrix.RealMatrix;

/// Solve the eigenproblem for a real symmetric system.
pub fn eighReal(comptime T: type, A_eigenvalues: *RealMatrix(T), A_eigenvectors: *RealMatrix(T), A: RealMatrix(T)) !void {
    if (A.rows == 1) {
        @memcpy(A_eigenvalues.data, A.data); A_eigenvectors.fill(1);
    }

    else if (A.rows == 2) {
        eighReal2x2(T, A_eigenvalues, A_eigenvectors, A);
    }

    else return error.NotImplemented;
}

/// Eigenproblem for a real symmetric 2x2 system using an analytical formula.
pub fn eighReal2x2(comptime T: type, A_eigenvalues: *RealMatrix(T), A_eigenvectors: *RealMatrix(T), A: RealMatrix(T)) void {
    const a = A.at(0, 0); const b = A.at(0, 1); const c = A.at(1, 1);

    A_eigenvalues.fill(0);

    const sqrt = std.math.sqrt(a * a + 4 * b * b - 2 * a * c + c * c);

    A_eigenvalues.ptr(0, 0).* = 0.5 * (a + c - sqrt);
    A_eigenvalues.ptr(1, 1).* = 0.5 * (a + c + sqrt);

    A_eigenvectors.ptr(0, 0).* = 0.5 * (a - c - sqrt) / b;
    A_eigenvectors.ptr(0, 1).* = 0.5 * (a - c + sqrt) / b;

    const norm1 = std.math.sqrt(A_eigenvectors.at(0, 0) * A_eigenvectors.at(0, 0) + 1);
    const norm2 = std.math.sqrt(A_eigenvectors.at(0, 1) * A_eigenvectors.at(0, 1) + 1);

    A_eigenvectors.ptr(0, 0).* /= norm1;
    A_eigenvectors.ptr(0, 1).* /= norm2;

    A_eigenvectors.ptr(1, 0).* = 1 / norm1;
    A_eigenvectors.ptr(1, 1).* = 1 / norm2;
}

/// Multiply two matrices A and B and store the result in C. The function returns an error if the dimensions are incompatible.
pub fn mmReal(comptime T: type, C: *RealMatrix(T), A: RealMatrix(T), B: RealMatrix(T)) void {
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

    mmReal(T, &C, A, B);

    return C;
}

test "mmRealAlloc" {
    var A = try RealMatrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try RealMatrix(f64).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try RealMatrix(f64).init(2, 2, std.testing.allocator); defer C.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2; A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;
    B.ptr(0, 0).* = 5; B.ptr(0, 1).* = 6; B.ptr(1, 0).* = 7; B.ptr(1, 1).* = 8;

    C.ptr(0, 0).* = 19; C.ptr(0, 1).* = 22; C.ptr(1, 0).* = 43; C.ptr(1, 1).* = 50;

    var D = try mmRealAlloc(f64, A, B); defer D.deinit();

    try std.testing.expect(C.eq(D));
}
