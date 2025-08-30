//! Collection of more advanced linear algebra operations.

const std = @import("std");

const global_variables = @import("global_variables.zig");
const real_matrix = @import("real_matrix.zig");
const complex_matrix = @import("complex_matrix.zig");

const Complex = std.math.complex.Complex;
const RealMatrix = real_matrix.RealMatrix;
const ComplexMatrix = complex_matrix.ComplexMatrix;

const TEST_TOLERANCE = global_variables.TEST_TOLERANCE;

/// Diagonalize a real symmetric matrix A. The provided matrix is overwritten by the diagonal form.
pub fn diagonalizeSymmetric(comptime T: type, A: *RealMatrix(T)) !void {
    if (A.rows == 1) return;

    if (A.rows == 2) {
        diagonalizeSymmetric2x2(T, A);
    }

    else return error.NotImplemented;
}

/// Diagonalize a real symmetric 2x2 matrix A using an analytical formula. The provided matrix is overwritten by the diagonal form.
pub fn diagonalizeSymmetric2x2(comptime T: type, A: *RealMatrix(T)) void {
    const a = A.at(0, 0); const b = A.at(0, 1); const c = A.at(1, 1);

    const sqrt = std.math.sqrt(a * a + 4 * b * b - 2 * a * c + c * c);

    A.ptr(0, 0).* = 0.5 * (a + c - sqrt);
    A.ptr(1, 1).* = 0.5 * (a + c + sqrt);

    A.ptr(0, 1).* = 0;
    A.ptr(1, 0).* = 0;
}

/// Solve the eigenproblem for a real symmetric system.
pub fn eigensystemSymmetric(comptime T: type, A_eigenvalues: *RealMatrix(T), A_eigenvectors: *RealMatrix(T), A: RealMatrix(T)) !void {
    if (A.rows == 1) {
        @memcpy(A_eigenvalues.data, A.data); A_eigenvectors.fill(1);
    }

    else if (A.rows == 2) {
        eigensystemSymmetric2x2(T, A_eigenvalues, A_eigenvectors, A);
    }

    else return error.NotImplemented;
}

/// Eigenproblem for a real symmetric 2x2 system using an analytical formula.
pub fn eigensystemSymmetric2x2(comptime T: type, A_eigenvalues: *RealMatrix(T), A_eigenvectors: *RealMatrix(T), A: RealMatrix(T)) void {
    const a = A.at(0, 0); const b = A.at(0, 1); const c = A.at(1, 1); const d = std.math.hypot(a - c, 2 * b);

    A_eigenvalues.ptr(0, 0).* = 0.5 * (a + c - d);
    A_eigenvalues.ptr(1, 1).* = 0.5 * (a + c + d);

    A_eigenvalues.ptr(0, 1).* = 0;
    A_eigenvalues.ptr(1, 0).* = 0;

    A_eigenvectors.ptr(0, 0).* = b; A_eigenvectors.ptr(1, 0).* = A_eigenvalues.at(0, 0) - a;

    var norm = std.math.hypot(A_eigenvectors.at(0, 0), A_eigenvectors.at(1, 0));

    if (norm == 0) {

        A_eigenvectors.ptr(0, 0).* = c - A_eigenvalues.at(0, 0); A_eigenvectors.ptr(1, 0).* = -b;

        norm = std.math.hypot(A_eigenvectors.at(0, 0), A_eigenvectors.at(1, 0));
    }

    A_eigenvectors.ptr(0, 0).* /= norm;
    A_eigenvectors.ptr(1, 0).* /= norm;

    A_eigenvectors.ptr(0, 1).* = -A_eigenvectors.at(1, 0);
    A_eigenvectors.ptr(1, 1).* =  A_eigenvectors.at(0, 0);
}

/// Multiply two complex matrices A and B and store the result in C.
pub fn mmComplex(comptime T: type, C: *ComplexMatrix(T), A: ComplexMatrix(T), B: ComplexMatrix(T)) void {
    for (0..A.rows) |i| for (0..B.cols) |j| {

        var sum = Complex(T).init(0, 0);

        for (0..A.cols) |k| {
            sum = sum.add(A.at(i, k).mul(B.at(k, j)));
        }

        C.ptr(i, j).* = sum;
    };
}

/// Multiply two complex matrices A and B and return the result. The function returns an error if the allocation fails.
pub fn mmComplexAlloc(comptime T: type, A: ComplexMatrix(T), B: ComplexMatrix(T)) !ComplexMatrix(T) {
    var C = try ComplexMatrix(T).init(A.rows, B.cols, A.allocator);

    mmComplex(T, &C, A, B);

    return C;
}

/// Multiply two matrices A and B and store the result in C.
pub fn mmReal(comptime T: type, C: *RealMatrix(T), A: RealMatrix(T), B: RealMatrix(T)) void {
    for (0..A.rows) |i| for (0..B.cols) |j| {

        var sum: T = 0;

        for (0..A.cols) |k| {
            sum += A.at(i, k) * B.at(k, j);
        }

        C.ptr(i, j).* = sum;
    };
}

/// Multiply two matrices A and B and return the result. The function returns an error if the allocation fails.
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

    try std.testing.expect(C.eq(D, TEST_TOLERANCE));
}
