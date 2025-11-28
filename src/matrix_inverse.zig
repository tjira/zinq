//! File with linear solvers.

const std = @import("std");

const complex_matrix = @import("complex_matrix.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const errro_handling = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const Complex = std.math.complex.Complex;
const ComplexMatrix = complex_matrix.ComplexMatrix;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const eigensystemJacobiAlloc = eigenproblem_solver.eigensystemJacobiAlloc;
const eigensystemSymmetricAlloc = eigenproblem_solver.eigensystemSymmetricAlloc;
const mm = matrix_multiplication.mm;
const mmAlloc = matrix_multiplication.mmAlloc;
const throw = errro_handling.throw;

const SINGULARITY_TOLERANCE = global_variables.SINGULARITY_TOLERANCE;
const TEST_TOLERANCE = global_variables.TEST_TOLERANCE;

/// Form the inverse of a hermitian matrix A using its eigenvalue decomposition. Ainv = C * J_inv * C^T
pub fn inverseHermitian(comptime T: type, Ainv: *ComplexMatrix(T), A: ComplexMatrix(T), AJ: ComplexMatrix(T), AC: ComplexMatrix(T), Atmp: *ComplexMatrix(T)) !void {
    return try inverseSymmetricOrHermitian(ComplexMatrix, T, Ainv, A, AJ, AC, Atmp);
}

/// Form the inverse of a hermitian matrix A using its eigenvalue decomposition. Ainv = C * J_inv * C^T
pub fn inverseHermitianAlloc(comptime T: type, A: ComplexMatrix(T), allocator: std.mem.Allocator) !ComplexMatrix(T) {
    return try inverseSymmetricOrHermitianAlloc(ComplexMatrix, T, A, allocator);
}

/// Form the inverse of a symmetric matrix A using its eigenvalue decomposition. Ainv = C * J_inv * C^T
pub fn inverseSymmetric(comptime T: type, Ainv: *RealMatrix(T), A: RealMatrix(T), AJ: RealMatrix(T), AC: RealMatrix(T), Atmp: *RealMatrix(T)) !void {
    return try inverseSymmetricOrHermitian(RealMatrix, T, Ainv, A, AJ, AC, Atmp);
}

/// Form the inverse of a symmetric matrix A using its eigenvalue decomposition. Ainv = C * J_inv * C^T
pub fn inverseSymmetricAlloc(comptime T: type, A: RealMatrix(T), allocator: std.mem.Allocator) !RealMatrix(T) {
    return try inverseSymmetricOrHermitianAlloc(RealMatrix, T, A, allocator);
}

/// Form the inverse of a symmetric or hermitian matrix A using its eigenvalue decomposition. Ainv = C * J_inv * C^T
pub fn inverseSymmetricOrHermitian(comptime M: fn (comptime type) type, comptime T: type, Ainv: *M(T), A: M(T), AJ: M(T), AC: M(T), Atmp: *M(T)) !void {
    if (!AJ.isSquare() or !AC.isSquare()) return throw(void, "EIGENVALUE MATRIX OR EIGENVECTOR MATRIX IS NOT SQUARE AND THE INVERSE CANNOT BE FORMED", .{});
    if (AJ.rows != AC.rows or AJ.cols != AC.cols) return throw(void, "EIGENVALUE MATRIX AND EIGENVECTOR MATRIX MUST HAVE THE SAME DIMENSIONS", .{});

    if (comptime M == RealMatrix) {
        if (!A.isSymmetric(0)) return throw(void, "THE MATRIX YOU ARE PASSING TO THE INVERSE SYMMETRIC FUNCTION IS NOT SYMMETRIC", .{});
        for (0..AJ.rows) |i| if (@abs(AJ.at(i, i)) < SINGULARITY_TOLERANCE) return throw(void, "THE MATRIX IS SINGULAR AND THE INVERSE CANNOT BE FORMED", .{});
    } else {
        if (!A.isHermitian(0)) return throw(void, "THE MATRIX YOU ARE PASSING TO THE JACOBI EIGENSOLVER IS NOT HERMITIAN", .{});
        for (0..AJ.rows) |i| if (std.math.complex.abs(AJ.at(i, i)) < SINGULARITY_TOLERANCE) return throw(void, "THE MATRIX IS SINGULAR AND THE INVERSE CANNOT BE FORMED", .{});
    }

    try AJ.copyTo(Ainv);

    for (0..Ainv.rows) |i| Ainv.ptr(i, i).* = if (comptime M == RealMatrix) 1 / Ainv.at(i, i) else Complex(T).init(1, 0).div(Ainv.at(i, i));
    
    try mm(T, Atmp, AC, false, Ainv.*, false); try mm(T, Ainv, Atmp.*, false, AC, true);
}

/// Form the inverse of a symmetric or hermitian matrix A using its eigenvalue decomposition. Ainv = C * J_inv * C^T
pub fn inverseSymmetricOrHermitianAlloc(comptime M: fn (comptime type) type, comptime T: type, A: M(T), allocator: std.mem.Allocator) !M(T) {
    const AJC = try eigensystemJacobiAlloc(M, T, A, allocator); defer AJC.J.deinit(); defer AJC.C.deinit();

    var Ainv = try M(T).init(A.rows, A.cols, allocator);
    var Atmp = try M(T).init(A.rows, A.cols, allocator);

    try inverseSymmetricOrHermitian(M, T, &Ainv, A, AJC.J, AJC.C, &Atmp);

    Atmp.deinit();

    return Ainv;
}

test "Symmetric 3x3 Inverse" {
    var A = try RealMatrix(f64).init(3, 3, std.testing.allocator); defer A.deinit();

    var Ainv_expected = try RealMatrix(f64).init(3, 3, std.testing.allocator); defer Ainv_expected.deinit();

    A.ptr(0, 0).* = 1.0; A.ptr(0, 1).* = 2.0; A.ptr(0, 2).* = 3.0;
    A.ptr(1, 0).* = 2.0; A.ptr(1, 1).* = 4.0; A.ptr(1, 2).* = 5.0;
    A.ptr(2, 0).* = 3.0; A.ptr(2, 1).* = 5.0; A.ptr(2, 2).* = 6.0;

    Ainv_expected.ptr(0, 0).* =  1.0; Ainv_expected.ptr(0, 1).* = -3.0; Ainv_expected.ptr(0, 2).* =  2.0;
    Ainv_expected.ptr(1, 0).* = -3.0; Ainv_expected.ptr(1, 1).* =  3.0; Ainv_expected.ptr(1, 2).* = -1.0;
    Ainv_expected.ptr(2, 0).* =  2.0; Ainv_expected.ptr(2, 1).* = -1.0; Ainv_expected.ptr(2, 2).* =  0.0;

    const Ainv = try inverseSymmetricAlloc(f64, A, std.testing.allocator); defer Ainv.deinit();

    try std.testing.expect(Ainv.eq(Ainv_expected, TEST_TOLERANCE));
}
