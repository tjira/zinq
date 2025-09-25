//! File with linear solvers.

const std = @import("std");

const eigenproblem_solver = @import("eigenproblem_solver.zig");
const errro_handling = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const eigensystemSymmetricAlloc = eigenproblem_solver.eigensystemSymmetricAlloc;
const mm = matrix_multiplication.mm;
const mmAlloc = matrix_multiplication.mmAlloc;
const throw = errro_handling.throw;

const SINGULARITY_TOLERANCE = global_variables.SINGULARITY_TOLERANCE;
const TEST_TOLERANCE = global_variables.TEST_TOLERANCE;

/// Solve the linear system Ax = b using the eigenvalue decomposition of A. The matrix A must be symmetric.
pub fn linearSolveSymmetric(comptime T: type, x: *RealVector(T), A: RealMatrix(T), AJ: RealMatrix(T), AC: RealMatrix(T), b: RealVector(T), y: *RealVector(T)) !void {
    if (!A.isSymmetric(0)) return throw(void, "THE MATRIX YOU ARE PASSING TO THE SYMMETRIC LINEAR SYSTEM SOLVER IS NOT SYMMETRIC", .{});
    if (!AJ.isSquare() or !AC.isSquare()) return throw(void, "EIGENVALUE MATRIX OR EIGENVECTOR MATRIX IS NOT SQUARE AND THE LINEAR SYSTEM CAN'T BE SOLVED", .{});
    if (AJ.rows != AC.rows or AJ.cols != AC.cols) return throw(void, "EIGENVALUE MATRIX AND EIGENVECTOR MATRIX MUST HAVE THE SAME DIMENSIONS", .{});
    if (AJ.rows != b.len) return throw(void, "THE LENGTH OF THE RIGHT-HAND SIDE VECTOR MUST BE EQUAL TO THE DIMENSIONS OF THE EIGENVALUE MATRIX", .{});
    if (x.len != b.len) return throw(void, "THE LENGTH OF THE SOLUTION VECTOR MUST BE EQUAL TO THE LENGTH OF THE RIGHT-HAND SIDE VECTOR", .{});
    if (y.len != b.len) return throw(void, "THE LENGTH OF THE TEMPORARY VECTOR MUST BE EQUAL TO THE LENGTH OF THE RIGHT-HAND SIDE VECTOR", .{});
    for (0..AJ.rows) |i| if (@abs(AJ.at(i, i)) < SINGULARITY_TOLERANCE) return throw(void, "THE MATRIX IS SINGULAR AND THE LINEAR SYSTEM CAN'T BE SOLVED", .{});

    var x_matrix = x.asMatrix();
    var y_matrix = y.asMatrix();
    
    try mm(T, &y_matrix, AC, true, b.asMatrix(), false);

    for (0..y.len) |i| y.ptr(i).* /= AJ.at(i, i);

    try mm(T, &x_matrix, AC, false, y_matrix, false);
}

/// Solve the linear system Ax = b using the eigenvalue decomposition of A. the result is returned as a new vector.
pub fn linearSolveSymmetricAlloc(comptime T: type, A: RealMatrix(T), b: RealVector(T), allocator: std.mem.Allocator) !RealVector(T) {
    const AJC = try eigensystemSymmetricAlloc(T, A, allocator); defer AJC.J.deinit(); defer AJC.C.deinit();

    var x = try RealVector(T).init(b.len, allocator);
    var y = try RealVector(T).init(b.len, allocator);

    try linearSolveSymmetric(T, &x, A, AJC.J, AJC.C, b, &y);

    y.deinit();

    return x;
}

test "Symmetric 3x3 Linear System" {
    var A = try real_matrix.RealMatrix(f64).init(3, 3, std.testing.allocator); defer A.deinit();
    var b = try real_vector.RealVector(f64).init(3, std.testing.allocator); defer b.deinit();

    var x_expected = try real_vector.RealVector(f64).init(3, std.testing.allocator); defer x_expected.deinit();

    A.ptr(0, 0).* = 4; A.ptr(0, 1).* = 1; A.ptr(0, 2).* = 2;
    A.ptr(1, 0).* = 1; A.ptr(1, 1).* = 3; A.ptr(1, 2).* = 0;
    A.ptr(2, 0).* = 2; A.ptr(2, 1).* = 0; A.ptr(2, 2).* = 5;

    b.ptr(0).* = 4; b.ptr(1).* = 5; b.ptr(2).* = 6;

    x_expected.ptr(0).* = -0.02325581395349; x_expected.ptr(1).* = 1.67441860465116; x_expected.ptr(2).* = 1.20930232558139;

    const x = try linearSolveSymmetricAlloc(f64, A, b, std.testing.allocator); defer x.deinit();

    try std.testing.expect(x.eq(x_expected, TEST_TOLERANCE));
}
