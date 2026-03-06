//! File with linear solvers.

const std = @import("std");

const eigenproblem_solver = @import("eigenproblem_solver.zig");
const global_variables = @import("global_variables.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const eigensystemHermitianAlloc = eigenproblem_solver.eigensystemHermitianAlloc;
const mm = matrix_multiplication.mm;
const mmAlloc = matrix_multiplication.mmAlloc;

const TEST_TOLERANCE = global_variables.TEST_TOLERANCE;

/// Solve the linear system Ax = b using the eigenvalue decomposition of A. The matrix A must be symmetric.
pub fn linearSolveSymmetric(comptime T: type, x: *RealVector(T), AJ: RealMatrix(T), AC: RealMatrix(T), b: RealVector(T), y: *RealVector(T)) !void {
    const tolerance = @as(T, @floatFromInt(AJ.rows)) * std.math.floatEps(T) * AJ.maxAbsDiagonal();

    for (0..AJ.rows) |i| if (@abs(AJ.at(i, i)) <= tolerance) return error.SingularMatrix;

    var x_matrix = x.asMatrix();
    var y_matrix = y.asMatrix();
    
    try mm(T, &y_matrix, AC, true, b.asMatrix(), false);

    for (0..y.len) |i| y.ptr(i).* /= AJ.at(i, i);

    try mm(T, &x_matrix, AC, false, y_matrix, false);
}

/// Solve the linear system Ax = b using the eigenvalue decomposition of A. the result is returned as a new vector.
pub fn linearSolveSymmetricAlloc(comptime T: type, A: RealMatrix(T), b: RealVector(T), allocator: std.mem.Allocator) !RealVector(T) {
    const AJC = try eigensystemHermitianAlloc(T, A, allocator); defer AJC.J.deinit(allocator); defer AJC.C.deinit(allocator);

    var x = try RealVector(T).init(b.len, allocator);
    var y = try RealVector(T).init(b.len, allocator);

    try linearSolveSymmetric(T, &x, AJC.J, AJC.C, b, &y);

    y.deinit(allocator);

    return x;
}

test "Symmetric 3x3 Linear System" {
    var A = try real_matrix.RealMatrix(f64).init(3, 3, std.testing.allocator); defer A.deinit(std.testing.allocator);
    var b = try real_vector.RealVector(f64).init(3, std.testing.allocator); defer b.deinit(std.testing.allocator);

    var x_expected = try real_vector.RealVector(f64).init(3, std.testing.allocator); defer x_expected.deinit(std.testing.allocator);

    A.ptr(0, 0).* = 4; A.ptr(0, 1).* = 1; A.ptr(0, 2).* = 2;
    A.ptr(1, 0).* = 1; A.ptr(1, 1).* = 3; A.ptr(1, 2).* = 0;
    A.ptr(2, 0).* = 2; A.ptr(2, 1).* = 0; A.ptr(2, 2).* = 5;

    b.ptr(0).* = 4; b.ptr(1).* = 5; b.ptr(2).* = 6;

    x_expected.ptr(0).* = -0.02325581395349; x_expected.ptr(1).* = 1.67441860465116; x_expected.ptr(2).* = 1.20930232558139;

    const x = try linearSolveSymmetricAlloc(f64, A, b, std.testing.allocator); defer x.deinit(std.testing.allocator);

    try std.testing.expect(x.eq(x_expected, TEST_TOLERANCE));
}
