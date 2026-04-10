//! File with linear solvers.

const std = @import("std");

const complex_matrix = @import("complex_matrix.zig");
const complex_vector = @import("complex_vector.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const global_variables = @import("global_variables.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const ComplexMatrix = complex_matrix.ComplexMatrix;
const ComplexVector = complex_vector.ComplexVector;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const eigensystemHermitianAlloc = eigenproblem_solver.eigensystemHermitianAlloc;
const mm = matrix_multiplication.mm;

const TEST_TOLERANCE = global_variables.TEST_TOLERANCE;

/// Solve the linear system Ax = b using the eigenvalue decomposition of A. The matrix A must be symmetric.
pub fn linearSolveHermitian(comptime T: type, x: anytype, AJ: anytype, AC: anytype, b: anytype, y: anytype) !void {
    var x_matrix = x.asMatrix();
    var y_matrix = y.asMatrix();

    if (comptime @TypeOf(AJ, AC) == RealMatrix(T) and @TypeOf(x.*, b, y.*) == RealVector(T)) {return try linearSolveHermitianSpectral(RealMatrix, T, &x_matrix, AJ, AC, b.asMatrix(), &y_matrix);}
    else if (comptime @TypeOf(AJ, AC) == ComplexMatrix(T) and @TypeOf(x.*, b, y.*) == ComplexVector(T)) {return try linearSolveHermitianSpectral(ComplexMatrix, T, &x_matrix, AJ, AC, b.asMatrix(), &y_matrix);}
    else @compileError("UNSUPPORTED MATRIX OR VECTOR TYPE IN HERMITIAN LINEAR SOLVER");
}

/// Solve the linear system Ax = b using the eigenvalue decomposition of A. the result is returned as a new vector.
pub fn linearSolveHermitianAlloc(comptime T: type, A: anytype, b: anytype, allocator: std.mem.Allocator) !@TypeOf(b) {
    const AJC = try eigensystemHermitianAlloc(T, A, allocator); defer AJC.J.deinit(allocator); defer AJC.C.deinit(allocator);

    var x = try @TypeOf(b).init(b.len, allocator); errdefer x.deinit(allocator);
    var y = try @TypeOf(b).init(b.len, allocator);    defer y.deinit(allocator);

    try linearSolveHermitian(T, &x, AJC.J, AJC.C, b, &y);

    return x;
}

/// Solve the linear system Ax = b using the eigenvalue decomposition of A. The matrix A must be symmetric.
pub fn linearSolveHermitianSpectral(comptime M: fn (comptime type) type, comptime T: type, x: *M(T), AJ: M(T), AC: M(T), b: M(T), y: *M(T)) !void {
    const tolerance = @as(T, @floatFromInt(AJ.rows)) * std.math.floatEps(T) * AJ.maxAbsDiagonal();

    for (0..AJ.rows) |i| if ((if (comptime M == ComplexMatrix) AJ.at(i, i).magnitude() else @abs(AJ.at(i, i))) <= tolerance) {
        return error.NumericalError;
    };
    
    try mm(T, y, AC, true, b, false);

    if (comptime M == RealMatrix) {
        for (0..y.rows) |i| y.ptr(i, 0).* /= AJ.at(i, i);
    }
    else {
        for (0..y.rows) |i| y.ptr(i, 0).* = y.at(i, 0).div(AJ.at(i, i));
    }

    try mm(T, x, AC, false, y.*, false);
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

    const x = try linearSolveHermitianAlloc(f64, A, b, std.testing.allocator); defer x.deinit(std.testing.allocator);

    try std.testing.expect(x.eq(x_expected, TEST_TOLERANCE));
}
