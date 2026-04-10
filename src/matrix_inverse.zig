//! File with linear solvers.

const std = @import("std");

const complex_matrix = @import("complex_matrix.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const global_variables = @import("global_variables.zig");
const real_matrix = @import("real_matrix.zig");

const Complex = std.math.complex.Complex;
const ComplexMatrix = complex_matrix.ComplexMatrix;
const RealMatrix = real_matrix.RealMatrix;

const eigensystemHermitianAlloc = eigenproblem_solver.eigensystemHermitianAlloc;

const TEST_TOLERANCE = global_variables.TEST_TOLERANCE;

/// Form the inverse of a symmetric or hermitian matrix A using its eigenvalue decomposition. Ainv = C * J_inv * C^T
pub fn inverseHermitian(comptime T: type, Ainv: anytype, AJ: anytype, AC: anytype) !void {
    const tolerance = @as(T, @floatFromInt(AJ.rows)) * std.math.floatEps(T) * AJ.maxAbsDiagonal();

    for (0..AJ.rows) |i| if (@abs(AJ.at(i, i)) < tolerance) {

        std.log.err("MATRIX PASSED TO INVERSE FUNCTION IS SINGULAR OR NEARLY SINGULAR", .{});

        return error.NumericalError;
    };

    if (comptime @TypeOf(AJ, AC) == RealMatrix(T)) {return try pseudoInverseHermitianSpectral(RealMatrix, T, Ainv, AJ, AC, tolerance);} 
    else if (comptime @TypeOf(AJ, AC) == ComplexMatrix(T)) {return try pseudoInverseHermitianSpectral(ComplexMatrix, T, Ainv, AJ, AC, tolerance);} 
    else @compileError("UNSUPPORTED MATRIX TYPE IN PSEUDO INVERSE FUNCTION");
}

/// Form the inverse of a symmetric or hermitian matrix A using its eigenvalue decomposition. Ainv = C * J_inv * C^T
pub fn inverseHermitianAlloc(comptime T: type, A: anytype, allocator: std.mem.Allocator) !@TypeOf(A) {
    var Ainv = try @TypeOf(A).init(A.rows, A.cols, allocator);

    const AJC = try eigensystemHermitianAlloc(T, A, allocator); defer AJC.J.deinit(allocator); defer AJC.C.deinit(allocator);

    try inverseHermitian(T, &Ainv, AJC.J, AJC.C);

    return Ainv;
}

/// Form the pseudo inverse of a hermitian matrix A using its eigenvalue decomposition. Ainv = C * J_inv * C^T
pub fn pseudoInverseHermitian(comptime T: type, Ainv: anytype, AJ: anytype, AC: anytype, thresh: T) !void {
    if (comptime @TypeOf(AJ, AC) == RealMatrix(T)) {return try pseudoInverseHermitianSpectral(RealMatrix, T, Ainv, AJ, AC, thresh);} 
    else if (comptime @TypeOf(AJ, AC) == ComplexMatrix(T)) {return try pseudoInverseHermitianSpectral(ComplexMatrix, T, Ainv, AJ, AC, thresh);} 
    else @compileError("UNSUPPORTED MATRIX TYPE IN PSEUDO INVERSE FUNCTION");
}

/// Form the pseudo inverse of a symmetric or hermitian matrix A using its eigenvalue decomposition. Ainv = C * J_inv * C^T
pub fn pseudoInverseHermitianSpectral(comptime M: fn (comptime type) type, comptime T: type, Ainv: *M(T), AJ: M(T), AC: M(T), thresh: T) !void {
    for (0..Ainv.rows) |i| {

        if ((if (comptime M == RealMatrix) @abs(AJ.at(i, i)) else std.math.complex.abs(AJ.at(i, i))) < thresh) {
            Ainv.ptr(Ainv.rows - 1, i).* = if (comptime M == RealMatrix) 0 else Complex(T).init(0, 0); continue;
        }

        Ainv.ptr(Ainv.rows - 1, i).* = if (comptime M == RealMatrix) 1 / AJ.at(i, i) else Complex(T).init(1, 0).div(AJ.at(i, i));
    }

    for(0..Ainv.rows) |i|  for(i..Ainv.cols) |j|  {

        var sum = if (comptime M == RealMatrix) @as(T, @floatCast(0)) else Complex(T).init(0, 0);

        for(0..Ainv.rows) |k|  {

            if (comptime M == RealMatrix) {
                sum += AC.at(i, k) * Ainv.at(Ainv.rows - 1, k) * AC.at(j, k); 
            } else {
                sum = sum.add(AC.at(i, k).mul(Ainv.at(Ainv.rows - 1, k)).mul(AC.at(j, k).conjugate()));
            }
        }

        Ainv.ptr(i, j).* = sum;
    };

    for (0..Ainv.rows) |i| for (0..i) |j| {
        Ainv.ptr(i, j).* = if (comptime M == RealMatrix) Ainv.at(j, i) else Ainv.at(j, i).conjugate();
    };
}

test "Symmetric 3x3 Inverse" {
    var A = try RealMatrix(f64).init(3, 3, std.testing.allocator); defer A.deinit(std.testing.allocator);

    var Ainv_expected = try RealMatrix(f64).init(3, 3, std.testing.allocator); defer Ainv_expected.deinit(std.testing.allocator);

    A.ptr(0, 0).* = 1.0; A.ptr(0, 1).* = 2.0; A.ptr(0, 2).* = 3.0;
    A.ptr(1, 0).* = 2.0; A.ptr(1, 1).* = 4.0; A.ptr(1, 2).* = 5.0;
    A.ptr(2, 0).* = 3.0; A.ptr(2, 1).* = 5.0; A.ptr(2, 2).* = 6.0;

    Ainv_expected.ptr(0, 0).* =  1.0; Ainv_expected.ptr(0, 1).* = -3.0; Ainv_expected.ptr(0, 2).* =  2.0;
    Ainv_expected.ptr(1, 0).* = -3.0; Ainv_expected.ptr(1, 1).* =  3.0; Ainv_expected.ptr(1, 2).* = -1.0;
    Ainv_expected.ptr(2, 0).* =  2.0; Ainv_expected.ptr(2, 1).* = -1.0; Ainv_expected.ptr(2, 2).* =  0.0;

    const Ainv = try inverseHermitianAlloc(f64, A, std.testing.allocator); defer Ainv.deinit(std.testing.allocator);

    try std.testing.expect(Ainv.eq(Ainv_expected, TEST_TOLERANCE));
}
