//! Algorithms for matrix multiplication.

const std = @import("std");

const error_handling = @import("error_handling.zig");
const real_matrix = @import("real_matrix.zig");
const complex_matrix = @import("complex_matrix.zig");
const global_variables = @import("global_variables.zig");

const RealMatrix = real_matrix.RealMatrix;
const ComplexMatrix = complex_matrix.ComplexMatrix;
const Complex = std.math.complex.Complex;

const throw = error_handling.throw;

const TEST_TOLERANCE = global_variables.TEST_TOLERANCE;

/// General matrix multiplication function.
pub fn mm(comptime T: type, C: anytype, A: anytype, comptime at: bool, B: anytype, comptime bt: bool) !void {
    if (comptime @TypeOf(C.*) != mmResultType(T, @TypeOf(A), @TypeOf(B))) {
        @compileError("OUTPUT MATRIX TYPE DOESN'T MATCH THE EXPECTED TYPE FROM MULTIPLYING THE INPUT MATRICES");
    }

    if (comptime !at and !bt) if (A.cols != B.rows) {
        return throw(void, "CAN'T MULTIPLY A MATRIX WITH {d} COLUMNS BY A MATRIX WITH {d} ROWS", .{A.cols, B.rows});
    };
    if (comptime at and !bt) if (A.rows != B.rows) {
        return throw(void, "CAN'T MULTIPLY A TRANSPOSED MATRIX WITH {d} COLUMNS BY A MATRIX WITH {d} ROWS", .{A.rows, B.rows});
    };
    if (comptime !at and bt) if (A.cols != B.cols) {
        return throw(void, "CAN'T MULTIPLY A MATRIX WITH {d} COLUMNS BY A MATRIX WITH {d} ROWS", .{A.cols, B.cols});
    };
    if (comptime at and bt) if (A.rows != B.cols) {
        return throw(void, "CAN'T MULTIPLY A TRANSPOSED MATRIX WITH {d} COLUMNS BY A MATRIX WITH {d} ROWS", .{A.rows, B.cols});
    };

    if (comptime @TypeOf(C) == *RealMatrix(T)) return mmReal(T, C, A, at, B, bt);
    if (comptime @TypeOf(C) == *ComplexMatrix(T)) return mmComplex(T, C, A, at, B, bt);

    return throw(void, "UNSUPPORTED OUTPUT MATRIX TYPE IN MATRIX MULTIPLICATION", .{});
}

/// General matrix multiplication function that returns a new matrix.
pub fn mmAlloc(comptime T: type, A: anytype, comptime at: bool, B: anytype, comptime bt: bool, allocator: std.mem.Allocator) !mmResultType(T, @TypeOf(A), @TypeOf(B)) {
    var C = try mmResultType(T, @TypeOf(A), @TypeOf(B)).init(A.rows, B.cols, allocator);

    try mm(T, &C, A, at, B, bt);

    return C;
}

/// Matrix multiplication function between two complex matrices A and B, with options to transpose A and/or B. The result is stored in C.
pub fn mmComplex(comptime T: type, C: *ComplexMatrix(T), A: anytype, comptime at: bool, B: anytype, comptime bt: bool) void {
    for (0..A.rows) |i| for (0..B.cols) |j| {

        var sum = Complex(T).init(0, 0);

        for (0..A.cols) |k| {

            const a = if (comptime at) A.at(k, i) else A.at(i, k);
            const b = if (comptime bt) B.at(j, k) else B.at(k, j);

            const ac = if (comptime @TypeOf(a) == Complex(T)) a else Complex(T).init(a, 0);
            const bc = if (comptime @TypeOf(b) == Complex(T)) b else Complex(T).init(b, 0);

            sum = sum.add(ac.mul(bc));
        }

        C.ptr(i, j).* = sum;
    };
}

/// Matrix multiplication function between two real matrices A and B, with options to transpose A and/or B. The result is stored in C.
pub fn mmReal(comptime T: type, C: *RealMatrix(T), A: RealMatrix(T), comptime at: bool, B: RealMatrix(T), comptime bt: bool) void {
    for (0..A.rows) |i| for (0..B.cols) |j| {

        var sum: T = 0;

        for (0..A.cols) |k| {

            const a = if (comptime at) A.at(k, i) else A.at(i, k);
            const b = if (comptime bt) B.at(j, k) else B.at(k, j);

            sum += a * b;
        }

        C.ptr(i, j).* = sum;
    };
}

/// Returns the type of the result matrix from multiplying two matrices A and B.
pub fn mmResultType(comptime T: type, comptime A_type: type, comptime B_type: type) type {
    if (comptime (A_type != RealMatrix(T) and A_type != ComplexMatrix(T)) or (B_type != RealMatrix(T) and B_type != ComplexMatrix(T))) {
        @compileError("UNSUPPORTED INPUT MATRIX TYPE IN MATRIX MULTIPLICATION");
    }

    return if (comptime A_type == RealMatrix(T) and B_type == RealMatrix(T)) RealMatrix(T) else ComplexMatrix(T);
}

test "3x3 Real Matrix Multiplication" {
    var A = try RealMatrix(f64).init(3, 3, std.testing.allocator); defer A.deinit();
    var B = try RealMatrix(f64).init(3, 3, std.testing.allocator); defer B.deinit();

    var C_expected = try RealMatrix(f64).init(3, 3, std.testing.allocator); defer C_expected.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2; A.ptr(0, 2).* = 3;
    A.ptr(1, 0).* = 4; A.ptr(1, 1).* = 5; A.ptr(1, 2).* = 6;
    A.ptr(2, 0).* = 7; A.ptr(2, 1).* = 8; A.ptr(2, 2).* = 9;

    B.ptr(0, 0).* = 9; B.ptr(0, 1).* = 8; B.ptr(0, 2).* = 7;
    B.ptr(1, 0).* = 6; B.ptr(1, 1).* = 5; B.ptr(1, 2).* = 4;
    B.ptr(2, 0).* = 3; B.ptr(2, 1).* = 2; B.ptr(2, 2).* = 1;

    C_expected.ptr(0, 0).* =  30; C_expected.ptr(0, 1).* =  24; C_expected.ptr(0, 2).* = 18;
    C_expected.ptr(1, 0).* =  84; C_expected.ptr(1, 1).* =  69; C_expected.ptr(1, 2).* = 54;
    C_expected.ptr(2, 0).* = 138; C_expected.ptr(2, 1).* = 114; C_expected.ptr(2, 2).* = 90;

    const C = try mmAlloc(f64, A, false, B, false, std.testing.allocator); defer C.deinit();

    try std.testing.expect(C.eq(C_expected, TEST_TOLERANCE));
}
