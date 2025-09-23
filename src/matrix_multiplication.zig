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
    if (comptime (@TypeOf(A) != RealMatrix(T) and @TypeOf(A) != ComplexMatrix(T)) or (@TypeOf(B) != RealMatrix(T) and @TypeOf(B) != ComplexMatrix(T))) {
        return throw(void, "UNSUPPORTED INPUT MATRIX TYPE A IN MATRIX MULTIPLICATION", .{});
    }

    if (comptime !at and !bt) if (A.cols != B.rows) {
        return throw(void, "CAN'T MULTIPLY A MATRIX WITH {d} COLUMNS BY A MATRIX WITH {d} ROWS", .{A.cols, B.rows});
    } else if (comptime at and !bt) if (A.rows != B.rows) {
        return throw(void, "CAN'T MULTIPLY A TRANSPOSED MATRIX WITH {d} COLUMNS BY A MATRIX WITH {d} ROWS", .{A.rows, B.rows});
    } else if (comptime !at and bt) if (A.cols != B.cols) {
        return throw(void, "CAN'T MULTIPLY A MATRIX WITH {d} COLUMNS BY A MATRIX WITH {d} ROWS", .{A.cols, B.cols});
    } else if (comptime at and bt) if (A.rows != B.cols) {
        return throw(void, "CAN'T MULTIPLY A TRANSPOSED MATRIX WITH {d} COLUMNS BY A MATRIX WITH {d} ROWS", .{A.rows, B.cols});
    };

    if (comptime @TypeOf(C) == *RealMatrix(T)) {
        return mmReal(T, C, A, at, B, bt);
    } else if (comptime @TypeOf(C) == *ComplexMatrix(T)) {
        return mmComplex(T, C, A, at, B, bt);
    }

    else comptime return throw(void, "UNSUPPORTED OUTPUT MATRIX TYPE IN MATRIX MULTIPLICATION", .{});
}

/// Matrix multiplication function between two complex matrices A and B, with options to transpose A and/or B. The result is stored in C.
pub fn mmComplex(comptime T: type, C: *ComplexMatrix(T), A: anytype, comptime at: bool, B: anytype, comptime bt: bool) !void {
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
pub fn mmReal(comptime T: type, C: *RealMatrix(T), A: RealMatrix(T), comptime at: bool, B: RealMatrix(T), comptime bt: bool) !void {
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
