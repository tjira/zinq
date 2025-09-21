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

/// Multiply two complex matrices A and B and store the result in C.
pub fn mmComplexComplex(comptime T: type, C: *ComplexMatrix(T), A: ComplexMatrix(T), B: ComplexMatrix(T)) !void {
    if (A.cols != B.rows) {
        return throw(void, "CAN'T MULTIPLY A MATRIX WITH {d} COLUMNS BY A MATRIX WITH {d} ROWS", .{A.cols, B.rows});
    }

    for (0..A.rows) |i| for (0..B.cols) |j| {

        var sum = Complex(T).init(0, 0);

        for (0..A.cols) |k| {
            sum = sum.add(A.at(i, k).mul(B.at(k, j)));
        }

        C.ptr(i, j).* = sum;
    };
}

/// Multiply two complex matrices A and B and return the result. The function returns an error if the allocation fails.
pub fn mmComplexComplexAlloc(comptime T: type, A: ComplexMatrix(T), B: ComplexMatrix(T)) !ComplexMatrix(T) {
    var C = try ComplexMatrix(T).init(A.rows, B.cols, A.allocator);

    try mmComplexComplex(T, &C, A, B);

    return C;
}

/// Multiply one complex matrix A by a real matrix B and store the result in C.
pub fn mmComplexReal(comptime T: type, C: *ComplexMatrix(T), A: ComplexMatrix(T), B: RealMatrix(T)) !void {
    if (A.cols != B.rows) {
        return throw(void, "CAN'T MULTIPLY A MATRIX WITH {d} COLUMNS BY A MATRIX WITH {d} ROWS", .{A.cols, B.rows});
    }

    for (0..A.rows) |i| for (0..B.cols) |j| {

        var sum = Complex(T).init(0, 0);

        for (0..A.cols) |k| {
            sum = sum.add(A.at(i, k).mul(Complex(T).init(B.at(k, j), 0)));
        }

        C.ptr(i, j).* = sum;
    };
}

/// Multiply one complex matrix A by a transposed real matrix BT and store the result in C.
pub fn mmComplexRealTrans(comptime T: type, C: *ComplexMatrix(T), A: ComplexMatrix(T), BT: RealMatrix(T)) !void {
    if (A.cols != BT.cols) {
        return throw(void, "CAN'T MULTIPLY A MATRIX WITH {d} COLUMNS BY A MATRIX WITH {d} ROWS", .{A.cols, BT.cols});
    }

    for (0..A.rows) |i| for (0..BT.rows) |j| {

        var sum = Complex(T).init(0, 0);

        for (0..A.cols) |k| {
            sum = sum.add(A.at(i, k).mul(Complex(T).init(BT.at(j, k), 0)));
        }

        C.ptr(i, j).* = sum;
    };
}

/// Multiply one real matrix A by a complex matrix B and store the result in C.
pub fn mmRealComplex(comptime T: type, C: *ComplexMatrix(T), A: RealMatrix(T), B: ComplexMatrix(T)) !void {
    if (A.cols != B.rows) {
        return throw(void, "CAN'T MULTIPLY A MATRIX WITH {d} COLUMNS BY A MATRIX WITH {d} ROWS", .{A.cols, B.rows});
    }

    for (0..A.rows) |i| for (0..B.cols) |j| {

        var sum = Complex(T).init(0, 0);

        for (0..A.cols) |k| {
            sum = sum.add(Complex(T).init(A.at(i, k), 0).mul(B.at(k, j)));
        }

        C.ptr(i, j).* = sum;
    };
}

/// Multiply two matrices A and B and store the result in C.
pub fn mmRealReal(comptime T: type, C: *RealMatrix(T), A: RealMatrix(T), B: RealMatrix(T)) !void {
    if (A.cols != B.rows) {
        return throw(void, "CAN'T MULTIPLY A MATRIX WITH {d} COLUMNS BY A MATRIX WITH {d} ROWS", .{A.cols, B.rows});
    }

    for (0..A.rows) |i| for (0..B.cols) |j| {

        var sum: T = 0;

        for (0..A.cols) |k| {
            sum += A.at(i, k) * B.at(k, j);
        }

        C.ptr(i, j).* = sum;
    };
}

/// Multiply two matrices A and B and return the result. The function returns an error if the allocation fails.
pub fn mmRealRealAlloc(comptime T: type, A: RealMatrix(T), B: RealMatrix(T)) !RealMatrix(T) {
    var C = try RealMatrix(T).init(A.rows, B.cols, A.allocator);

    try mmRealReal(T, &C, A, B);

    return C;
}

/// Multiply one real matrix A by a transposed real matrix BT and store the result in C.
pub fn mmRealRealTrans(comptime T: type, C: *RealMatrix(T), A: RealMatrix(T), BT: RealMatrix(T)) !void {
    if (A.cols != BT.cols) {
        return throw(void, "CAN'T MULTIPLY A MATRIX WITH {d} COLUMNS BY A MATRIX WITH {d} ROWS", .{A.cols, BT.cols});
    }

    for (0..A.rows) |i| for (0..BT.rows) |j| {

        var sum: T = 0;

        for (0..A.cols) |k| {
            sum += A.at(i, k) * BT.at(j, k);
        }

        C.ptr(i, j).* = sum;
    };
}

/// Multiply two matrices A and B^T and return the result. The function returns an error if the allocation fails.
pub fn mmRealRealTransAlloc(comptime T: type, A: RealMatrix(T), B: RealMatrix(T)) !RealMatrix(T) {
    var C = try RealMatrix(T).init(A.rows, B.rows, A.allocator);

    try mmRealRealTrans(T, &C, A, B);

    return C;
}

/// Multiply one transposed real matrix A by a complex matrix B and store the result in C.
pub fn mmRealTransComplex(comptime T: type, C: *ComplexMatrix(T), AT: RealMatrix(T), B: ComplexMatrix(T)) !void {
    if (AT.rows != B.rows) {
        return throw(void, "CAN'T MULTIPLY A MATRIX WITH {d} COLUMNS BY A MATRIX WITH {d} ROWS", .{AT.rows, B.rows});
    }

    for (0..AT.cols) |i| for (0..B.cols) |j| {

        var sum = Complex(T).init(0, 0);

        for (0..AT.rows) |k| {
            sum = sum.add(Complex(T).init(AT.at(k, i), 0).mul(B.at(k, j)));
        }

        C.ptr(i, j).* = sum;
    };
}

/// Multiply one transposed real matrix A by a real matrix B and store the result in C.
pub fn mmRealTransReal(comptime T: type, C: *RealMatrix(T), AT: RealMatrix(T), B: RealMatrix(T)) !void {
    if (AT.rows != B.rows) {
        return throw(void, "CAN'T MULTIPLY A MATRIX WITH {d} COLUMNS BY A MATRIX WITH {d} ROWS", .{AT.rows, B.rows});
    }

    for (0..AT.cols) |i| for (0..B.cols) |j| {

        var sum: T = 0;

        for (0..AT.rows) |k| {
            sum += AT.at(k, i) * B.at(k, j);
        }

        C.ptr(i, j).* = sum;
    };
}

/// Multiply two matrices A^T and B and return the result. The function returns an error if the allocation fails.
pub fn mmRealTransRealAlloc(comptime T: type, A: RealMatrix(T), B: RealMatrix(T)) !RealMatrix(T) {
    var C = try RealMatrix(T).init(A.cols, B.cols, A.allocator);

    try mmRealTransReal(T, &C, A, B);

    return C;
}

test "mmRealRealAlloc" {
    var A = try RealMatrix(f64).init(2, 2, std.testing.allocator); defer A.deinit();
    var B = try RealMatrix(f64).init(2, 2, std.testing.allocator); defer B.deinit();
    var C = try RealMatrix(f64).init(2, 2, std.testing.allocator); defer C.deinit();

    A.ptr(0, 0).* = 1; A.ptr(0, 1).* = 2; A.ptr(1, 0).* = 3; A.ptr(1, 1).* = 4;
    B.ptr(0, 0).* = 5; B.ptr(0, 1).* = 6; B.ptr(1, 0).* = 7; B.ptr(1, 1).* = 8;

    C.ptr(0, 0).* = 19; C.ptr(0, 1).* = 22; C.ptr(1, 0).* = 43; C.ptr(1, 1).* = 50;

    var D = try mmRealRealAlloc(f64, A, B); defer D.deinit();

    try std.testing.expect(C.eq(D, TEST_TOLERANCE));
}
