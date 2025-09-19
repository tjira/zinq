//! Matrix struct that uses complex numbers for its elements.

const std = @import("std");

const complex_vector = @import("complex_vector.zig");

const Complex = std.math.complex.Complex;
const ComplexVector = complex_vector.ComplexVector;

/// Complex matrix class. The matrix is stored in a flat array in row-major order.
pub fn ComplexMatrix(comptime T: type) type {
    return struct {
        data: []Complex(T),
        rows: usize,
        cols: usize,

        allocator: ?std.mem.Allocator,

        /// Initialize a matrix with a given number of rows and columns and specify an allocator. The function returns an error if the allocation fails.
        pub fn init(rows: usize, cols: usize, allocator: ?std.mem.Allocator) !@This() {
            return @This(){
                .data = try allocator.?.alloc(Complex(T), rows * cols),
                .rows = rows,
                .cols = cols,
                .allocator = allocator
            };
        }

        /// Initialize a matrix and fills it with zeros.
        pub fn initZero(rows: usize, cols: usize, allocator: ?std.mem.Allocator) !@This() {
            var A = try @This().init(rows, cols, allocator); A.zero();

            return A;
        }

        /// Free the memory allocated for the matrix.
        pub fn deinit(self: @This()) void {
            if (self.allocator) |allocator| allocator.free(self.data);
        }

        /// Divide the matrix by a scalar value.
        pub fn divs(self: *@This(), value: Complex(T)) void {
            for (self.data) |*element| {
                element.* = element.div(value);
            }
        }

        /// Get the element at (i, j).
        pub fn at(self: @This(), i: usize, j: usize) Complex(T) {
            return self.data[i * self.cols + j];
        }

        /// Equality operator for matrices, given a tolerance.
        pub fn eq(self: @This(), other: @This(), tol: T) bool {
            if (self.rows != other.rows or self.cols != other.cols) return false;

            for (self.data, 0..) |element, index| {
                if (@abs(element.re - other.data[index].re) > tol) return false;
                if (@abs(element.im - other.data[index].im) > tol) return false;
            }

            return true;
        }

        /// Fill the matrix with a given value.
        pub fn fill(self: *@This(), value: Complex(T)) void {
            for (self.data) |*element| element.* = value;
        }

        /// Multiply the matrix by a scalar value.
        pub fn muls(self: *@This(), value: Complex(T)) void {
            for (self.data) |*element| {
                element.* = element.mul(value);
            }
        }

        /// Get the pointer to the element at (i, j).
        pub fn ptr(self: *@This(), i: usize, j: usize) *Complex(T) {
            return &self.data[i * self.cols + j];
        }

        /// Returns the row of the matrix as a complex vector.
        pub fn row(self: @This(), i: usize) ComplexVector(T) {
            return ComplexVector(T){
                .data = self.data[i * self.cols .. (i + 1) * self.cols],
                .len = self.cols,
                .allocator = null
            };
        }

        /// Fills the matrix with zeros.
        pub fn zero(self: *@This()) void {
            self.fill(Complex(T).init(0, 0));
        }
    };
}

test "init, deinit" {
    var A = try ComplexMatrix(f64).init(345, 753, std.testing.allocator); defer A.deinit();

    try std.testing.expect(A.rows == 345);
    try std.testing.expect(A.cols == 753);
    try std.testing.expect(A.data.len == 259785);
}
