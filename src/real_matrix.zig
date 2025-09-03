//! Matrix struct that uses real numbers (f64, i32, usize, ...) for its elements.

const std = @import("std");

const real_vector = @import("real_vector.zig");

const RealVector = real_vector.RealVector;

/// Real matrix class. The matrix is stored in a flat array in row-major order.
pub fn RealMatrix(comptime T: type) type {
    return struct {
        data: []T,
        rows: usize,
        cols: usize,

        allocator: std.mem.Allocator,

        /// Initialize a matrix with a given number of rows and columns and specify an allocator. The function returns an error if the allocation fails.
        pub fn init(rows: usize, cols: usize, allocator: std.mem.Allocator) !@This() {
            return @This(){
                .data = try allocator.alloc(T, rows * cols),
                .rows = rows,
                .cols = cols,
                .allocator = allocator
            };
        }

        /// Initialize a matrix and fills it with zeros.
        pub fn initZero(rows: usize, cols: usize, allocator: std.mem.Allocator) !@This() {
            var A = try @This().init(rows, cols, allocator); A.zero();

            return A;
        }

        /// Free the memory allocated for the matrix.
        pub fn deinit(self: @This()) void {
            self.allocator.free(self.data);
        }

        /// Add another matrix to this matrix.
        pub fn add(self: *@This(), other: @This()) void {
            std.debug.assert(self.rows == other.rows and self.cols == other.cols);

            for (self.data, 0..) |*element, index| {
                element.* += other.data[index];
            }
        }

        /// Get the element at (i, j).
        pub fn at(self: @This(), i: usize, j: usize) T {
            return self.data[i * self.cols + j];
        }

        /// Divide the matrix by a scalar value.
        pub fn divs(self: *@This(), value: T) void {
            for (self.data) |*element| element.* /= value;
        }

        /// Equality operator for matrices, given a tolerance.
        pub fn eq(self: @This(), other: @This(), tol: T) bool {
            if (self.rows != other.rows or self.cols != other.cols) return false;

            for (self.data, 0..) |element, index| {
                if (@abs(element - other.data[index]) > tol) return false;
            }

            return true;
        }

        /// Fill the matrix with a given value.
        pub fn fill(self: *@This(), value: T) void {
            @memset(self.data, value);
        }

        /// Calculate the Frobenius norm of the matrix.
        pub fn frobeniusNorm(self: @This()) T {
            return std.math.sqrt(self.frobeniusNormSquared());
        }

        /// Calculate square of the Frobenius norm of the matrix.
        pub fn frobeniusNormSquared(self: @This()) T {
            var sum: T = 0;

            for (self.data) |element| {
                sum += element * element;
            }

            return sum;
        }

        /// Set the matrix to the identity matrix. The matrix must be square.
        pub fn identity(self: *@This()) void {
            std.debug.assert(self.isSquare());

            self.zero();

            for (0..self.rows) |i| {
                self.ptr(i, i).* = 1;
            }
        }

        /// Square checker.
        pub fn isSquare(self: @This()) bool {
            return self.rows == self.cols;
        }

        /// Calculate the Frobenius norm of the off-diagonal elements of the matrix.
        pub fn offDiagonalFrobeniusNorm(self: @This()) T {
            return std.math.sqrt(self.offDiagonalFrobeniusNormSquared());
        }

        /// Calculate square of the Frobenius norm of the off-diagonal elements of the matrix.
        pub fn offDiagonalFrobeniusNormSquared(self: @This()) T {
            var sum: T = 0;

            for (0..self.rows) |i| for (0..self.cols) |j| if (i != j) {
                sum += self.at(i, j) * self.at(i, j);
            };

            return sum;
        }

        /// Get the pointer to the element at (i, j).
        pub fn ptr(self: *@This(), i: usize, j: usize) *T {
            return &self.data[i * self.cols + j];
        }

        /// Get a row as view to a a real vector.
        pub fn row(self: @This(), i: usize) RealVector(T) {
            return RealVector(T){
                .data = self.data[i * self.cols..(i + 1) * self.cols],
                .len = self.cols,
                .allocator = self.allocator
            };
        }

        /// Fills the matrix with zeros.
        pub fn zero(self: *@This()) void {
            self.fill(0);
        }
    };
}

test "init, deinit" {
    var A = try RealMatrix(f64).init(345, 753, std.testing.allocator); defer A.deinit();

    try std.testing.expect(A.rows == 345);
    try std.testing.expect(A.cols == 753);
    try std.testing.expect(A.data.len == 259785);
}
