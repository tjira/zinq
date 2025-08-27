//! Matrix struct that uses real numbers (f64, i32, usize, ...) for its elements.

const std = @import("std");

/// Real matrix class. The matrix is stored in a flat array in row-major order.
pub fn RealMatrix(comptime T: type) type {
    return struct {
        data: []T,
        rows: usize,
        cols: usize,
        allocator: std.mem.Allocator,

        /// Initialize a matrix with a given number of rows and columns and specify an allocator. The function returns an error if the allocation fails.
        pub fn init(rows: usize, cols: usize, allocator: std.mem.Allocator) !RealMatrix(T) {
            return RealMatrix(T){
                .data = try allocator.alloc(T, rows * cols),
                .rows = rows,
                .cols = cols,
                .allocator = allocator
            };
        }

        /// Initialize a matrix and fills it with zeros.
        pub fn initZero(rows: usize, cols: usize, allocator: std.mem.Allocator) !RealMatrix(T) {
            var A = try RealMatrix(T).init(rows, cols, allocator); A.zero();

            return A;
        }

        /// Free the memory allocated for the matrix.
        pub fn deinit(self: RealMatrix(T)) void {
            self.allocator.free(self.data);
        }

        /// Get the element at (i, j).
        pub fn at(self: RealMatrix(T), i: usize, j: usize) T {
            return self.data[i * self.cols + j];
        }

        /// Equality operator for matrices.
        pub fn eq(self: RealMatrix(T), other: RealMatrix(T)) bool {
            if (self.rows != other.rows or self.cols != other.cols) return false;

            for (self.data, 0..) |element, index| {
                if (element != other.data[index]) return false;
            }

            return true;
        }

        /// Fill the matrix with a given value.
        pub fn fill(self: *RealMatrix(T), value: T) void {
            @memset(self.data, value);
        }

        /// Square checker.
        pub fn isSquare(self: RealMatrix(T)) bool {
            return self.rows == self.cols;
        }

        /// Get the pointer to the element at (i, j).
        pub fn ptr(self: *RealMatrix(T), i: usize, j: usize) *T {
            return &self.data[i * self.cols + j];
        }

        /// Fills the matrix with zeros.
        pub fn zero(self: *RealMatrix(T)) void {
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
