//! Matrix struct that uses real numbers (f64, i32, usize, ...) for its elements.

const std = @import("std");

const error_handling = @import("error_handling.zig");
const real_vector = @import("real_vector.zig");
const strided_real_vector = @import("strided_real_vector.zig");

const RealVector = real_vector.RealVector;
const StridedRealVector = strided_real_vector.StridedRealVector;

const throw = error_handling.throw;

/// Real matrix class. The matrix is stored in a flat array in row-major order.
pub fn RealMatrix(comptime T: type) type {
    return struct {
        data: []T,
        rows: usize,
        cols: usize,

        allocator: ?std.mem.Allocator,

        /// Initialize a matrix with a given number of rows and columns and specify an allocator. The function returns an error if the allocation fails.
        pub fn init(rows: usize, cols: usize, allocator: ?std.mem.Allocator) !@This() {
            return @This(){
                .data = try allocator.?.alloc(T, rows * cols),
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

        /// Add another matrix to this matrix.
        pub fn add(self: *@This(), other: @This()) !void {
            if (self.rows != other.rows or self.cols != other.cols) {
                return throw(void, "CAN'T ADD A {d}x{d} MATRIX TO A {d}x{d} MATRIX", .{other.rows, other.cols, self.rows, self.cols});
            }

            for (self.data, 0..) |*element, index| {
                element.* += other.data[index];
            }
        }

        /// Get the element at (i, j).
        pub fn at(self: @This(), i: usize, j: usize) T {
            return self.data[i * self.cols + j];
        }

        /// Returns the matrix as a view to a real vector.
        pub fn asVector(self: @This()) RealVector(T) {
            return RealVector(T){
                .data = self.data,
                .len = self.rows * self.cols,
                .allocator = null
            };
        }

        /// Clone the matrix.
        pub fn clone(self: @This()) !@This() {
            var B = try @This().init(self.rows, self.cols, self.allocator);

            try self.copyTo(&B);

            return B;
        }

        /// Returns the column as a view to a strided real vector.
        pub fn column(self: @This(), j: usize) StridedRealVector(T) {
            return StridedRealVector(T){
                .data = self.data,
                .len = self.rows,
                .stride = self.cols,
                .zero = j
            };
        }

        /// Copy the contents of this matrix to another matrix.
        pub fn copyTo(self: @This(), other: *@This()) !void {
            if (self.rows != other.rows or self.cols != other.cols) {
                return throw(void, "CAN'T COPY A {d}x{d} MATRIX TO A {d}x{d} MATRIX", .{self.rows, self.cols, other.rows, other.cols});
            }

            for (self.data, 0..) |element, index| {
                other.data[index] = element;
            }
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
            for (self.data) |*element| element.* = value;
        }

        /// Calculate the Frobenius norm of the matrix.
        pub fn frobeniusNorm(self: @This()) T {
            var sum: T = 0;

            for (self.data) |element| {
                sum += element * element;
            }

            return std.math.sqrt(sum);
        }

        /// Set the matrix to the identity matrix. The matrix must be square.
        pub fn identity(self: *@This()) void {
            self.fill(0);

            for (0..@min(self.rows, self.cols)) |i| {
                self.ptr(i, i).* = 1;
            }
        }

        /// Square checker.
        pub fn isSquare(self: @This()) bool {
            return self.rows == self.cols;
        }

        /// Symmetric checker.
        pub fn isSymmetric(self: @This(), tol: T) bool {
            if (!self.isSquare()) return false;

            for (0..self.rows) |i| for (i + 1..self.cols) |j| {
                if (@abs(self.at(i, j) - self.at(j, i)) > tol) return false;
            };

            return true;
        }

        /// Multiply the matrix by a scalar value.
        pub fn muls(self: *@This(), value: T) void {
            for (self.data) |*element| element.* *= value;
        }

        /// Calculate the Frobenius norm of the off-diagonal elements of the matrix.
        pub fn offDiagonalFrobeniusNorm(self: @This()) T {
            var sum: T = 0;

            for (0..self.rows) |i| for (0..self.cols) |j| if (i != j) {
                sum += self.at(i, j) * self.at(i, j);
            };

            return std.math.sqrt(sum);
        }

        /// Get the pointer to the element at (i, j).
        pub fn ptr(self: *@This(), i: usize, j: usize) *T {
            return &self.data[i * self.cols + j];
        }

        /// Reshape the matrix.
        pub fn reshape(self: *@This(), m: usize, n: usize) !void {
            if (m * n != self.rows * self.cols) return throw(void, "CAN'T RESHAPE A {d}x{d} MATRIX INTO A {d}x{d} MATRIX", .{self.rows, self.cols, m, n});

            self.rows = m;
            self.cols = n;
        }

        /// Get a row as view to a a real vector.
        pub fn row(self: @This(), i: usize) RealVector(T) {
            return RealVector(T){
                .data = self.data[i * self.cols..(i + 1) * self.cols],
                .len = self.cols,
                .allocator = null
            };
        }

        /// Symmetrize the matrix: A = 0.5 * (A + A^T)
        pub fn symmetrize(self: *@This()) !void {
            if (!self.isSquare()) {
                return throw(void, "CAN'T SYMMETRIZE A NON-SQUARE {d}x{d} MATRIX", .{self.rows, self.cols});
            }

            for (0..self.rows) |i| for (i + 1..self.cols) |j| {

                const avg = 0.5 * (self.at(i, j) + self.at(j, i));

                self.ptr(i, j).* = avg; self.ptr(j, i).* = avg;
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
