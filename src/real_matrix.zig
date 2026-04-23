//! Matrix struct that uses real numbers (f64, i32, usize, ...) for its elements.

const std = @import("std");

const real_vector = @import("real_vector.zig");
const strided_real_vector = @import("strided_real_vector.zig");

const RealVector = real_vector.RealVector;
const StridedRealVector = strided_real_vector.StridedRealVector;

/// Real matrix class. The matrix is stored in a flat array in row-major order.
pub fn RealMatrix(comptime T: type) type {
    return struct {
        data: []T,
        rows: usize,
        cols: usize,

        /// Initialize a matrix with a given number of rows and columns and specify an allocator. The function returns an error if the allocation fails.
        pub fn init(rows: usize, cols: usize, allocator: std.mem.Allocator) !@This() {
            return @This(){ .data = try allocator.alloc(T, rows * cols), .rows = rows, .cols = cols };
        }

        /// Initialize a matrix and fills it with zeros.
        pub fn initZero(rows: usize, cols: usize, allocator: std.mem.Allocator) !@This() {
            var A = try @This().init(rows, cols, allocator);
            A.zero();

            return A;
        }

        /// Free the memory allocated for the matrix.
        pub fn deinit(self: @This(), allocator: std.mem.Allocator) void {
            allocator.free(self.data);
        }

        /// Add another matrix to this matrix.
        pub fn add(self: *@This(), other: @This()) !void {
            if (self.rows != other.rows or self.cols != other.cols) {
                std.log.err("CANNOT ADD MATRICES OF DIFFERENT DIMENSIONS, FIRST WITH DIMENSIONS {d}x{d}, SECOND WITH DIMENSIONS {d}x{d}", .{ self.rows, self.cols, other.rows, other.cols });

                return error.ProgrammingError;
            }

            for (self.data, 0..) |*element, index| {
                element.* += other.data[index];
            }
        }

        /// Adds an empty row at the end of the matrix. The number of columns will be the same as the current number of columns.
        pub fn addRow(self: *@This(), allocator: std.mem.Allocator) !void {
            self.data = try allocator.realloc(self.data, (self.rows + 1) * self.cols);
            self.rows += 1;
        }

        /// Returns the matrix as a view to a real vector.
        pub fn asVector(self: @This()) RealVector(T) {
            return RealVector(T){ .data = self.data, .len = self.rows * self.cols };
        }

        /// Get the element at (i, j).
        pub fn at(self: @This(), i: usize, j: usize) T {
            return self.data[i * self.cols + j];
        }

        /// Clone the matrix.
        pub fn clone(self: @This(), allocator: std.mem.Allocator) !@This() {
            var B = try @This().init(self.rows, self.cols, allocator);

            try self.copyTo(&B);

            return B;
        }

        /// Returns the column as a view to a strided real vector.
        pub fn column(self: @This(), j: usize) StridedRealVector(T) {
            return StridedRealVector(T){ .data = self.data, .len = self.rows, .stride = self.cols, .zero = j };
        }

        /// Copy the contents of this matrix to another matrix.
        pub fn copyTo(self: @This(), other: *@This()) !void {
            if (self.rows != other.rows or self.cols != other.cols) {
                std.log.err("CANNOT COPY MATRICES OF DIFFERENT DIMENSIONS, FIRST WITH DIMENSIONS {d}x{d}, SECOND WITH DIMENSIONS {d}x{d}", .{ self.rows, self.cols, other.rows, other.cols });

                return error.ProgrammingError;
            }

            for (self.data, 0..) |element, index| {
                other.data[index] = element;
            }
        }

        /// Covariance matrix where each column is a variable and each row is an observation. The covariance is calculated using the unbiased estimator (dividing by n - 1).
        pub fn cov(self: @This(), allocator: std.mem.Allocator) !RealMatrix(T) {
            var result = try @This().init(self.cols, self.cols, allocator);
            errdefer result.deinit(allocator);

            for (0..self.cols) |i| for (i..self.cols) |j| {
                result.ptr(i, j).* = try self.column(i).covariance(self.column(j));
                result.ptr(j, i).* = result.at(i, j);
            };

            return result;
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

        /// Expand the matrix and keep the index structure.
        pub fn expand(self: *@This(), m: usize, n: usize, allocator: std.mem.Allocator) !void {
            if (m < self.rows or n < self.cols) {
                std.log.err("CANNOT EXPAND MATRIX TO SMALLER DIMENSIONS, CURRENT DIMENSIONS {d}x{d}, NEW DIMENSIONS {d}x{d}", .{ self.rows, self.cols, m, n });

                return error.ProgrammingError;
            }

            self.data = try allocator.realloc(self.data, m * n);

            var cursor = m * n;

            for (0..self.rows) |i| {
                const ii = self.rows - i - 1;

                const src_start = ii * self.cols;
                const dest_start = ii * n;

                @memmove(self.data[dest_start .. dest_start + self.cols], self.data[src_start .. src_start + self.cols]);

                @memset(self.data[dest_start + self.cols .. cursor], 0);

                cursor = dest_start;
            }

            self.rows = m;
            self.cols = n;
        }

        /// Fill the matrix with a given value.
        pub fn fill(self: *@This(), value: T) void {
            for (self.data) |*element| element.* = value;
        }

        /// Trace of a matrix.
        pub fn trace(self: @This()) T {
            var result: T = 0;

            for (0..@min(self.rows, self.cols)) |i| result += self.at(i, i);

            return result;
        }

        /// Transpose the matrix.
        pub fn transpose(self: *@This()) !void {
            if (!self.isSquare()) {
                std.log.err("CANNOT TRANSPOSE A NON-SQUARE MATRIX IN PLACE, CURRENT DIMENSIONS {d}x{d}", .{ self.rows, self.cols });

                return error.ProgrammingError;
            }

            for (0..self.rows) |i| for (i + 1..self.cols) |j| {
                const temp = self.at(i, j);

                self.ptr(i, j).* = self.at(j, i);
                self.ptr(j, i).* = temp;
            };
        }

        /// Maximum absolute value of the diagonal elements of the matrix.
        pub fn maxAbsDiagonal(self: @This()) T {
            var max: T = 0;

            for (0..@min(self.rows, self.cols)) |i| {
                const abs_value = @abs(self.at(i, i));

                if (abs_value > max) max = abs_value;
            }

            return max;
        }

        /// Minimum absolute value of the diagonal elements of the matrix.
        pub fn minAbsDiagonal(self: @This()) T {
            var min: T = std.math.inf(T);

            for (0..@min(self.rows, self.cols)) |i| {
                const abs_value = @abs(self.at(i, i));

                if (abs_value < min) min = abs_value;
            }

            return min;
        }

        /// Calculate the Frobenius norm of the matrix.
        pub fn frobeniusNorm(self: @This()) T {
            var sum: T = 0;

            for (self.data) |element| {
                sum += element * element;
            }

            return std.math.sqrt(sum);
        }

        /// Set the matrix to the identity matrix.
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

        /// Fill the matrix with random values using a given seed.
        pub fn randomize(self: *@This(), seed: u64) void {
            var split_mix = std.Random.SplitMix64.init(seed);
            var rng = std.Random.DefaultPrng.init(split_mix.next());

            var random = rng.random();

            for (self.data) |*e| e.* = random.float(T);
        }

        /// Reshape the matrix.
        pub fn reshape(self: *@This(), m: usize, n: usize) !void {
            if (m * n != self.rows * self.cols) {
                std.log.err("CANNOT RESHAPE MATRIX TO INCOMPATIBLE DIMENSIONS, CURRENT DIMENSIONS {d}x{d}, NEW DIMENSIONS {d}x{d}", .{ self.rows, self.cols, m, n });

                return error.ProgrammingError;
            }

            self.rows = m;
            self.cols = n;
        }

        /// Get a row as view to a a real vector.
        pub fn row(self: @This(), i: usize) RealVector(T) {
            return RealVector(T){ .data = self.data[i * self.cols .. (i + 1) * self.cols], .len = self.cols };
        }

        /// Shrink the matrix to a provided number of rows.
        pub fn shrinkRows(self: *@This(), m: usize, allocator: std.mem.Allocator) !void {
            if (m > self.rows) {
                std.log.err("CANNOT SHRINK MATRIX TO A NUMBER OF ROWS GREATER THAN THE CURRENT NUMBER OF ROWS, CURRENT NUMBER OF ROWS {d}, NEW NUMBER OF ROWS {d}", .{ self.rows, m });

                return error.ProgrammingError;
            }

            self.data = try allocator.realloc(self.data, m * self.cols);
            self.rows = m;
        }

        /// Shrink the matrix to a provided number of columns.
        pub fn shrinkCols(self: *@This(), n: usize, allocator: std.mem.Allocator) !void {
            if (n > self.cols) {
                std.log.err("CANNOT SHRINK MATRIX TO A NUMBER OF COLUMNS GREATER THAN THE CURRENT NUMBER OF COLUMNS, CURRENT NUMBER OF COLUMNS {d}, NEW NUMBER OF COLUMNS {d}", .{ self.cols, n });

                return error.ProgrammingError;
            }

            for (0..self.rows) |i| {
                const src_start = i * self.cols;
                const dest_start = i * n;

                std.mem.copyForwards(T, self.data[dest_start .. dest_start + n], self.data[src_start .. src_start + n]);
            }

            self.data = try allocator.realloc(self.data, self.rows * n);
            self.cols = n;
        }

        /// Symmetrize the matrix: A = 0.5 * (A + A^T)
        pub fn symmetrize(self: *@This()) !void {
            if (!self.isSquare()) {
                std.log.err("CANNOT SYMMETRIZE A NON-SQUARE MATRIX, CURRENT DIMENSIONS {d}x{d}", .{ self.rows, self.cols });

                return error.ProgrammingError;
            }

            for (0..self.rows) |i| for (i + 1..self.cols) |j| {
                const avg = 0.5 * (self.at(i, j) + self.at(j, i));

                self.ptr(i, j).* = avg;
                self.ptr(j, i).* = avg;
            };
        }

        /// Fills the matrix with zeros.
        pub fn zero(self: *@This()) void {
            self.fill(0);
        }
    };
}

test "init, deinit" {
    var A = try RealMatrix(f64).init(345, 753, std.testing.allocator);
    defer A.deinit(std.testing.allocator);

    try std.testing.expectEqual(A.rows, 345);
    try std.testing.expectEqual(A.cols, 753);
    try std.testing.expectEqual(A.data.len, 259785);
}
