//! Matrix struct that uses complex numbers for its elements.

const std = @import("std");

const complex_vector = @import("complex_vector.zig");
const strided_complex_vector = @import("strided_complex_vector.zig");

const Complex = std.math.complex.Complex;
const ComplexVector = complex_vector.ComplexVector;
const StridedComplexVector = strided_complex_vector.StridedComplexVector;

/// Complex matrix class. The matrix is stored in a flat array in row-major order.
pub fn ComplexMatrix(comptime T: type) type {
    return struct {
        data: []Complex(T),
        rows: usize,
        cols: usize,

        /// Initialize a matrix with a given number of rows and columns and specify an allocator. The function returns an error if the allocation fails.
        pub fn init(rows: usize, cols: usize, allocator: std.mem.Allocator) !@This() {
            return @This(){ .data = try allocator.alloc(Complex(T), rows * cols), .rows = rows, .cols = cols };
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
                std.log.err("CANNOT ADD MATRICES OF DIFFERENT SIZES, FIRST WITH SIZE {d}x{d}, SECOND WITH SIZE {d}x{d}", .{ self.rows, self.cols, other.rows, other.cols });

                return error.ProgrammingError;
            }

            for (self.data, 0..) |*element, index| {
                element.* = element.add(other.data[index]);
            }
        }

        /// Returns the matrix as a view to a complex vector.
        pub fn asVector(self: @This()) ComplexVector(T) {
            return ComplexVector(T){ .data = self.data, .len = self.rows * self.cols };
        }

        /// Get the element at (i, j).
        pub fn at(self: @This(), i: usize, j: usize) Complex(T) {
            return self.data[i * self.cols + j];
        }

        /// Returns the column as a view to a strided real vector.
        pub fn column(self: @This(), j: usize) StridedComplexVector(T) {
            return StridedComplexVector(T){ .data = self.data, .len = self.rows, .stride = self.cols, .zero = j };
        }

        /// Copy the contents of this matrix to another matrix.
        pub fn copyTo(self: @This(), other: *@This()) !void {
            if (self.rows != other.rows or self.cols != other.cols) {
                std.log.err("CANNOT COPY MATRICES OF DIFFERENT SIZES, FIRST WITH SIZE {d}x{d}, SECOND WITH SIZE {d}x{d}", .{ self.rows, self.cols, other.rows, other.cols });

                return error.ProgrammingError;
            }

            for (self.data, 0..) |element, index| {
                other.data[index] = element;
            }
        }

        /// Divide the matrix by a scalar value.
        pub fn divs(self: *@This(), value: Complex(T)) void {
            for (self.data) |*element| {
                element.* = element.div(value);
            }
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

        /// Calculate the Frobenius norm of the matrix.
        pub fn frobeniusNorm(self: @This()) T {
            var sum: T = 0;

            for (self.data) |element| {
                sum += element.magnitude() * element.magnitude();
            }

            return std.math.sqrt(sum);
        }

        /// Maximum absolute value on the diagonal of the matrix.
        pub fn maxAbsDiagonal(self: @This()) T {
            var max: T = 0;

            for (0..@min(self.rows, self.cols)) |i| {
                const abs_value = self.at(i, i).magnitude();

                if (abs_value > max) max = abs_value;
            }

            return max;
        }

        /// Minimum absolute value on the diagonal of the matrix.
        pub fn minAbsDiagonal(self: @This()) T {
            var min: T = std.math.inf(T);

            for (0..@min(self.rows, self.cols)) |i| {
                const abs_value = self.at(i, i).magnitude();

                if (abs_value < min) min = abs_value;
            }

            return min;
        }

        /// Set the matrix to the identity matrix.
        pub fn identity(self: *@This()) void {
            self.fill(Complex(T).init(0, 0));

            for (0..@min(self.rows, self.cols)) |i| {
                self.ptr(i, i).* = Complex(T).init(1, 0);
            }
        }

        /// Is the metrix hermitian.
        pub fn isHermitian(self: @This(), tol: T) bool {
            if (!self.isSquare()) return false;

            for (0..self.rows) |i| for (0..self.cols) |j| {
                if (@abs(self.at(i, j).re - self.at(j, i).re) > tol) return false;
                if (@abs(self.at(i, j).im + self.at(j, i).im) > tol) return false;
            };

            return true;
        }

        /// Square checker.
        pub fn isSquare(self: @This()) bool {
            return self.rows == self.cols;
        }

        /// Multiply the matrix by a scalar value.
        pub fn muls(self: *@This(), value: Complex(T)) void {
            for (self.data) |*element| {
                element.* = element.mul(value);
            }
        }

        /// Calculate the Frobenius norm of the off-diagonal elements of the matrix.
        pub fn offDiagonalFrobeniusNorm(self: @This()) T {
            var sum: T = 0;

            for (0..self.rows) |i| for (0..self.cols) |j| if (i != j) {
                sum += self.at(i, j).magnitude() * self.at(i, j).magnitude();
            };

            return std.math.sqrt(sum);
        }

        /// Get the pointer to the element at (i, j).
        pub fn ptr(self: *@This(), i: usize, j: usize) *Complex(T) {
            return &self.data[i * self.cols + j];
        }

        /// Returns the row of the matrix as a complex vector.
        pub fn row(self: @This(), i: usize) ComplexVector(T) {
            return ComplexVector(T){ .data = self.data[i * self.cols .. (i + 1) * self.cols], .len = self.cols };
        }

        /// Returns the trace of the matrix.
        pub fn trace(self: @This()) Complex(T) {
            var result = Complex(T).init(0, 0);

            for (0..@min(self.rows, self.cols)) |i| {
                result = result.add(self.at(i, i));
            }

            return result;
        }

        /// Fills the matrix with zeros.
        pub fn zero(self: *@This()) void {
            self.fill(Complex(T).init(0, 0));
        }
    };
}

test "init, deinit" {
    var A = try ComplexMatrix(f64).init(345, 753, std.testing.allocator);
    defer A.deinit(std.testing.allocator);

    try std.testing.expectEqual(A.rows, 345);
    try std.testing.expectEqual(A.cols, 753);
    try std.testing.expectEqual(A.data.len, 259785);
}
