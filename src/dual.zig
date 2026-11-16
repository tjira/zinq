//! Implements dual numbers to enable automatic differentiation for scientific and numerical computations.

const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

/// Returns a dual number type for computing value and first derivative of a scalar function.
pub fn ScalarDual(comptime T: type) type {
    return struct {
        val: T,
        der: T,

        /// Initializes a dual number with a real value and its corresponding first-order derivative.
        pub fn init(val: T, der: T) @This() {
            return .{ .val = val, .der = der };
        }

        /// Computes the absolute value and its derivative using the sign function derivative step.
        pub fn abs(self: @This()) @This() {
            return .{ .val = @abs(self.val), .der = if (self.val >= 0) self.der else -self.der };
        }

        /// Computes the sum of two dual numbers using the linearity of differentiation: d(u+v) = du + dv.
        pub fn add(self: @This(), other: @This()) @This() {
            return .{ .val = self.val + other.val, .der = self.der + other.der };
        }

        /// Adds a constant scalar to a dual number, leaving the derivative term unchanged.
        pub fn adds(self: @This(), scalar: T) @This() {
            return .{ .val = self.val + scalar, .der = self.der };
        }

        /// Divides two dual numbers using the quotient rule: d(u/v) = (du*v - u*dv) / v^2.
        pub fn div(self: @This(), other: @This()) @This() {
            const der = (self.der * other.val - self.val * other.der) / (other.val * other.val);

            return .{ .val = self.val / other.val, .der = der };
        }

        /// Divides a dual number by a constant scalar, scaling both value and derivative.
        pub fn divs(self: @This(), scalar: T) @This() {
            return .{ .val = self.val / scalar, .der = self.der / scalar };
        }

        /// Computes the exponential of a dual number using the chain rule: d(e^u) = e^u * du.
        pub fn exp(self: @This()) @This() {
            const val = std.math.exp(self.val);

            return .{ .val = val, .der = val * self.der };
        }

        /// Multiplies two dual numbers using the product rule: d(uv) = u*dv + v*du.
        pub fn mul(self: @This(), other: @This()) @This() {
            return .{ .val = self.val * other.val, .der = self.val * other.der + other.val * self.der };
        }

        /// Multiplies a dual number by a constant scalar, scaling both value and derivative.
        pub fn muls(self: @This(), scalar: T) @This() {
            return .{ .val = self.val * scalar, .der = self.der * scalar };
        }

        /// Computes the difference of two dual numbers: d(u-v) = du - dv.
        pub fn sub(self: @This(), other: @This()) @This() {
            return .{ .val = self.val - other.val, .der = self.der - other.der };
        }

        /// Subtracts a constant scalar from a dual number, leaving the derivative term unchanged.
        pub fn subs(self: @This(), scalar: T) @This() {
            return .{ .val = self.val - scalar, .der = self.der };
        }
    };
}

/// Holds separate contiguous value and derivative vectors enabling vectorization in automatic differentiation.
pub fn DualVector(comptime T: type) type {
    return struct {
        val: Vector(T),
        der: Vector(T),

        /// Creates a dual vector from existing value and derivative vectors.
        pub fn fromVectors(val: Vector(T), der: Vector(T)) @This() {
            std.debug.assert(val.length() == der.length());

            return .{ .val = val, .der = der };
        }

        /// Allocates value and derivative storage for a dual vector.
        pub fn init(len: usize, gpa: Allocator) !@This() {
            var val = try Vector(T).init(len, gpa);
            errdefer val.deinit(gpa);

            var der = try Vector(T).init(len, gpa);
            errdefer der.deinit(gpa);

            return .{ .val = val, .der = der };
        }

        /// Allocates zero-initialized value and derivative storage for a dual vector.
        pub fn initZero(len: usize, gpa: Allocator) !@This() {
            var val = try Vector(T).initZero(len, gpa);
            errdefer val.deinit(gpa);

            var der = try Vector(T).initZero(len, gpa);
            errdefer der.deinit(gpa);

            return .{ .val = val, .der = der };
        }

        /// Deallocates memory associated with value and derivative vectors.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.val.deinit(gpa);
            self.der.deinit(gpa);
        }

        /// Returns the number of elements in the vector.
        pub fn length(self: @This()) usize {
            return self.val.length();
        }

        /// Returns the dual scalar element at the given linear index.
        pub fn at(self: @This(), i: usize) ScalarDual(T) {
            return ScalarDual(T).init(self.val.at(i), self.der.at(i));
        }

        /// Assigns dual scalar element at linear coordinate i.
        pub fn set(self: *@This(), i: usize, val: ScalarDual(T)) void {
            self.val.ptr(i).* = val.val;
            self.der.ptr(i).* = val.der;
        }

        /// Assigns real value and first derivative at the specified index.
        pub fn setValDer(self: *@This(), i: usize, val: T, der: T) void {
            self.val.ptr(i).* = val;
            self.der.ptr(i).* = der;
        }

        /// Zeros all value and derivative elements of the dual vector.
        pub fn zero(self: *@This()) void {
            self.val.zero();
            self.der.zero();
        }
    };
}

/// Holds separate contiguous value and derivative matrices enabling BLAS calls in automatic differentiation.
pub fn DualMatrix(comptime T: type) type {
    return struct {
        val: Matrix(T),
        der: Matrix(T),

        shape: [2]usize,

        /// Creates a dual matrix from existing value and derivative matrices.
        pub fn fromMatrices(val: Matrix(T), der: Matrix(T)) @This() {
            std.debug.assert(val.nrow() == der.nrow());
            std.debug.assert(val.ncol() == der.ncol());

            return .{ .val = val, .der = der, .shape = val.shape };
        }

        /// Allocates value and derivative storage for a dual matrix.
        pub fn init(rows: usize, cols: usize, gpa: Allocator) !@This() {
            var val = try Matrix(T).init(rows, cols, gpa);
            errdefer val.deinit(gpa);

            var der = try Matrix(T).init(rows, cols, gpa);
            errdefer der.deinit(gpa);

            return .{ .val = val, .der = der, .shape = .{ rows, cols } };
        }

        /// Allocates zero-initialized value and derivative storage for a dual matrix.
        pub fn initZero(rows: usize, cols: usize, gpa: Allocator) !@This() {
            var val = try Matrix(T).initZero(rows, cols, gpa);
            errdefer val.deinit(gpa);

            var der = try Matrix(T).initZero(rows, cols, gpa);
            errdefer der.deinit(gpa);

            return .{ .val = val, .der = der, .shape = .{ rows, cols } };
        }

        /// Deallocates memory associated with value and derivative matrices.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.val.deinit(gpa);
            self.der.deinit(gpa);
        }

        /// Returns the number of rows in the matrix.
        pub fn nrow(self: @This()) usize {
            return self.val.nrow();
        }

        /// Returns the number of columns in the matrix.
        pub fn ncol(self: @This()) usize {
            return self.val.ncol();
        }

        /// Returns the dual scalar element at the specified row and column.
        pub fn at(self: @This(), i: usize, j: usize) ScalarDual(T) {
            return ScalarDual(T).init(self.val.at(i, j), self.der.at(i, j));
        }

        /// Assigns dual scalar element at row i and column j.
        pub fn set(self: *@This(), i: usize, j: usize, val: ScalarDual(T)) void {
            self.val.ptr(i, j).* = val.val;
            self.der.ptr(i, j).* = val.der;
        }

        /// Assigns real value and first derivative at the specified row and column.
        pub fn setValDer(self: *@This(), i: usize, j: usize, val: T, der: T) void {
            self.val.ptr(i, j).* = val;
            self.der.ptr(i, j).* = der;
        }

        /// Zeros all value and derivative elements of the dual matrix.
        pub fn zero(self: *@This()) void {
            self.val.zero();
            self.der.zero();
        }

        /// Reinterprets dual matrix as flat dual vector.
        pub fn asVector(self: @This()) DualVector(T) {
            return .{ .val = self.val.asVector(), .der = self.der.asVector() };
        }
    };
}

/// Holds separate contiguous value and derivative rank-N tensors enabling BLAS contractions in dual arithmetic.
pub fn DualTensor(comptime T: type, comptime N: usize) type {
    return struct {
        val: Tensor(T, N),
        der: Tensor(T, N),

        shape: [N]usize,

        /// Creates a rank-N dual tensor from existing value and derivative tensors.
        pub fn fromTensors(val: Tensor(T, N), der: Tensor(T, N)) @This() {
            std.debug.assert(std.mem.eql(usize, &val.shape, &der.shape));

            return .{ .val = val, .der = der, .shape = val.shape };
        }

        /// Allocates value and derivative storage for a rank-N dual tensor.
        pub fn init(shape: [N]usize, gpa: Allocator) !@This() {
            var val = try Tensor(T, N).init(shape, gpa);
            errdefer val.deinit(gpa);

            var der = try Tensor(T, N).init(shape, gpa);
            errdefer der.deinit(gpa);

            return .{ .val = val, .der = der, .shape = shape };
        }

        /// Allocates zero-initialized value and derivative storage for a rank-N dual tensor.
        pub fn initZero(shape: [N]usize, gpa: Allocator) !@This() {
            var val = try Tensor(T, N).initZero(shape, gpa);
            errdefer val.deinit(gpa);

            var der = try Tensor(T, N).initZero(shape, gpa);
            errdefer der.deinit(gpa);

            return .{ .val = val, .der = der, .shape = shape };
        }

        /// Deallocates memory associated with value and derivative tensors.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.val.deinit(gpa);
            self.der.deinit(gpa);
        }

        /// Returns the dual scalar element at the given multi-dimensional indices.
        pub fn at(self: @This(), indices: [N]usize) ScalarDual(T) {
            return ScalarDual(T).init(self.val.at(indices), self.der.at(indices));
        }

        /// Assigns dual scalar element at multi-dimensional coordinate indices.
        pub fn set(self: *@This(), indices: [N]usize, val: ScalarDual(T)) void {
            self.val.ptr(indices).* = val.val;
            self.der.ptr(indices).* = val.der;
        }

        /// Assigns real value and first derivative at the specified tensor indices.
        pub fn setValDer(self: *@This(), indices: [N]usize, val: T, der: T) void {
            self.val.ptr(indices).* = val;
            self.der.ptr(indices).* = der;
        }

        /// Zeros all value and derivative elements of the dual tensor.
        pub fn zero(self: *@This()) void {
            self.val.zero();
            self.der.zero();
        }

        /// Reshapes dual tensor as a 2D dual matrix with specified dimensions.
        pub fn matrixView(self: @This(), rows: usize, cols: usize) DualMatrix(T) {
            return .{ .val = self.val.matrixView(rows, cols), .der = self.der.matrixView(rows, cols), .shape = .{ rows, cols } };
        }

        /// Slices the first tensor dimension, returning a rank-(N-1) dual tensor.
        pub fn slice(self: @This(), i: usize) DualTensor(T, N - 1) {
            const s_val = self.val.slice(i);
            const s_der = self.der.slice(i);

            return .{ .val = s_val, .der = s_der, .shape = s_val.shape };
        }

        /// Converts a 2D dual tensor to a dual matrix.
        pub fn asMatrix(self: @This()) DualMatrix(T) {
            return .{ .val = self.val.asMatrix(), .der = self.der.asMatrix(), .shape = self.val.shape };
        }
    };
}
