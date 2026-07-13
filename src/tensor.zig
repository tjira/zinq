//! Multi-dimensional array (tensor, matrix, vector) containers for numerical physics and linear algebra.

const std = @import("std");

const Value = @import("value.zig").Value;

/// Returns a 2D matrix type representing linear operators or grids in coordinate space.
pub fn Matrix(comptime T: type) type {
    return struct {
        data: []T,
        shape: [2]usize,

        /// Allocates memory for a matrix of size rows * cols.
        pub fn init(rows: usize, cols: usize, gpa: std.mem.Allocator) !@This() {
            return .{ .data = try gpa.alloc(T, rows * cols), .shape = .{ rows, cols } };
        }

        /// Wraps an existing slice into a 2D matrix view of size rows * cols.
        pub fn fromSlice(rows: usize, cols: usize, data: []T) @This() {
            std.debug.assert(data.len == rows * cols);

            return .{ .data = data, .shape = .{ rows, cols } };
        }

        /// Frees the allocated memory of the matrix.
        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            gpa.free(self.data);
        }

        /// Flattens the 2D matrix into a 1D vector representation.
        pub fn asVector(self: @This()) Vector(T) {
            return .{ .data = self.data, .shape = .{self.data.len} };
        }

        /// Returns the matrix element at row i and column j.
        pub fn at(self: @This(), i: usize, j: usize) T {
            std.debug.assert(i < self.shape[0]);
            std.debug.assert(j < self.shape[1]);

            return self.data[i * self.shape[1] + j];
        }

        /// Creates a deep copy of the matrix using the provided allocator.
        pub fn clone(self: @This(), gpa: std.mem.Allocator) !@This() {
            var A = try @This().init(self.shape[0], self.shape[1], gpa);

            for (0..self.data.len) |i| {
                A.data[i] = self.data[i];
            }

            return A;
        }

        /// Divides all elements of the matrix by a scalar in-place.
        pub fn divs(self: *@This(), scalar: T) void {
            for (0..self.data.len) |i| {
                self.data[i] = Value(T).init(self.data[i]).div(Value(T).init(scalar)).val;
            }
        }

        /// Fills the matrix with a uniform scalar value.
        pub fn fill(self: *@This(), scalar: T) void {
            for (0..self.data.len) |i| {
                self.data[i] = scalar;
            }
        }

        /// Allocates a matrix and initializes all elements to zero.
        pub fn initZero(rows: usize, cols: usize, gpa: std.mem.Allocator) !@This() {
            var A = try @This().init(rows, cols, gpa);

            A.zero();

            return A;
        }

        /// Returns the maximum absolute value (infinity norm candidate) of the matrix elements.
        pub fn max(self: @This()) T {
            var max_val: T = 0;

            for (self.data) |x| {
                const abs_val = Value(T).init(x).abs().val;

                if (abs_val > max_val) {
                    max_val = abs_val;
                }
            }

            return max_val;
        }

        /// Returns the number of columns in the matrix.
        pub fn ncol(self: @This()) usize {
            return self.shape[1];
        }

        /// Returns the number of rows in the matrix.
        pub fn nrow(self: @This()) usize {
            return self.shape[0];
        }

        /// Returns a pointer to the matrix element at row i and column j.
        pub fn ptr(self: *@This(), i: usize, j: usize) *T {
            std.debug.assert(i < self.shape[0]);
            std.debug.assert(j < self.shape[1]);

            return &self.data[i * self.shape[1] + j];
        }

        /// Returns a 1D vector view of the specified row.
        pub fn row(self: @This(), i: usize) Vector(T) {
            std.debug.assert(i < self.shape[0]);

            return .{ .data = self.data[i * self.shape[1] .. (i + 1) * self.shape[1]], .shape = .{self.shape[1]} };
        }

        /// Returns a slice view of the specified row.
        pub fn rowSlice(self: @This(), i: usize) []T {
            std.debug.assert(i < self.shape[0]);

            return self.data[i * self.shape[1] .. (i + 1) * self.shape[1]];
        }

        /// Computes the root-mean-square value of the matrix elements.
        pub fn rms(self: @This()) T {
            var sum_sq: T = 0;

            for (self.data) |x| {
                const abs_val = Value(T).init(x).abs().val;

                sum_sq += abs_val * abs_val;
            }

            return @sqrt(sum_sq / @as(T, @floatFromInt(self.data.len)));
        }

        /// Returns a submatrix view consisting of the first n rows.
        pub fn takeRows(self: @This(), n: usize) @This() {
            std.debug.assert(n <= self.shape[0]);

            return .{ .data = self.data[0 .. n * self.shape[1]], .shape = .{ n, self.shape[1] } };
        }

        /// Sets all elements of the matrix to zero.
        pub fn zero(self: *@This()) void {
            self.fill(std.mem.zeroes(T));
        }
    };
}

/// Returns a 1D vector type representing physical coordinates or state vectors.
pub fn Vector(comptime T: type) type {
    return struct {
        data: []T,
        shape: [1]usize,

        /// Allocates memory for a vector of the specified size.
        pub fn init(size: usize, gpa: std.mem.Allocator) !@This() {
            return .{ .data = try gpa.alloc(T, size), .shape = .{size} };
        }

        /// Wraps an existing slice into a 1D vector view.
        pub fn fromSlice(data: []T) @This() {
            return .{ .data = data, .shape = .{data.len} };
        }

        /// Frees the allocated memory of the vector.
        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            gpa.free(self.data);
        }

        /// Promotes the 1D vector to a 2D column matrix (N x 1).
        pub fn asMatrix(self: @This()) Matrix(T) {
            return .{ .data = self.data, .shape = .{ self.shape[0], 1 } };
        }

        /// Returns the vector element at index i.
        pub fn at(self: @This(), i: usize) T {
            std.debug.assert(i < self.shape[0]);

            return self.data[i];
        }

        /// Divides all elements of the vector by a scalar in-place.
        pub fn divs(self: *@This(), scalar: T) void {
            for (0..self.data.len) |i| {
                self.data[i] = Value(T).init(self.data[i]).div(Value(T).init(scalar)).val;
            }
        }

        /// Allocates a vector and initializes all elements to zero.
        pub fn initZero(size: usize, gpa: std.mem.Allocator) !@This() {
            var v = try @This().init(size, gpa);

            v.zero();

            return v;
        }

        /// Returns the number of elements in the vector.
        pub fn length(self: @This()) usize {
            return self.shape[0];
        }

        /// Multiplies all elements of the vector by a scalar in-place.
        pub fn muls(self: *@This(), scalar: T) void {
            for (0..self.data.len) |i| {
                self.data[i] = Value(T).init(self.data[i]).mul(Value(T).init(scalar)).val;
            }
        }

        /// Returns a pointer to the vector element at index i.
        pub fn ptr(self: *@This(), i: usize) *T {
            std.debug.assert(i < self.shape[0]);

            return &self.data[i];
        }

        /// Returns a subvector view consisting of the first n elements.
        pub fn takeRows(self: @This(), n: usize) @This() {
            std.debug.assert(n <= self.shape[0]);

            return .{ .data = self.data[0..n], .shape = .{n} };
        }

        /// Sets all elements of the vector to zero.
        pub fn zero(self: *@This()) void {
            for (0..self.data.len) |i| {
                self.data[i] = std.mem.zeroes(T);
            }
        }
    };
}

/// Returns a multi-dimensional array type of rank N for physics tensors.
pub fn Tensor(comptime T: type, comptime N: usize) type {
    return struct {
        data: []T,
        shape: [N]usize,

        /// Allocates memory for a tensor with the given multi-index shape.
        pub fn init(shape: [N]usize, gpa: std.mem.Allocator) !@This() {
            var size: usize = 1;

            inline for (0..N) |i| {
                size *= shape[i];
            }

            return .{ .data = try gpa.alloc(T, size), .shape = shape };
        }

        /// Frees the allocated memory of the tensor.
        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            gpa.free(self.data);
        }

        /// Reshapes the N-dimensional tensor into a 2D matrix representation.
        pub fn asMatrix(self: @This()) Matrix(T) {
            var rows: usize = 1;
            var cols: usize = 1;

            for (0..(N + 1) / 2) |i| {
                rows *= self.shape[i];
            }

            for ((N + 1) / 2..N) |i| {
                cols *= self.shape[i];
            }

            return .{ .data = self.data, .shape = .{ rows, cols } };
        }

        /// Returns the tensor element at the specified multi-index coordinate.
        pub fn at(self: @This(), indx: [N]usize) T {
            var idx: usize = 0;
            var str: usize = 1;

            inline for (0..N) |k| {
                const j = N - 1 - k;

                std.debug.assert(indx[j] < self.shape[j]);

                idx += indx[j] * str;
                str *= self.shape[j];
            }

            return self.data[idx];
        }

        /// Allocates a tensor and initializes all elements to zero.
        pub fn initZero(shape: [N]usize, gpa: std.mem.Allocator) !@This() {
            var U = try @This().init(shape, gpa);

            U.zero();

            return U;
        }

        /// Returns a pointer to the tensor element at the specified multi-index coordinate.
        pub fn ptr(self: *@This(), indx: [N]usize) *T {
            var idx: usize = 0;
            var str: usize = 1;

            inline for (0..N) |k| {
                const j = N - 1 - k;

                std.debug.assert(indx[j] < self.shape[j]);

                idx += indx[j] * str;
                str *= self.shape[j];
            }

            return &self.data[idx];
        }

        /// Sets all elements of the tensor to zero.
        pub fn zero(self: *@This()) void {
            for (0..self.data.len) |i| {
                self.data[i] = std.mem.zeroes(T);
            }
        }
    };
}
