//! 4th order tensor struct that uses real numbers (f64, i32, usize, ...) for its elements.

const std = @import("std");

/// Real tensor class. The tensor is stored in a flat array in row-major order.
pub fn RealTensor4(comptime T: type) type {
    return struct {
        data: []T,
        shape: [4]usize,

        allocator: ?std.mem.Allocator,

        /// Initialize a 4th order tensor with a given shape and specify an allocator. The function returns an error if the allocation fails.
        pub fn init(shape: [4]usize, allocator: ?std.mem.Allocator) !@This() {
            return @This(){
                .data = try allocator.?.alloc(T, shape[0] * shape[1] * shape[2] * shape[3]),
                .shape = shape,
                .allocator = allocator
            };
        }

        /// Initialize a 4th order tensor and fills it with zeros.
        pub fn initZero(shape: [4]usize, allocator: ?std.mem.Allocator) !@This() {
            var A = try @This().init(shape, allocator); A.zero();

            return A;
        }

        /// Free the memory allocated for the tensor.
        pub fn deinit(self: @This()) void {
            if (self.allocator) |allocator| allocator.free(self.data);
        }

        /// Get the element at (i, j, k, l).
        pub fn at(self: @This(), i: usize, j: usize, k: usize, l: usize) T {
            return self.data[((i * self.shape[1] + j) * self.shape[2] + k) * self.shape[3] + l];
        }

        /// Fills the tensor with constants.
        pub fn fill(self: *@This(), value: T) void {
            for (self.data) |*element| element.* = value;
        }

        /// Get a pointer to the element at (i, j, k, l).
        pub fn ptr(self: *@This(), i: usize, j: usize, k: usize, l: usize) *T {
            return &self.data[((i * self.shape[1] + j) * self.shape[2] + k) * self.shape[3] + l];
        }

        /// Set all elements to zero.
        pub fn zero(self: *@This()) void {
            self.fill(0);
        }
    };
}
