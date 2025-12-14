//! 4th order tensor struct that uses real numbers (f64, i32, usize, ...) for its elements.

const std = @import("std");

const array_functions = @import("array_functions.zig");
const error_handling = @import("error_handling.zig");

const prod = array_functions.prod;
const throw = error_handling.throw;

/// Real tensor class. The tensor is stored in a flat array in row-major order.
pub fn RealTensor4(comptime T: type) type {
    return struct {
        data: []T,
        shape: [4]usize,

        /// Initialize a 4th order tensor with a given shape and specify an allocator. The function returns an error if the allocation fails.
        pub fn init(shape: [4]usize, allocator: std.mem.Allocator) !@This() {
            return @This(){
                .data = try allocator.alloc(T, shape[0] * shape[1] * shape[2] * shape[3]),
                .shape = shape
            };
        }

        /// Initialize a 4th order tensor and fills it with zeros.
        pub fn initZero(shape: [4]usize, allocator: std.mem.Allocator) !@This() {
            var A = try @This().init(shape, allocator); A.zero();

            return A;
        }

        /// Free the memory allocated for the tensor.
        pub fn deinit(self: @This(), allocator: std.mem.Allocator) void {
            allocator.free(self.data);
        }

        /// Get the element at (i, j, k, l).
        pub fn at(self: @This(), i: usize, j: usize, k: usize, l: usize) T {
            return self.data[((i * self.shape[1] + j) * self.shape[2] + k) * self.shape[3] + l];
        }

        /// Expands the tensor, keeping the correct indices.
        pub fn expand(self: *@This(), new_shape: [4]usize, allocator: std.mem.Allocator) !void {
            if (new_shape[0] < self.shape[0] or new_shape[1] < self.shape[1] or new_shape[2] < self.shape[2] or new_shape[3] < self.shape[3]) {
                return throw(void, "CAN'T EXPAND TENSOR FROM {d}x{d}x{d}x{d} TO {d}x{d}x{d}x{d}", .{self.shape[0], self.shape[1], self.shape[2], self.shape[3], new_shape[0], new_shape[1], new_shape[2], new_shape[3]});
            }

            self.data = try allocator.realloc(self.data, prod(usize, &new_shape));

            var cursor = prod(usize, &new_shape);

            for (0..self.shape[0]) |i| for (0..self.shape[1]) |j| for (0..self.shape[2]) |k| {

                const ii = self.shape[0] - i - 1; const jj = self.shape[1] - j - 1; const kk = self.shape[2] - k - 1;

                const src_idx  = ((ii * self.shape[1] + jj) * self.shape[2] + kk) * self.shape[3];
                const dest_idx = ((ii * new_shape[1] + jj) * new_shape[2] + kk) * new_shape[3];

                @memmove(self.data[dest_idx..dest_idx + self.shape[3]], self.data[src_idx..src_idx + self.shape[3]]);

                @memset(self.data[dest_idx + self.shape[3]..cursor], 0);

                cursor = dest_idx;
            };

            self.shape = new_shape;
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
