//! Real vector class. The vector is stored in a flat array.

const std = @import("std");

/// Real vector class.
pub fn RealVector(comptime T: type) type {
    return struct {
        data: []T,
        len: usize,
        allocator: std.mem.Allocator,

        /// Initialize a vector with a given length and specify an allocator. The function returns an error if the allocation fails.
        pub fn init(len: usize, allocator: std.mem.Allocator) !@This() {
            return @This(){
                .data = try allocator.alloc(T, len),
                .len = len,
                .allocator = allocator
            };
        }

        /// Initialize a vector and fills it with zeros.
        pub fn initZero(len: usize, allocator: std.mem.Allocator) !@This() {
            var v = try @This().init(len, allocator); v.zero();

            return v;
        }

        /// Free the memory allocated for the vector.
        pub fn deinit(self: @This()) void {
            self.allocator.free(self.data);
        }

        /// Get the element at index i.
        pub fn at(self: @This(), i: usize) T {
            return self.data[i];
        }

        /// Equality operator for vectors.
        pub fn eq(self: @This(), other: @This()) bool {
            if (self.len != other.len) return false;

            for (self.data, 0..) |element, index| {
                if (element != other.data[index]) return false;
            }

            return true;
        }

        /// Fill the vector with a given value.
        pub fn fill(self: *@This(), value: T) void {
            @memset(self.data, value);
        }

        /// Get the pointer to the element at index i.
        pub fn ptr(self: *@This(), i: usize) *T {
            return &self.data[i];
        }

        /// Fills the vector with zeros.
        pub fn zero(self: *@This()) void {
            self.fill(0);
        }
    };
}
