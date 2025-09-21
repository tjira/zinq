//! Implementation of a ring buffer (circular buffer) in Zig.

const std = @import("std");

/// A ring buffer (circular buffer) implementation.
pub fn RingBuffer(comptime T: type) type {
    return struct {
        data: []T,
        len: usize,
        head: usize,

        allocator: std.mem.Allocator,

        /// Initialize the ring buffer.
        pub fn init(max_len: usize, allocator: std.mem.Allocator) !@This() {
            return @This(){
                .data = try allocator.alloc(T, max_len),
                .len = 0,
                .head = 0,
                .allocator = allocator,
            };
        }

        /// Deinitialize the ring buffer.
        pub fn deinit(self: @This()) void {
            self.allocator.free(self.data);
        }

        /// Push a new element to the ring buffer.
        pub fn append(self: *@This(), value: T) void {
            self.data[self.head] = value;
            self.head = (self.head + 1) % self.data.len;
            self.len = @min(self.len + 1, self.data.len);
        }

        /// Fill the ring buffer with a specific value.
        pub fn fill(self: *@This(), value: T) void {
            for (self.data) |*elem| elem.* = value;
        }

        /// Element access operator to get the i-th last element added (0 is the most recent).
        pub fn last(self: @This(), i: usize) T {
            return self.data[(self.head + self.data.len - i - 1) % self.data.len];
        }

        /// Fill the ring buffer with zeros.
        pub fn zero(self: *@This()) void {
            self.fill(0);
        }
    };
}
