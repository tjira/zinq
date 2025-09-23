//! A strided vector is a view into a vector with a non-1 stride.

const std = @import("std");

/// Real strided vector class.
pub fn StridedRealVector(comptime T: type) type {
    return struct {
        data: []T,
        len: usize,
        stride: usize,
        zero: usize,

        /// Get the element at the specified index as a value.
        pub fn at(self: @This(), i: usize) T {
            return self.data[self.zero + i * self.stride];
        }

        /// Binary search for a value in the strided array. If non found, return the index where it can be inserted so that the array remains sorted.
        pub fn bisectRight(self: @This(), value: T) usize {
            var low: usize = 0; var high = self.len; var mid = (low + high) / 2;

            while (low < high) : (mid = (low + high) / 2) {
                if (self.at(mid) < value) low = mid + 1 else high = mid;
            }

            return high;
        }

        /// Fill the strided array elements with a constant value.
        pub fn fill(self: @This(), value: T) void {
            for (0..self.len) |i| {
                self.ptr(i).* = value;
            }
        }

        /// Multiply the strided array elements with a constant value.
        pub fn muls(self: @This(), value: T) void {
            for (0..self.len) |i| {
                self.ptr(i).* *= value;
            }
        }

        /// Return the element at the specified index as a mutable reference.
        pub fn ptr(self: @This(), i: usize) *T {
            return &self.data[self.zero + i * self.stride];
        }

        /// Takes a slice of the strided vector and return a new strided vector.
        pub fn slice(self: @This(), start: usize, end: usize) @This() {
            return .{
                .data = self.data,
                .len = end - start,
                .stride = self.stride,
                .zero = self.zero + start * self.stride,
            };
        }

        /// Subtract a constant value from the strided array elements.
        pub fn subs(self: @This(), value: T) void {
            for (0..self.len) |i| {
                self.ptr(i).* = self.at(i).sub(value);
            }
        }
    };
}
