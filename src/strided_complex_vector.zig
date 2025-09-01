//! A strided vector is a view into a vector with a non-1 stride.

const std = @import("std");

const Complex = std.math.Complex;

/// Complex strided vector class.
pub fn StridedComplexVector(comptime T: type) type {
    return struct {
        data: []Complex(T),
        len: usize,
        stride: usize,
        zero: usize,

        /// Get the element at the specified index as a value.
        pub fn at(self: @This(), i: usize) Complex(T) {
            return self.data[self.zero + i * self.stride];
        }

        /// Fill the strided array elements with a constant value.
        pub fn fill(self: @This(), value: Complex(T)) void {
            for (0..self.len) |i| {
                self.ptr(i).* = value;
            }
        }

        /// Multiply the strided array elements with a constant value.
        pub fn muls(self: @This(), value: T) void {
            for (0..self.len) |i| {
                self.ptr(i).* = self.at(i).mul(value);
            }
        }

        /// Return the element at the specified index as a mutable reference.
        pub fn ptr(self: @This(), i: usize) *Complex(T) {
            return &self.data[self.zero + i * self.stride];
        }

        /// Subtract a constant value from the strided array elements.
        pub fn subs(self: @This(), value: T) void {
            for (0..self.len) |i| {
                self.ptr(i).* = self.at(i).sub(value);
            }
        }
    };
}
