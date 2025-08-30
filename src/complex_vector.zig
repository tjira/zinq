//! Complex vector class. The vector is stored in a flat array.

const std = @import("std");

const strided_complex_vector = @import("strided_complex_vector.zig");

const Complex = std.math.Complex;
const StridedComplexVector = strided_complex_vector.StridedComplexVector;

/// Complex vector class.
pub fn ComplexVector(comptime T: type) type {
    return struct {
        data: []Complex(T),
        len: usize,
        allocator: std.mem.Allocator,

        /// Initialize a vector with a given length and specify an allocator. The function returns an error if the allocation fails.
        pub fn init(len: usize, allocator: std.mem.Allocator) !@This() {
            return @This(){
                .data = try allocator.alloc(Complex(T), len),
                .len = len,
                .allocator = allocator
            };
        }

        /// Free the memory allocated for the vector.
        pub fn deinit(self: @This()) void {
            self.allocator.free(self.data);
        }

        /// Create a strided view of the vector.
        pub fn asStrided(self: @This()) StridedComplexVector(T) {
            return StridedComplexVector(T){
                .data = self.data,
                .len = self.data.len,
                .stride = 1,
                .zero = 0,
            };
        }
        
        /// Get the element at the specified index as a value.
        pub fn at(self: @This(), i: usize) Complex(T) {
            return self.data[i];
        }

        /// Check for equality with another vector, given a tolerance.
        pub fn eq(self: @This(), other: @This(), tol: T) bool {
            if (self.len != other.len) return false;

            for (self.data, 0..) |element, i| {
                if (@abs(element.re - other.data[i].re) > tol) return false;
                if (@abs(element.im - other.data[i].im) > tol) return false;
            }

            return true;
        }

        /// Get the pointer to the element at the specified index.
        pub fn ptr(self: @This(), i: usize) *Complex(T) {
            return &self.data[i];
        }
    };
}
