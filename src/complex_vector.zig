//! Complex vector class. The vector is stored in a flat array.

const std = @import("std");

const complex_matrix = @import("complex_matrix.zig");
const error_handling = @import("error_handling.zig");
const strided_complex_vector = @import("strided_complex_vector.zig");

const Complex = std.math.Complex;
const ComplexMatrix = complex_matrix.ComplexMatrix;
const StridedComplexVector = strided_complex_vector.StridedComplexVector;

const throw = error_handling.throw;

/// Complex vector class.
pub fn ComplexVector(comptime T: type) type {
    return struct {
        data: []Complex(T),
        len: usize,

        allocator: ?std.mem.Allocator,

        /// Initialize a vector with a given length and specify an allocator. The function returns an error if the allocation fails.
        pub fn init(len: usize, allocator: ?std.mem.Allocator) !@This() {
            return @This(){
                .data = try allocator.?.alloc(Complex(T), len),
                .len = len,
                .allocator = allocator
            };
        }

        /// Initialize a vector and fills it with zeros.
        pub fn initZero(len: usize, allocator: ?std.mem.Allocator) !@This() {
            var v = try @This().init(len, allocator); v.zero();

            return v;
        }

        /// Free the memory allocated for the vector.
        pub fn deinit(self: @This()) void {
            if (self.allocator) |allocator| allocator.free(self.data);
        }

        /// Create a matrix view of the vector with one column.
        pub fn asMatrix(self: @This()) ComplexMatrix(T) {
            return ComplexMatrix(T){
                .data = self.data,
                .rows = self.len,
                .cols = 1,
                .allocator = null,
            };
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

        /// Clone the vector.
        pub fn clone(self: @This()) !@This() {
            var v = try @This().init(self.len, self.allocator);

            for (self.data, 0..) |element, i| {
                v.data[i] = element;
            }

            return v;
        }

        /// Copy the contents of this matrix to another matrix.
        pub fn copyTo(self: @This(), other: *@This()) !void {
            if (self.len != other.len) {
                return throw(void, "CAN'T COPY A {d}-ELEMENT VECTOR TO A {d}-ELEMENT VECTOR", .{self.len, other.len});
            }

            for (self.data, 0..) |element, index| {
                other.data[index] = element;
            }
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

        /// Fill the vector with a given value.
        pub fn fill(self: *@This(), value: Complex(T)) void {
            for (self.data) |*element| element.* = value;
        }

        /// Get the pointer to the element at the specified index.
        pub fn ptr(self: @This(), i: usize) *Complex(T) {
            return &self.data[i];
        }

        /// Get a slice of the internal data in the form of a vector.
        pub fn slice(self: @This(), start: usize, end: usize) @This() {
            return @This(){
                .data = self.data[start..end],
                .len = end - start,
                .allocator = null,
            };
        }

        /// Fills the vector with zeros.
        pub fn zero(self: *@This()) void {
            self.fill(Complex(T).init(0, 0));
        }
    };
}
