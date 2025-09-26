//! Real vector class. The vector is stored in a flat array.

const std = @import("std");

const error_handling = @import("error_handling.zig");
const real_matrix = @import("real_matrix.zig");

const RealMatrix = real_matrix.RealMatrix;

const throw = error_handling.throw;

/// Real vector class.
pub fn RealVector(comptime T: type) type {
    return struct {
        data: []T,
        len: usize,

        allocator: ?std.mem.Allocator,

        /// Initialize a vector with a given length and specify an allocator. The function returns an error if the allocation fails.
        pub fn init(len: usize, allocator: ?std.mem.Allocator) !@This() {
            return @This(){
                .data = try allocator.?.alloc(T, len),
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

        /// Add another vector to this vector.
        pub fn add(self: *@This(), other: @This()) void {
            std.debug.assert(self.len == other.len);

            for (self.data, 0..) |*element, index| {
                element.* += other.data[index];
            }
        }

        /// Return a view of the internal data as a column matrix.
        pub fn asMatrix(self: @This()) RealMatrix(T) {
            return RealMatrix(T){
                .data = self.data,
                .rows = self.len,
                .cols = 1,
                .allocator = null,
            };
        }

        /// Get the element at index i.
        pub fn at(self: @This(), i: usize) T {
            return self.data[i];
        }

        /// Clone the vector.
        pub fn clone(self: @This()) !@This() {
            var v = try @This().init(self.len, self.allocator);

            try self.copyTo(&v);

            return v;
        }

        /// Copy the contents of this vector to another vector.
        pub fn copyTo(self: @This(), other: *@This()) !void {
            if (self.len != other.len) return throw(void, "CAN'T COPY A VECTOR OF LENGTH {d} TO A VECTOR OF LENGTH {d}", .{self.len, other.len});

            for (self.data, 0..) |element, index| {
                other.data[index] = element;
            }
        }

        /// Divide the vector by a scalar value.
        pub fn divs(self: *@This(), value: T) void {
            for (self.data) |*element| element.* /= value;
        }

        /// Equality operator for vectors.
        pub fn eq(self: @This(), other: @This(), tol: T) bool {
            if (self.len != other.len) return false;

            for (self.data, 0..) |element, index| {
                if (@abs(element - other.data[index]) > tol) return false;
            }

            return true;
        }

        /// Fill the vector with a given value.
        pub fn fill(self: *@This(), value: T) void {
            for (self.data) |*element| element.* = value;
        }

        /// Calculate the mean of all elements in the vector.
        pub fn mean(self: @This()) T {
            return self.sum() / @as(T, @floatFromInt(self.len));
        }

        /// Multiply the vector by a scalar value.
        pub fn muls(self: *@This(), value: T) void {
            for (self.data) |*element| element.* *= value;
        }

        /// Calculate the norm of the vector.
        pub fn norm(self: @This()) T {
            var total: T = 0;

            for (self.data) |element| total += element * element;

            return std.math.sqrt(total);
        }

        /// Get the pointer to the element at index i.
        pub fn ptr(self: *@This(), i: usize) *T {
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

        /// Subtract another vector from this vector.
        pub fn sub(self: *@This(), other: @This()) !void {
            if (self.len != other.len) return throw(void, "CAN'T SUBTRACT A VECTOR OF LENGTH {d} FROM A VECTOR OF LENGTH {d}", .{other.len, self.len});

            for (self.data, 0..) |*element, index| {
                element.* -= other.data[index];
            }
        }

        /// Calculates the sum of all elements in the vector.
        pub fn sum(self: @This()) T {
            var total: T = 0;

            for (self.data) |element| total += element;

            return total;
        }

        /// Fills the vector with zeros.
        pub fn zero(self: *@This()) void {
            self.fill(0);
        }
    };
}
