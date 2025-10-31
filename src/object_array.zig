//! Array for storing object arrays.

const std = @import("std");

const complex_vector = @import("complex_vector.zig");
const complex_matrix = @import("complex_matrix.zig");
const device_write = @import("device_write.zig");
const error_handling = @import("error_handling.zig");
const ring_buffer = @import("ring_buffer.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const ComplexVector = complex_vector.ComplexVector;
const ComplexMatrix = complex_matrix.ComplexMatrix;
const RingBuffer = ring_buffer.RingBuffer;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const throw = error_handling.throw;

/// Array for storing object in an array.
pub fn ObjectArray(comptime O: fn (comptime type) type, comptime T: type) type {
    return struct {
        data: []O(T),
        len: usize,
        allocator: std.mem.Allocator,

        /// Initialize the object array.
        pub fn init(len: usize, params: anytype, allocator: std.mem.Allocator) !@This() {
            const data = try allocator.alloc(O(T), len);

            for (data) |*element| switch (O(T)) {
                ComplexMatrix(T) => |object| element.* = try object.init(params.rows, params.cols, allocator),
                ComplexVector(T) => |object| element.* = try object.init(params.len, allocator),
                RingBuffer(T) => |object| element.* = try object.init(params.max_len, allocator),
                RealMatrix(T) => |object| element.* = try object.init(params.rows, params.cols, allocator),
                RealVector(T) => |object| element.* = try object.init(params.rows, allocator),
                else => return throw(@This(), "UNSUPPORTED OBJECT TYPE FOR OBJECTARRAY", .{}),
            };

            return @This(){
                .data = data,
                .len = len,
                .allocator = allocator,
            };
        }

        /// Initialize the object array to zero.
        pub fn initZero(len: usize, params: anytype, allocator: std.mem.Allocator) !@This() {
            const data = try allocator.alloc(O(T), len);

            for (data) |*element| switch (O(T)) {
                ComplexMatrix(T) => |object| element.* = try object.initZero(params.rows, params.cols, allocator),
                ComplexVector(T) => |object| element.* = try object.initZero(params.len, allocator),
                RingBuffer(T) => |object| element.* = try object.initZero(params.max_len, allocator),
                RealMatrix(T) => |object| element.* = try object.initZero(params.rows, params.cols, allocator),
                RealVector(T) => |object| element.* = try object.initZero(params.rows, allocator),
                else => return throw(@This(), "UNSUPPORTED OBJECT TYPE FOR OBJECTARRAY", .{}),
            };

            return @This(){
                .data = data,
                .len = len,
                .allocator = allocator,
            };
        }

        /// Deinitialize the array array.
        pub fn deinit(self: @This()) void {
            for (self.data) |object| {
                object.deinit();
            }

            self.allocator.free(self.data);
        }

        /// Element access operator to get the i-th object.
        pub fn at(self: @This(), i: usize) O(T) {
            return self.data[i];
        }

        /// Element pointer access operator to get the i-th object.
        pub fn ptr(self: @This(), i: usize) *O(T) {
            return &self.data[i];
        }
    };
}

/// Array for storing complex matrices.
pub fn ComplexMatrixArray(comptime T: type) type {
    return ObjectArray(ComplexMatrix, T);
}

/// Array for storing complex vectors.
pub fn ComplexVectorArray(comptime T: type) type {
    return ObjectArray(ComplexVector, T);
}

/// Array for storing ring buffers in an array.
pub fn RingBufferArray(comptime T: type) type {
    return ObjectArray(RingBuffer, T);
}

/// Array for storing real matrices.
pub fn RealMatrixArray(comptime T: type) type {
    return ObjectArray(RealMatrix, T);
}

/// Array for storing real vectors.
pub fn RealVectorArray(comptime T: type) type {
    return ObjectArray(RealVector, T);
}
