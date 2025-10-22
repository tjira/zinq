//! Array for storing object arrays.

const std = @import("std");

const device_write = @import("device_write.zig");
const ring_buffer = @import("ring_buffer.zig");
const real_matrix = @import("real_matrix.zig");

const RingBuffer = ring_buffer.RingBuffer;
const RealMatrix = real_matrix.RealMatrix;

const throw = device_write.throw;

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
                RingBuffer(T) => |object| element.* = try object.init(params.max_len, allocator),
                RealMatrix(T) => |object| element.* = try object.init(params.rows, params.cols, allocator),
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
                RingBuffer(T) => |object| element.* = try object.initZero(params.max_len, allocator),
                RealMatrix(T) => |object| element.* = try object.initZero(params.rows, params.cols, allocator),
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

/// Array for storing ring buffers in an array.
pub fn RingBufferArray(comptime T: type) type {
    return ObjectArray(RingBuffer, T);
}

/// Array for storing real matrices.
pub fn RealMatrixArray(comptime T: type) type {
    return ObjectArray(RealMatrix, T);
}
