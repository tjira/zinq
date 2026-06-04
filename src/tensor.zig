const std = @import("std");

const div = @import("generic_op.zig").div;
const mul = @import("generic_op.zig").mul;

// MATRIX ==============================================================================================================

pub fn Matrix(comptime T: type) type {
    return struct {
        // zig fmt: off
        data: []T, shape: [2]usize,
        // zig fmt: on

        pub fn init(rows: usize, cols: usize, gpa: std.mem.Allocator) !@This() {
            return .{ .data = try gpa.alloc(T, rows * cols), .shape = .{ rows, cols } };
        }

        pub fn initZero(rows: usize, cols: usize, gpa: std.mem.Allocator) !@This() {
            var A = try @This().init(rows, cols, gpa);

            A.zero();

            return A;
        }

        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            gpa.free(self.data);
        }

        pub fn clone(self: @This(), gpa: std.mem.Allocator) !@This() {
            var A = try @This().init(self.shape[0], self.shape[1], gpa);

            for (self.data, 0..) |e, i| {
                A.data[i] = e;
            }

            return A;
        }

        pub fn at(self: @This(), i: usize, j: usize) T {
            return self.data[i * self.shape[1] + j];
        }

        pub fn divs(self: *@This(), scalar: T) void {
            for (self.data) |*e| {
                e.* = div(e.*, scalar);
            }
        }

        pub fn fill(self: *@This(), scalar: T) void {
            for (self.data) |*e| {
                e.* = scalar;
            }
        }

        pub fn ncol(self: @This()) usize {
            return self.shape[1];
        }

        pub fn nrow(self: @This()) usize {
            return self.shape[0];
        }

        pub fn ptr(self: *@This(), i: usize, j: usize) *T {
            return &self.data[i * self.shape[1] + j];
        }

        pub fn rowSlice(self: @This(), i: usize) []T {
            return self.data[i * self.shape[1] .. (i + 1) * self.shape[1]];
        }

        pub fn slice(self: @This(), start: usize, end: usize) []T {
            return self.data[start..end];
        }

        pub fn zero(self: *@This()) void {
            self.fill(std.mem.zeroes(T));
        }
    };
}

// VECTOR ==============================================================================================================

pub fn Vector(comptime T: type) type {
    return struct {
        // zig fmt: off
        data: []T, shape: [1]usize,
        // zig fmt: on

        pub fn init(size: usize, gpa: std.mem.Allocator) !@This() {
            return .{ .data = try gpa.alloc(T, size), .shape = .{size} };
        }

        pub fn initZero(size: usize, gpa: std.mem.Allocator) !@This() {
            var v = try @This().init(size, gpa);

            v.zero();

            return v;
        }

        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            gpa.free(self.data);
        }

        pub fn at(self: @This(), i: usize) T {
            return self.data[i];
        }

        pub fn divs(self: *@This(), scalar: T) void {
            for (self.data) |*e| {
                e.* = div(e.*, scalar);
            }
        }

        pub fn length(self: @This()) usize {
            return self.shape[0];
        }

        pub fn muls(self: *@This(), scalar: T) void {
            for (self.data) |*e| {
                e.* = mul(e.*, scalar);
            }
        }

        pub fn ptr(self: *@This(), i: usize) *T {
            return &self.data[i];
        }

        pub fn zero(self: *@This()) void {
            for (self.data) |*e| {
                e.* = std.mem.zeroes(T);
            }
        }
    };
}
