const std = @import("std");

const Value = @import("value.zig").Value;

// MATRIX ==============================================================================================================

pub fn Matrix(comptime T: type) type {
    return struct {
        data: []T,
        shape: [2]usize,

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
            std.debug.assert(i < self.shape[0]);
            std.debug.assert(j < self.shape[1]);

            return self.data[i * self.shape[1] + j];
        }

        pub fn divs(self: *@This(), scalar: T) void {
            for (self.data) |*e| {
                e.* = Value(T).init(e.*).div(Value(T).init(scalar)).val;
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
            std.debug.assert(i < self.shape[0]);
            std.debug.assert(j < self.shape[1]);

            return &self.data[i * self.shape[1] + j];
        }

        pub fn rowSlice(self: @This(), i: usize) []T {
            std.debug.assert(i < self.shape[0]);

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
        data: []T,
        shape: [1]usize,

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
            std.debug.assert(i < self.shape[0]);

            return self.data[i];
        }

        pub fn divs(self: *@This(), scalar: T) void {
            for (self.data) |*e| {
                e.* = Value(T).init(e.*).div(Value(T).init(scalar)).val;
            }
        }

        pub fn length(self: @This()) usize {
            return self.shape[0];
        }

        pub fn muls(self: *@This(), scalar: T) void {
            for (self.data) |*e| {
                e.* = Value(T).init(e.*).mul(Value(T).init(scalar)).val;
            }
        }

        pub fn ptr(self: *@This(), i: usize) *T {
            std.debug.assert(i < self.shape[0]);

            return &self.data[i];
        }

        pub fn zero(self: *@This()) void {
            for (self.data) |*e| {
                e.* = std.mem.zeroes(T);
            }
        }
    };
}
