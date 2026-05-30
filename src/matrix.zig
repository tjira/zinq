const std = @import("std");

const Vector = @import("vector.zig").Vector;

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

        pub fn at(self: @This(), i: usize, j: usize) T {
            return self.data[j * self.shape[0] + i];
        }

        pub fn colSlice(self: @This(), j: usize) []T {
            return self.data[j * self.shape[0] .. (j + 1) * self.shape[0]];
        }

        pub fn divs(self: *@This(), scalar: T) void {
            if (comptime @typeInfo(T) == .@"struct") for (self.data) |*e| {
                e.* = e.div(scalar);
            };

            if (comptime @typeInfo(T) != .@"struct") for (self.data) |*e| {
                e.* /= scalar;
            };
        }

        pub fn ncol(self: @This()) usize {
            return self.shape[1];
        }

        pub fn nrow(self: @This()) usize {
            return self.shape[0];
        }

        pub fn ptr(self: *@This(), i: usize, j: usize) *T {
            return &self.data[j * self.shape[0] + i];
        }

        pub fn slice(self: @This(), start: usize, end: usize) []T {
            return self.data[start..end];
        }

        pub fn zero(self: *@This()) void {
            @memset(self.data, std.mem.zeroes(T));
        }
    };
}
