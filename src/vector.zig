const std = @import("std");

pub fn Vector(comptime T: type) type {
    return struct {
        data: []T,
        shape: [1]usize,

        pub fn init(size: usize, gpa: std.mem.Allocator) !@This() {
            return .{ .data = try gpa.alloc(T, size), .shape = .{size} };
        }

        pub fn initZero(size: usize, gpa: std.mem.Allocator) !@This() {
            const v = try @This().init(size, gpa);

            @memset(v.data, std.mem.zeroes(T));

            return v;
        }

        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            gpa.free(self.data);
        }

        pub fn at(self: @This(), i: usize) T {
            return self.data[i];
        }

        pub fn divs(self: *@This(), scalar: T) void {
            if (comptime @typeInfo(T) == .@"struct") for (self.data) |*e| {
                e.* = e.div(scalar);
            };

            if (comptime @typeInfo(T) != .@"struct") for (self.data) |*e| {
                e.* /= scalar;
            };
        }

        pub fn muls(self: *@This(), scalar: T) void {
            if (comptime @typeInfo(T) == .@"struct") for (self.data) |*e| {
                e.* = e.mul(scalar);
            };

            if (comptime @typeInfo(T) != .@"struct") for (self.data) |*e| {
                e.* *= scalar;
            };
        }

        pub fn ptr(self: *@This(), i: usize) *T {
            return &self.data[i];
        }

        pub fn length(self: @This()) usize {
            return self.shape[0];
        }
    };
}
