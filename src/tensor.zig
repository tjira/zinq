const std = @import("std");

const Value = @import("value.zig").Value;

pub fn Matrix(comptime T: type) type {
    return struct {
        data: []T,
        shape: [2]usize,

        pub fn init(rows: usize, cols: usize, gpa: std.mem.Allocator) !@This() {
            return .{ .data = try gpa.alloc(T, rows * cols), .shape = .{ rows, cols } };
        }

        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            gpa.free(self.data);
        }

        pub fn at(self: @This(), i: usize, j: usize) T {
            std.debug.assert(i < self.shape[0]);
            std.debug.assert(j < self.shape[1]);

            return self.data[i * self.shape[1] + j];
        }

        pub fn clone(self: @This(), gpa: std.mem.Allocator) !@This() {
            var A = try @This().init(self.shape[0], self.shape[1], gpa);

            for (0..self.data.len) |i| {
                A.data[i] = self.data[i];
            }

            return A;
        }

        pub fn divs(self: *@This(), scalar: T) void {
            for (0..self.data.len) |i| {
                self.data[i] = Value(T).init(self.data[i]).div(Value(T).init(scalar)).val;
            }
        }

        pub fn fill(self: *@This(), scalar: T) void {
            for (0..self.data.len) |i| {
                self.data[i] = scalar;
            }
        }

        pub fn initZero(rows: usize, cols: usize, gpa: std.mem.Allocator) !@This() {
            var A = try @This().init(rows, cols, gpa);

            A.zero();

            return A;
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

        pub fn takeRows(self: @This(), n: usize) @This() {
            std.debug.assert(n <= self.shape[0]);

            return .{ .data = self.data[0 .. n * self.shape[1]], .shape = .{ n, self.shape[1] } };
        }

        pub fn zero(self: *@This()) void {
            self.fill(std.mem.zeroes(T));
        }
    };
}

pub fn Vector(comptime T: type) type {
    return struct {
        data: []T,
        shape: [1]usize,

        pub fn init(size: usize, gpa: std.mem.Allocator) !@This() {
            return .{ .data = try gpa.alloc(T, size), .shape = .{size} };
        }

        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            gpa.free(self.data);
        }

        pub fn asMatrix(self: @This()) Matrix(T) {
            return .{ .data = self.data, .shape = .{ self.shape[0], 1 } };
        }

        pub fn at(self: @This(), i: usize) T {
            std.debug.assert(i < self.shape[0]);

            return self.data[i];
        }

        pub fn divs(self: *@This(), scalar: T) void {
            for (0..self.data.len) |i| {
                self.data[i] = Value(T).init(self.data[i]).div(Value(T).init(scalar)).val;
            }
        }

        pub fn initZero(size: usize, gpa: std.mem.Allocator) !@This() {
            var v = try @This().init(size, gpa);

            v.zero();

            return v;
        }

        pub fn length(self: @This()) usize {
            return self.shape[0];
        }

        pub fn muls(self: *@This(), scalar: T) void {
            for (0..self.data.len) |i| {
                self.data[i] = Value(T).init(self.data[i]).mul(Value(T).init(scalar)).val;
            }
        }

        pub fn ptr(self: *@This(), i: usize) *T {
            std.debug.assert(i < self.shape[0]);

            return &self.data[i];
        }

        pub fn takeRows(self: @This(), n: usize) @This() {
            std.debug.assert(n <= self.shape[0]);

            return .{ .data = self.data[0..n], .shape = .{n} };
        }

        pub fn zero(self: *@This()) void {
            for (0..self.data.len) |i| {
                self.data[i] = std.mem.zeroes(T);
            }
        }
    };
}

pub fn Tensor(comptime T: type, comptime N: usize) type {
    return struct {
        data: []T,
        shape: [N]usize,

        pub fn init(shape: [N]usize, gpa: std.mem.Allocator) !@This() {
            var size: usize = 1;

            inline for (0..N) |i| {
                size *= shape[i];
            }

            return .{ .data = try gpa.alloc(T, size), .shape = shape };
        }

        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            gpa.free(self.data);
        }

        pub fn asMatrix(self: @This()) Matrix(T) {
            var rows: usize = 1;
            var cols: usize = 1;

            for (0..(N + 1) / 2) |i| {
                rows *= self.shape[i];
            }

            for ((N + 1) / 2..N) |i| {
                cols *= self.shape[i];
            }

            return .{ .data = self.data, .shape = .{ rows, cols } };
        }

        pub fn at(self: @This(), indx: [N]usize) T {
            var idx: usize = 0;
            var str: usize = 1;

            inline for (0..N) |k| {
                const j = N - 1 - k;

                std.debug.assert(indx[j] < self.shape[j]);

                idx += indx[j] * str;
                str *= self.shape[j];
            }

            return self.data[idx];
        }

        pub fn initZero(shape: [N]usize, gpa: std.mem.Allocator) !@This() {
            var U = try @This().init(shape, gpa);

            U.zero();

            return U;
        }

        pub fn ptr(self: *@This(), indx: [N]usize) *T {
            var idx: usize = 0;
            var str: usize = 1;

            inline for (0..N) |k| {
                const j = N - 1 - k;

                std.debug.assert(indx[j] < self.shape[j]);

                idx += indx[j] * str;
                str *= self.shape[j];
            }

            return &self.data[idx];
        }

        pub fn zero(self: *@This()) void {
            for (0..self.data.len) |i| {
                self.data[i] = std.mem.zeroes(T);
            }
        }
    };
}
