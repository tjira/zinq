const std = @import("std");

const Allocator = std.mem.Allocator;

const Value = @import("value.zig").Value;

const isFloat = @import("value.zig").isFloat;

pub fn Integrator(comptime T: type) type {
    const U = if (isFloat(T)) T else @typeInfo(T).@"struct".fields[0].type;

    return struct {
        method: Method,

        pub const Method = union(enum) {
            euler: Euler(T),
        };

        pub fn init(method: Method) @This() {
            return .{ .method = method };
        }

        pub fn step(self: *@This(), y: []T, dt: U, ctx: anytype, comptime dFn: anytype) void {
            switch (self.method) {
                .euler => |*euler| euler.step(y, dt, ctx, dFn),
            }
        }
    };
}

pub fn Euler(comptime T: type) type {
    const U = if (isFloat(T)) T else @typeInfo(T).@"struct".fields[0].type;

    return struct {
        dy: []T,

        pub fn init(nstate: usize, gpa: Allocator) !@This() {
            return .{
                .dy = try gpa.alloc(T, nstate),
            };
        }

        pub fn step(self: *@This(), y: []T, dt: U, ctx: anytype, comptime dFn: anytype) void {
            dFn(ctx, y, self.dy);

            for (0..y.len) |i| {
                y[i] = Value(T).init(y[i]).add(Value(T).init(self.dy[i]).muls(dt)).val;
            }
        }
    };
}
