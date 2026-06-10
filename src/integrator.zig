const std = @import("std");

const Allocator = std.mem.Allocator;

const Value = @import("value.zig").Value;

const primType = @import("value.zig").primType;

pub fn ButcherTableau(comptime T: type, comptime STAGES: usize) type {
    return struct {
        a: [STAGES][STAGES]T,

        b: [STAGES]T,
        c: [STAGES]T,
    };
}

pub fn Integrator(comptime T: type) type {
    const U = primType(T);

    return struct {
        method: Method,

        pub const Method = union(enum) {
            rk1: Rk1(T),
            rk4: Rk4(T),
        };

        pub fn init(tag: std.meta.Tag(Method), nstate: usize, gpa: Allocator) !@This() {
            inline for (std.meta.fields(Method)) |field| if (tag == @field(std.meta.Tag(Method), field.name)) {
                return .{ .method = @unionInit(Method, field.name, try field.type.init(nstate, gpa)) };
            };

            unreachable;
        }

        pub fn deinit(self: @This(), gpa: Allocator) void {
            switch (self.method) {
                inline else => |*method| method.deinit(gpa),
            }
        }

        pub fn step(self: *@This(), y: []T, dt: U, ctx: anytype, comptime dFn: anytype) void {
            switch (self.method) {
                inline else => |*method| method.step(y, dt, ctx, dFn),
            }
        }
    };
}

fn RungeKutta(comptime T: type, comptime tab: anytype) type {
    const U = primType(T);

    return struct {
        k: [tab.b.len][]T,

        tmp: []T,

        pub fn init(nstate: usize, gpa: Allocator) !@This() {
            var self: @This() = undefined;

            const k_mem = try gpa.alloc(T, tab.b.len * nstate);
            errdefer gpa.free(k_mem);

            for (0..self.k.len) |i| {
                self.k[i] = k_mem[i * nstate .. (i + 1) * nstate];
            }

            self.tmp = try gpa.alloc(T, nstate);

            return self;
        }

        pub fn deinit(self: @This(), gpa: Allocator) void {
            gpa.free(self.k[0].ptr[0 .. self.k.len * self.tmp.len]);

            gpa.free(self.tmp);
        }

        pub fn step(self: *@This(), y: []T, dt: U, ctx: anytype, comptime dFn: anytype) void {
            std.debug.assert(self.tmp.len == y.len);

            for (self.k) |ki| {
                std.debug.assert(ki.len == y.len);
            }

            inline for (0..self.k.len) |i| {
                for (0..y.len) |j| {
                    var sum = Value(T).init(y[j]);

                    inline for (0..i) |k| if (tab.a[i][k] != 0) {
                        sum = sum.add(Value(T).init(self.k[k][j]).muls(tab.a[i][k] * dt));
                    };

                    self.tmp[j] = sum.val;
                }

                dFn(ctx, self.tmp, self.k[i]);
            }

            for (0..y.len) |i| {
                var sum = Value(T).init(y[i]);

                inline for (0..self.k.len) |j| if (tab.b[j] != 0.0) {
                    sum = sum.add(Value(T).init(self.k[j][i]).muls(tab.b[j] * dt));
                };

                y[i] = sum.val;
            }
        }
    };
}

fn Rk1(comptime T: type) type {
    return RungeKutta(T, rk1Tableau(primType(T)));
}

fn rk1Tableau(comptime U: type) ButcherTableau(U, 1) {
    return .{
        .a = .{
            .{0.0},
        },

        .b = .{1.0},
        .c = .{0.0},
    };
}

fn Rk4(comptime T: type) type {
    return RungeKutta(T, rk4Tableau(primType(T)));
}

fn rk4Tableau(comptime U: type) ButcherTableau(U, 4) {
    return .{
        .a = .{
            .{ 0.0, 0.0, 0.0, 0.0 },
            .{ 0.5, 0.0, 0.0, 0.0 },
            .{ 0.0, 0.5, 0.0, 0.0 },
            .{ 0.0, 0.0, 1.0, 0.0 },
        },

        .b = .{ 1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0 },
        .c = .{ 0.0 / 1.0, 1.0 / 2.0, 1.0 / 2.0, 1.0 / 1.0 },
    };
}
