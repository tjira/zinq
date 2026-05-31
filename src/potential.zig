const std = @import("std");

const Matrix = @import("tensor.zig").Matrix;
const Vector = @import("tensor.zig").Vector;

// OPTIONS =============================================================================================================

pub const Options = union(enum) {
    harmonic: struct {
        k: []const f64 = &.{1},
    },
    tully_1: struct {
        A: f64 = 0.01,
        B: f64 = 1.6,
        C: f64 = 0.005,
        D: f64 = 1.0,
    },
};

// GENERIC POTENTIAL ===================================================================================================

pub fn Potential(comptime T: type) type {
    return union(enum) {
        // zig fmt: off
        harmonic: Harmonic(T),
        tully_1:  Tully1  (T),
        // zig fmt: on

        pub fn init(options: Options) @This() {
            return switch (options) {
                // zig fmt: off
                .harmonic => |f| .{ .harmonic = Harmonic(T).init(f.k               ) },
                .tully_1  => |f| .{ .tully_1  = Tully1  (T).init(f.A, f.B, f.C, f.D) },
                // zig fmt: on
            };
        }

        pub fn evalBatch(self: @This(), V: *Matrix(T), r: Matrix(T)) void {
            switch (self) {
                inline else => |field| {
                    for (0..r.ncol()) |i| {
                        const val = field.eval(r.colSlice(i));

                        for (0..field.nstate()) |j| for (0..field.nstate()) |k| {
                            V.ptr(j * field.nstate() + k, i).* = val[j][k];
                        };
                    }
                },
            }
        }

        pub fn ndim(self: @This()) usize {
            return switch (self) {
                inline else => |field| field.ndim(),
            };
        }

        pub fn nstate(self: @This()) usize {
            return switch (self) {
                inline else => |field| field.nstate(),
            };
        }
    };
}

// SPECIFIC POTENTIALS =================================================================================================

pub fn Harmonic(comptime T: type) type {
    return struct {
        k: []const T,

        pub fn init(k: []const T) @This() {
            return .{ .k = k };
        }

        pub fn eval(self: @This(), r: []const T) [1][1]T {
            var V00: T = 0;

            for (r, 0..) |ri, i| {
                V00 += 0.5 * self.k[i] * ri * ri;
            }

            return .{
                .{V00},
            };
        }

        pub fn ndim(self: @This()) usize {
            return self.k.len;
        }

        pub fn nstate(_: @This()) usize {
            return 1;
        }
    };
}

pub fn Tully1(comptime T: type) type {
    return struct {
        A: T,
        B: T,
        C: T,
        D: T,

        pub fn init(A: T, B: T, C: T, D: T) @This() {
            return .{ .A = A, .B = B, .C = C, .D = D };
        }

        pub fn eval(self: @This(), r: []const T) [2][2]T {
            const V00 = std.math.sign(r[0]) * self.A * (1 - std.math.exp(-self.B * @abs(r[0])));
            const V01 = self.C * std.math.exp(-self.D * r[0] * r[0]);
            const V11 = -V00;

            return .{
                .{ V00, V01 },
                .{ V01, V11 },
            };
        }

        pub fn ndim(_: @This()) usize {
            return 1;
        }

        pub fn nstate(_: @This()) usize {
            return 2;
        }
    };
}
