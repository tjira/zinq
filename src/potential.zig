const std = @import("std");

// zig fmt: off
const Matrix     = @import("tensor.zig")    .Matrix;
const Vector     = @import("tensor.zig")    .Vector;
const ScalarDual = @import("dual.zig"  ).ScalarDual;
// zig fmt: on

// OPTIONS =============================================================================================================

pub const Options = union(enum) {
    harmonic: struct {
        k: []const f64 = &.{1},
    },
    time_linear: struct {
        // zig fmt: off
        a: f64 = 10,
        g: f64 =  2,
        // zig fmt: on
    },
    tully_1: struct {
        A: f64 = 0.010,
        B: f64 = 1.600,
        C: f64 = 0.005,
        D: f64 = 1.000,
    },
};

// GENERIC POTENTIAL ===================================================================================================

pub fn Potential(comptime T: type) type {
    return union(enum) {
        // zig fmt: off
        harmonic:    Harmonic  (T),
        time_linear: TimeLinear(T),
        tully_1:     Tully1    (T),
        // zig fmt: on

        pub fn init(options: Options) @This() {
            return switch (options) {
                // zig fmt: off
                .harmonic    => |f| .{ .harmonic    = Harmonic  (T).init(f.k               ) },
                .time_linear => |f| .{ .time_linear = TimeLinear(T).init(f.a, f.g          ) },
                .tully_1     => |f| .{ .tully_1     = Tully1    (T).init(f.A, f.B, f.C, f.D) },
                // zig fmt: on
            };
        }

        pub fn eval(self: @This(), V: []T, r: []const T, t: T) void {
            switch (self) {
                inline else => |field| field.eval(V, r, t),
            }
        }

        pub fn evalBatch(self: @This(), V: *Matrix(T), r: Matrix(T), t: T) void {
            switch (self) {
                inline else => |field| {
                    for (0..r.nrow()) |i| field.eval(V.rowSlice(i), r.rowSlice(i), t);
                },
            }
        }

        pub fn evalDualBatch(self: @This(), V: *Matrix(ScalarDual(T)), r: Matrix(ScalarDual(T)), t: T) void {
            const t_dual = ScalarDual(T).init(t, 0.0);

            switch (self) {
                inline else => |field| {
                    for (0..r.nrow()) |i| field.evalDual(V.rowSlice(i), r.rowSlice(i), t_dual);
                },
            }
        }

        pub fn is_td(self: @This()) bool {
            return switch (self) {
                .time_linear => true,
                inline else => false,
            };
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

        pub fn eval(self: @This(), V: []T, r: []const T, _: T) void {
            V[0] = 0;

            for (r, 0..) |e, i| {
                V[0] += 0.5 * self.k[i] * e * e;
            }
        }

        pub fn evalDual(self: @This(), V: []ScalarDual(T), r: []const ScalarDual(T), _: ScalarDual(T)) void {
            V[0] = ScalarDual(T).init(0, 0);

            for (r, 0..) |e, i| {
                V[0] = V[0].add(e.mul(e).muls(0.5 * self.k[i]));
            }
        }

        pub fn ndim(self: @This()) usize {
            return self.k.len;
        }

        pub fn nstate(_: @This()) usize {
            return 1;
        }
    };
}

pub fn TimeLinear(comptime T: type) type {
    return struct {
        a: T,
        g: T,

        pub fn init(a: T, g: T) @This() {
            return .{ .a = a, .g = g };
        }

        pub fn eval(self: @This(), V: []T, _: []const T, t: T) void {
            V[0] = self.a * (t - self.a);
            V[1] = self.g;
            V[2] = V[1];
            V[3] = -V[0];
        }

        pub fn evalDual(self: @This(), V: []ScalarDual(T), _: []const ScalarDual(T), t: ScalarDual(T)) void {
            V[0] = t.subs(self.a).muls(self.a);
            V[1] = ScalarDual(T).init(self.g, 0);
            V[2] = V[1];
            V[3] = V[0].muls(-1);
        }

        pub fn ndim(_: @This()) usize {
            return 1;
        }

        pub fn nstate(_: @This()) usize {
            return 2;
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

        pub fn eval(self: @This(), V: []T, r: []const T, _: T) void {
            V[0] = std.math.sign(r[0]) * self.A * (1 - std.math.exp(-self.B * @abs(r[0])));
            V[1] = self.C * std.math.exp(-self.D * r[0] * r[0]);
            V[2] = V[1];
            V[3] = -V[0];
        }

        pub fn evalDual(self: @This(), V: []ScalarDual(T), r: []const ScalarDual(T), _: ScalarDual(T)) void {
            V[0] = r[0].abs().muls(-self.B).exp().muls(-1).adds(1).muls(std.math.sign(r[0].val) * self.A);
            V[1] = r[0].mul(r[0]).muls(-self.D).exp().muls(self.C);
            V[2] = V[1];
            V[3] = V[0].muls(-1);
        }

        pub fn ndim(_: @This()) usize {
            return 1;
        }

        pub fn nstate(_: @This()) usize {
            return 2;
        }
    };
}
