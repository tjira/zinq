const std = @import("std");

// zig fmt: off
const Matrix     = @import("tensor.zig")    .Matrix;
const Vector     = @import("tensor.zig")    .Vector;
const ScalarDual = @import("dual.zig"  ).ScalarDual;
const Value      = @import("value.zig" )     .Value;
// zig fmt: on

const eighSlice = @import("openblas.zig").eighSlice;

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

        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, t: U) void {
            switch (self) {
                inline else => |field| field.eval(U, V, r, t),
            }
        }

        pub fn evalBatch(self: @This(), comptime U: type, V: *Matrix(U), r: Matrix(U), t: T) void {
            const t_val = Value(U).fromFloat(t);

            switch (self) {
                inline else => |field| {
                    for (0..r.nrow()) |i| field.eval(U, V.rowSlice(i), r.rowSlice(i), t_val.val);
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

        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, _: U) void {
            var sum = Value(U).fromFloat(0);

            for (r, 0..) |e, i| {
                sum = sum.add(Value(U).init(e).mul(Value(U).init(e)).muls(0.5 * self.k[i]));
            }

            V[0] = sum.val;
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

        pub fn eval(self: @This(), comptime U: type, V: []U, _: []const U, t: U) void {
            const a = Value(U).fromFloat(self.a);
            const g = Value(U).fromFloat(self.g);

            const V0 = Value(U).init(t).sub(a).mul(a);

            // zig fmt: off
            V[0] = V0.val; V[1] =        g.val;
            V[2] =   V[1]; V[3] = V0.neg().val;
            // zig fmt: on
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

        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, _: U) void {
            const r0 = Value(U).init(r[0]);

            const A = Value(U).fromFloat(self.A);
            const B = Value(U).fromFloat(self.B);
            const C = Value(U).fromFloat(self.C);
            const D = Value(U).fromFloat(self.D);

            const v0 = r0.abs().mul(B).neg().exp().neg().adds(1).mul(A).mul(r0.sign());
            const v1 = r0.mul(r0).mul(D).neg().exp().mul(C);

            // zig fmt: off
            V[0] = v0.val; V[1] =       v1.val;
            V[2] =   V[1]; V[3] = v0.neg().val;
            // zig fmt: on
        }

        pub fn ndim(_: @This()) usize {
            return 1;
        }

        pub fn nstate(_: @This()) usize {
            return 2;
        }
    };
}
