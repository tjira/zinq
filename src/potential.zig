const std = @import("std");

const Matrix = @import("matrix.zig").Matrix;
const Vector = @import("vector.zig").Vector;

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
        harmonic: Harmonic(T),
        tully_1: Tully1(T),

        pub fn init(options: Options) @This() {
            return switch (options) {
                .harmonic => |f| .{ .harmonic = Harmonic(T).init(f.k) },
                .tully_1 => |f| .{ .tully_1 = Tully1(T).init(f.A, f.B, f.C, f.D) },
            };
        }

        pub fn evalMany(self: @This(), V: *Matrix(T), r: Matrix(T)) void {
            switch (self) {
                inline else => |field| evalAny(field, V, r),
            }
        }

        pub fn nstate(self: @This()) usize {
            return switch (self) {
                inline else => |field| field.nstate(),
            };
        }

        fn evalAny(field: anytype, V: *Matrix(T), r: Matrix(T)) void {
            const ftype = @TypeOf(field);

            if (comptime @hasDecl(ftype, "evalCoupled")) {
                return evalCoupled(field, V, r);
            }

            if (comptime @hasDecl(ftype, "evalSeparable")) {
                return evalSeparable(field, V, r);
            }

            @compileError("POTENTIAL DOES NOT IMPLEMENT CORRECT EVAL METHOD");
        }

        fn evalCoupled(field: anytype, V: *Matrix(T), r: Matrix(T)) void {
            for (0..r.nrow()) |i| {
                field.evalCoupled(V, r, i);
            }
        }

        fn evalSeparable(field: anytype, V: *Matrix(T), r: Matrix(T)) void {
            V.zero();

            for (0..r.ncol()) |j| {
                for (0..r.nrow()) |i| {
                    field.evalSeparable(V, r, i, j);
                }
            }
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

        pub fn evalSeparable(self: @This(), V: *Matrix(T), r: Matrix(T), i: usize, j: usize) void {
            const rj = r.at(i, j);

            V.ptr(0, i).* += 0.5 * self.k[j] * rj * rj;
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

        pub fn evalCoupled(self: @This(), V: *Matrix(T), r: Matrix(T), i: usize) void {
            const r0 = r.at(i, 0);

            const V00 = std.math.sign(r0) * self.A * (1 - std.math.exp(-self.B * @abs(r0)));
            const V01 = self.C * std.math.exp(-self.D * r0 * r0);
            const V11 = -V00;

            V.ptr(0, i).* = V00;
            V.ptr(1, i).* = V01;
            V.ptr(2, i).* = V01;
            V.ptr(3, i).* = V11;
        }

        pub fn nstate(_: @This()) usize {
            return 2;
        }
    };
}
