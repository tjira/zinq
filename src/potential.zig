const std = @import("std");

const Allocator = std.mem.Allocator;

const Expression = @import("exprtk.zig").Expression;
const Matrix = @import("tensor.zig").Matrix;
const Value = @import("value.zig").Value;
const Vector = @import("tensor.zig").Vector;

const isDual = @import("value.zig").isDual;
const readMatrix = @import("read_write.zig").readMatrix;

pub const Options = union(enum) {
    file: struct {
        ndim: u32,
        path: []const u8,
    },
    harmonic: struct {
        k: []const f64 = &.{1},
    },
    time_linear: struct {
        a: f64 = 10,
        g: f64 = 2,
    },
    tully_1: struct {
        A: f64 = 0.010,
        B: f64 = 1.600,
        C: f64 = 0.005,
        D: f64 = 1.000,
    },
    custom: struct {
        matrix: []const []const []const u8,
        ndim: u32,
        time_dependent: bool = false,
    },
};

pub fn Potential(comptime T: type) type {
    return union(enum) {
        file: File(T),
        harmonic: Harmonic(T),
        time_linear: TimeLinear(T),
        tully_1: Tully1(T),
        custom: Custom(T),

        pub fn init(io: std.Io, options: Options, allocator: Allocator) !@This() {
            return switch (options) {
                .file => |f| .{ .file = try File(T).init(f.ndim, f.path, io, allocator) },
                .harmonic => |f| .{ .harmonic = Harmonic(T).init(f.k) },
                .time_linear => |f| .{ .time_linear = TimeLinear(T).init(f.a, f.g) },
                .tully_1 => |f| .{ .tully_1 = Tully1(T).init(f.A, f.B, f.C, f.D) },
                .custom => |f| .{ .custom = try Custom(T).init(f.ndim, f.matrix, f.time_dependent, allocator) },
            };
        }

        pub fn deinit(self: *@This(), allocator: Allocator) void {
            switch (self.*) {
                inline .file, .custom => |*pot| pot.deinit(allocator),

                inline else => {},
            }
        }

        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, t: U) void {
            std.debug.assert(V.len == self.nstate() * self.nstate());

            std.debug.assert(r.len == self.ndim());

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

        pub fn isTd(self: @This()) bool {
            return switch (self) {
                inline else => |field| field.isTd(),
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

fn Harmonic(comptime T: type) type {
    return struct {
        k: []const T,

        pub fn init(k: []const T) @This() {
            return .{ .k = k };
        }

        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, _: U) void {
            var sum = Value(U).fromFloat(0);

            for (0..r.len) |i| {
                sum = sum.add(Value(U).init(r[i]).mul(Value(U).init(r[i])).muls(0.5 * self.k[i]));
            }

            V[0] = sum.val;
        }

        pub fn isTd(_: @This()) bool {
            return false;
        }

        pub fn ndim(self: @This()) usize {
            return self.k.len;
        }

        pub fn nstate(_: @This()) usize {
            return 1;
        }
    };
}

fn TimeLinear(comptime T: type) type {
    return struct {
        a: T,
        g: T,

        pub fn init(a: T, g: T) @This() {
            return .{ .a = a, .g = g };
        }

        pub fn eval(self: @This(), comptime U: type, V: []U, _: []const U, t: U) void {
            const a = Value(U).fromFloat(self.a);
            const g = Value(U).fromFloat(self.g);

            const V00 = Value(U).init(t).sub(a).mul(a);
            const V01 = g;
            const V11 = V00.neg();

            V[0] = V00.val;
            V[1] = V01.val;
            V[2] = V01.val;
            V[3] = V11.val;
        }

        pub fn isTd(_: @This()) bool {
            return true;
        }

        pub fn ndim(_: @This()) usize {
            return 1;
        }

        pub fn nstate(_: @This()) usize {
            return 2;
        }
    };
}

fn Tully1(comptime T: type) type {
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

            const V00 = r0.abs().mul(B).neg().exp().neg().adds(1).mul(A).mul(r0.sign());
            const V01 = r0.mul(r0).mul(D).neg().exp().mul(C);
            const V11 = V00.neg();

            V[0] = V00.val;
            V[1] = V01.val;
            V[2] = V01.val;
            V[3] = V11.val;
        }

        pub fn isTd(_: @This()) bool {
            return false;
        }

        pub fn ndim(_: @This()) usize {
            return 1;
        }

        pub fn nstate(_: @This()) usize {
            return 2;
        }
    };
}

fn Custom(comptime T: type) type {
    return struct {
        expressions: []const Expression(T),

        is_td: bool,
        r_vals: []T,

        pub fn init(dim: usize, matrix: []const []const []const u8, is_td: bool, allocator: Allocator) !@This() {
            for (0..matrix.len) |i| {
                if (matrix[i].len != matrix.len) return error.NonSquareMatrix;
            }

            const exprs = try allocator.alloc(Expression(T), matrix.len * matrix.len);
            errdefer allocator.free(exprs);

            for (0..matrix.len) |i| for (0..matrix.len) |j| {
                exprs[i * matrix.len + j] = Expression(T).init(matrix[i][j], dim) catch |err| {
                    for (0..(i * matrix.len + j)) |k| {
                        exprs[k].deinit();
                    }

                    return err;
                };
            };

            errdefer for (0..matrix.len * matrix.len) |i| {
                exprs[i].deinit();
            };

            return .{ .expressions = exprs, .r_vals = try allocator.alloc(T, dim), .is_td = is_td };
        }

        pub fn deinit(self: *@This(), allocator: Allocator) void {
            allocator.free(self.r_vals);

            for (self.expressions) |expr| {
                expr.deinit();
            }

            allocator.free(self.expressions);
        }

        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, t: U) void {
            if (comptime U == T) for (0..self.nstate()) |i| for (0..self.nstate()) |j| {
                V[i * self.nstate() + j] = self.expressions[i * self.nstate() + j].evaluate_d0(r, t);
            };

            if (comptime isDual(U)) {
                for (0..r.len) |i| {
                    self.r_vals[i] = r[i].val;
                }

                for (0..self.nstate()) |i| for (0..self.nstate()) |j| {
                    const expr, var der: T = .{ self.expressions[i * self.nstate() + j], 0 };

                    for (0..r.len) |k| {
                        der += expr.evaluate_d1(self.r_vals, t.val, k) * r[k].der;
                    }

                    if (t.der != 0) {
                        der += expr.evaluate_d1(self.r_vals, t.val, r.len) * t.der;
                    }

                    V[i * self.nstate() + j] = U.init(expr.evaluate_d0(self.r_vals, t.val), der);
                };
            }

            if (comptime U != T and !isDual(U)) {
                @compileError("UNSUPPORTED NUMBER TYPE FOR CUSTOM POTENTIAL");
            }
        }

        pub fn isTd(self: @This()) bool {
            return self.is_td;
        }

        pub fn ndim(self: @This()) usize {
            return self.r_vals.len;
        }

        pub fn nstate(self: @This()) usize {
            return std.math.sqrt(self.expressions.len);
        }
    };
}

fn File(comptime T: type) type {
    return struct {
        U: Matrix(T),
        ndims: usize,

        pub fn init(ndims: usize, path: []const u8, io: std.Io, gpa: Allocator) !@This() {
            return .{ .U = try readMatrix(T, io, path, gpa), .ndims = ndims };
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.U.deinit(gpa);
        }

        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, _: U) void {
            for (0..self.nstate()) |i| for (0..self.nstate()) |j| {
                const col = i * self.nstate() + j;

                V[col] = lerp(T, U, self.U, self.ndim() + col, r);
            };
        }

        pub fn isTd(_: @This()) bool {
            return false;
        }

        pub fn ndim(self: @This()) usize {
            return self.ndims;
        }

        pub fn nstate(self: @This()) usize {
            return std.math.sqrt(self.U.ncol() - self.ndim());
        }
    };
}

fn lerp(comptime T: type, comptime U: type, grid: Matrix(T), column: usize, r: []const U) U {
    std.debug.assert(column < grid.ncol());

    const size: usize = @round(std.math.pow(T, @as(T, @floatFromInt(grid.nrow())), 1 / @as(T, @floatFromInt(r.len))));

    var result = Value(U).fromFloat(0);

    for (0..@as(usize, 1) << @as(u5, @intCast(r.len))) |i| {
        var w, var j: usize = .{ Value(U).fromFloat(1), 0 };

        for (0..r.len) |k| {
            const stride = std.math.pow(usize, size, r.len - k - 1);

            var low: usize, var high: usize, var mid: usize = .{ 0, size, size / 2 };

            while (low < high) : (mid = (low + high) / 2) {
                if (grid.at(mid * stride, k) <= if (comptime isDual(U)) r[k].val else r[k]) low = mid + 1 else high = mid;
            }

            const idx = @min(@max(low, 1), size - 1);

            const x0 = grid.at((idx - 1) * stride, k);
            const x1 = grid.at((idx - 0) * stride, k);

            const wgh = Value(U).init(r[k]).sub(Value(U).fromFloat(x0)).muls(1 / (x1 - x0));

            const use_upper = ((i >> @as(u5, @intCast(r.len - k - 1))) & 1) == 1;

            const wmul = if (use_upper) wgh else Value(U).fromFloat(1).sub(wgh);

            w, j = .{ w.mul(wmul), j + (if (use_upper) idx else idx - 1) * stride };
        }

        result = result.add(w.muls(grid.at(j, column)));
    }

    return result.val;
}
