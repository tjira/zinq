//! Defines potential energy surfaces (PES) including harmonic, time-linear coupling, Tully avoided crossing, custom, and file-interpolated types.

const std = @import("std");

const Allocator = std.mem.Allocator;

const Expression = @import("exprtk.zig").Expression;
const Matrix = @import("tensor.zig").Matrix;
const Value = @import("value.zig").Value;
const Vector = @import("tensor.zig").Vector;

const isDual = @import("value.zig").isDual;
const readMatrix = @import("read_write.zig").readMatrix;

const CM2EV = @import("constant.zig").CM2EV;
const EV2AU = @import("constant.zig").EV2AU;

/// Parameter union for potential energy surfaces supporting harmonic, coupling, and Tully 1, 2, and 3 models.
pub const Options = union(enum) {
    file: struct {
        ndim: u32,
        path: []const u8,
    },
    harmonic: struct {
        k: []const f64 = &.{1},
    },
    henon_heiles: struct {
        k: f64 = 1,
        l: f64 = 0.1,
    },
    jahn_teller: struct {
        k: f64 = 1,
        g: f64 = 1,
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
    tully_2: struct {
        A: f64 = 0.10,
        B: f64 = 0.28,
        C: f64 = 0.015,
        D: f64 = 0.06,
        E0: f64 = 0.05,
    },
    tully_3: struct {
        A: f64 = 6.0e-4,
        B: f64 = 0.10,
        C: f64 = 0.90,
    },
    lvc: struct {
        frequencies: []const f64,
        excitation_energies: []const f64,
        kappa: []const []const f64,
        lambda: []const []const []const f64,
    },
    custom: struct {
        matrix: []const []const []const u8,
        ndim: u32,
        time_dependent: bool = false,
    },
};

/// Returns a generic union representing a potential energy surface (PES) with coordinate and time evaluations.
pub fn Potential(comptime T: type) type {
    return union(enum) {
        file: File(T),
        harmonic: Harmonic(T),
        henon_heiles: HenonHeiles(T),
        jahn_teller: JahnTeller(T),
        time_linear: TimeLinear(T),
        tully_1: Tully1(T),
        tully_2: Tully2(T),
        tully_3: Tully3(T),
        lvc: Lvc(T),
        custom: Custom(T),

        /// Initializes the selected potential energy surface based on configuration options.
        pub fn init(io: std.Io, options: Options, allocator: Allocator) !@This() {
            return switch (options) {
                .file => |f| .{ .file = try File(T).init(f.ndim, f.path, io, allocator) },
                .harmonic => |f| .{ .harmonic = Harmonic(T).init(f.k) },
                .henon_heiles => |f| .{ .henon_heiles = HenonHeiles(T).init(f.k, f.l) },
                .jahn_teller => |f| .{ .jahn_teller = JahnTeller(T).init(f.k, f.g) },
                .time_linear => |f| .{ .time_linear = TimeLinear(T).init(f.a, f.g) },
                .tully_1 => |f| .{ .tully_1 = Tully1(T).init(f.A, f.B, f.C, f.D) },
                .tully_2 => |f| .{ .tully_2 = Tully2(T).init(f.A, f.B, f.C, f.D, f.E0) },
                .tully_3 => |f| .{ .tully_3 = Tully3(T).init(f.A, f.B, f.C) },
                .lvc => |f| .{ .lvc = Lvc(T).init(f.frequencies, f.excitation_energies, f.kappa, f.lambda) },
                .custom => |f| .{ .custom = try Custom(T).init(f.ndim, f.matrix, f.time_dependent, allocator) },
            };
        }

        /// Deallocates memory associated with the potential energy surface representation.
        pub fn deinit(self: *@This(), allocator: Allocator) void {
            switch (self.*) {
                inline .file, .custom => |*pot| pot.deinit(allocator),

                inline else => {},
            }
        }

        /// Evaluates the potential energy matrix elements at coordinates r and time t.
        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, t: U) void {
            std.debug.assert(V.len == self.nstate() * self.nstate());

            std.debug.assert(r.len == self.ndim());

            switch (self) {
                inline else => |field| field.eval(U, V, r, t),
            }
        }

        /// Evaluates the potential energy matrix elements for a batch of coordinate coordinates.
        pub fn evalBatch(self: @This(), comptime U: type, V: *Matrix(U), r: Matrix(U), t: T) void {
            const t_val = Value(U).fromFloat(t);

            switch (self) {
                inline else => |field| {
                    for (0..r.nrow()) |i| field.eval(U, V.rowSlice(i), r.rowSlice(i), t_val.val);
                },
            }
        }

        /// Returns true if the potential has explicit time dependence.
        pub fn isTd(self: @This()) bool {
            return switch (self) {
                inline else => |field| field.isTd(),
            };
        }

        /// Returns the number of spatial dimensions of the coordinates.
        pub fn ndim(self: @This()) usize {
            return switch (self) {
                inline else => |field| field.ndim(),
            };
        }

        /// Returns the number of electronic states in the potential representation.
        pub fn nstate(self: @This()) usize {
            return switch (self) {
                inline else => |field| field.nstate(),
            };
        }
    };
}

/// Returns a harmonic potential energy surface type parameterized by force constants k.
fn Harmonic(comptime T: type) type {
    return struct {
        k: []const T,

        /// Initializes a harmonic potential with specified force constants.
        pub fn init(k: []const T) @This() {
            return .{ .k = k };
        }

        /// Evaluates the harmonic potential energy: V = 0.5 * sum( k_i * r_i^2 ).
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

/// Returns a Henon-Heiles potential energy surface type parameterized by frequency and anharmonicity.
fn HenonHeiles(comptime T: type) type {
    return struct {
        k: T,
        l: T,

        /// Initializes a Henon-Heiles potential with frequency and anharmonicity parameters.
        pub fn init(k: T, l: T) @This() {
            return .{ .k = k, .l = l };
        }

        /// Evaluates the Henon-Heiles potential energy: V = 0.5 * omg^2 * (x^2 + y^2) + lmb * (x^2 * y - y^3 / 3).
        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, _: U) void {
            std.debug.assert(r.len == 2);

            const x = Value(U).init(r[0]);
            const y = Value(U).init(r[1]);

            const k = Value(U).fromFloat(self.k);
            const l = Value(U).fromFloat(self.l);

            const V0 = k.muls(0.5).mul(x.mul(x).add(y.mul(y)));
            const V1 = l.mul(x.mul(x).mul(y).sub(y.mul(y).mul(y).divs(3)));

            V[0] = V0.add(V1).val;
        }

        /// Returns false as this potential is independent of time.
        pub fn isTd(_: @This()) bool {
            return false;
        }

        pub fn ndim(_: @This()) usize {
            return 2;
        }

        pub fn nstate(_: @This()) usize {
            return 1;
        }
    };
}

/// Returns a Jahn-Teller potential energy surface type parameterized by frequency and coupling strength.
fn JahnTeller(comptime T: type) type {
    return struct {
        k: T,
        g: T,

        /// Initializes a Jahn-Teller potential with frequency and coupling strength parameters.
        pub fn init(k: T, g: T) @This() {
            return .{ .k = k, .g = g };
        }

        /// Evaluates the Jahn-Teller potential energy matrix: V = 0.5 * k * (x^2 + y^2) * I + g * [[x, y], [y, -x]].
        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, _: U) void {
            std.debug.assert(r.len == 2);

            const x = Value(U).init(r[0]);
            const y = Value(U).init(r[1]);

            const k = Value(U).fromFloat(self.k);
            const g = Value(U).fromFloat(self.g);

            const V00 = k.muls(0.5).mul(x.mul(x).add(y.mul(y))).add(g.mul(x));
            const V01 = g.mul(y);
            const V11 = k.muls(0.5).mul(x.mul(x).add(y.mul(y))).sub(g.mul(x));

            V[0] = V00.val;
            V[1] = V01.val;
            V[2] = V01.val;
            V[3] = V11.val;
        }

        pub fn isTd(_: @This()) bool {
            return false;
        }

        pub fn ndim(_: @This()) usize {
            return 2;
        }

        pub fn nstate(_: @This()) usize {
            return 2;
        }
    };
}

/// Returns a time-dependent linear two-state coupling potential.
fn TimeLinear(comptime T: type) type {
    return struct {
        a: T,
        g: T,

        /// Initializes the linear coupling potential with slope and coupling parameters.
        pub fn init(a: T, g: T) @This() {
            return .{ .a = a, .g = g };
        }

        /// Evaluates the two-state time-dependent linear coupling potential matrix.
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

/// Returns a Tully 1 simple avoided crossing potential model type.
fn Tully1(comptime T: type) type {
    return struct {
        A: T,
        B: T,
        C: T,
        D: T,

        /// Initializes Tully 1 model parameters representing simple avoided crossing.
        pub fn init(A: T, B: T, C: T, D: T) @This() {
            return .{ .A = A, .B = B, .C = C, .D = D };
        }

        /// Evaluates the Tully simple avoided crossing two-state potential matrix.
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

/// Returns a Tully 2 dual avoided crossing potential model type.
fn Tully2(comptime T: type) type {
    return struct {
        A: T,
        B: T,
        C: T,
        D: T,

        E0: T,

        /// Initializes Tully 2 model parameters representing dual avoided crossing.
        pub fn init(A: T, B: T, C: T, D: T, E0: T) @This() {
            return .{ .A = A, .B = B, .C = C, .D = D, .E0 = E0 };
        }

        /// Evaluates the Tully dual avoided crossing two-state potential matrix.
        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, _: U) void {
            const r0 = Value(U).init(r[0]);

            const A = Value(U).fromFloat(self.A);
            const B = Value(U).fromFloat(self.B);
            const C = Value(U).fromFloat(self.C);
            const D = Value(U).fromFloat(self.D);

            const E0 = Value(U).fromFloat(self.E0);

            const V00 = Value(U).fromFloat(0);
            const V01 = r0.mul(r0).mul(D).neg().exp().mul(C);
            const V11 = r0.mul(r0).mul(B).neg().exp().mul(A).neg().add(E0);

            V[0] = V00.val;
            V[1] = V01.val;
            V[2] = V01.val;
            V[3] = V11.val;
        }

        /// Returns false as this potential is independent of time.
        pub fn isTd(_: @This()) bool {
            return false;
        }

        /// Returns the single nuclear coordinate dimension of the model.
        pub fn ndim(_: @This()) usize {
            return 1;
        }

        /// Returns the number of electronic states in this model system.
        pub fn nstate(_: @This()) usize {
            return 2;
        }
    };
}

/// Returns a Tully 3 extended coupling with reflection potential model type.
fn Tully3(comptime T: type) type {
    return struct {
        A: T,
        B: T,
        C: T,

        /// Initializes Tully 3 model parameters representing extended coupling.
        pub fn init(A: T, B: T, C: T) @This() {
            return .{ .A = A, .B = B, .C = C };
        }

        /// Evaluates the Tully extended coupling two-state potential matrix.
        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, _: U) void {
            const r0 = Value(U).init(r[0]);

            const A = Value(U).fromFloat(self.A);
            const B = Value(U).fromFloat(self.B);
            const C = Value(U).fromFloat(self.C);

            const V00 = A;
            const V01 = Value(U).fromFloat(1).sub(r0.abs().mul(C).neg().exp()).mul(r0.sign()).adds(1).mul(B);
            const V11 = A.neg();

            V[0] = V00.val;
            V[1] = V01.val;
            V[2] = V01.val;
            V[3] = V11.val;
        }

        /// Returns false as this potential is independent of time.
        pub fn isTd(_: @This()) bool {
            return false;
        }

        /// Returns the single nuclear coordinate dimension of the model.
        pub fn ndim(_: @This()) usize {
            return 1;
        }

        /// Returns the number of electronic states in this model system.
        pub fn nstate(_: @This()) usize {
            return 2;
        }
    };
}

/// Returns a custom potential type parsed from analytical string expressions.
fn Custom(comptime T: type) type {
    return struct {
        expressions: []const Expression(T),

        is_td: bool,
        r_vals: []T,

        /// Parses custom potential matrix expressions for a given coordinate dimension.
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

        /// Frees parsed expression resources and internal coordinate buffer.
        pub fn deinit(self: *@This(), allocator: Allocator) void {
            allocator.free(self.r_vals);

            for (self.expressions) |expr| {
                expr.deinit();
            }

            allocator.free(self.expressions);
        }

        /// Evaluates custom potentials by parsing expressions using ExprTk.
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

/// Returns a potential type interpolated from grid data stored in a file.
fn File(comptime T: type) type {
    return struct {
        U: Matrix(T),
        ndims: usize,

        /// Reads the potential grid data from the given file path.
        pub fn init(ndims: usize, path: []const u8, io: std.Io, gpa: Allocator) !@This() {
            return .{ .U = try readMatrix(T, io, path, gpa), .ndims = ndims };
        }

        /// Deallocates the stored potential grid matrix.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.U.deinit(gpa);
        }

        /// Evaluates the potential via multilinear interpolation of grid values.
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

/// Returns a linear vibronic coupling (LVC) potential energy surface type parameterized in eV and cm^-1.
fn Lvc(comptime T: type) type {
    return struct {
        frequencies: []const T,
        kap: []const []const T,

        excitation_energies: []const T,
        lmb: []const []const []const T,

        /// Initializes the linear vibronic coupling potential with frequencies, excitations, and coupling constants.
        pub fn init(frequencies: []const T, exc_en: []const T, kappa: []const []const T, lambda: []const []const []const T) @This() {
            std.debug.assert(kappa.len == exc_en.len);

            for (kappa) |row| {
                std.debug.assert(row.len == frequencies.len);
            }

            std.debug.assert(lambda.len == exc_en.len);

            for (lambda) |row| {
                std.debug.assert(row.len == exc_en.len);

                for (row) |col| {
                    std.debug.assert(col.len == frequencies.len);
                }
            }

            return .{ .frequencies = frequencies, .excitation_energies = exc_en, .kap = kappa, .lmb = lambda };
        }

        /// Evaluates the vibronic coupling potential matrix in Hartree at dimensionless normal coordinates r.
        pub fn eval(self: @This(), comptime U: type, V: []U, r: []const U, _: U) void {
            std.debug.assert(V.len == self.nstate() * self.nstate());

            var v0 = Value(U).fromFloat(0);

            for (0..self.frequencies.len) |i| {
                const q = Value(U).init(r[i]);

                v0 = v0.add(q.mul(q).muls(0.5 * self.frequencies[i] * CM2EV));
            }

            for (0..self.excitation_energies.len) |n| {
                var term = Value(U).fromFloat(self.excitation_energies[n]);

                for (0..self.frequencies.len) |i| {
                    term = term.add(Value(U).init(r[i]).muls(self.kap[n][i]));
                }

                V[n * self.nstate() + n] = v0.add(term).val;
            }

            for (0..self.excitation_energies.len) |m| for (m + 1..self.excitation_energies.len) |n| {
                var term = Value(U).fromFloat(0);

                for (0..self.frequencies.len) |i| {
                    term = term.add(Value(U).init(r[i]).muls(self.lmb[m][n][i]));
                }

                V[m * self.nstate() + n] = term.val;
                V[n * self.nstate() + m] = term.val;
            };

            for (0..V.len) |i| {
                V[i] = Value(U).init(V[i]).muls(EV2AU).val;
            }
        }

        /// Returns false as this vibronic potential has no explicit time dependence.
        pub fn isTd(_: @This()) bool {
            return false;
        }

        /// Returns the number of normal modes (nuclear degrees of freedom) in the potential.
        pub fn ndim(self: @This()) usize {
            return self.frequencies.len;
        }

        /// Returns the number of electronic states in the non-adiabatic potential representation.
        pub fn nstate(self: @This()) usize {
            return self.excitation_energies.len;
        }
    };
}

/// Performs multilinear interpolation on a rectangular coordinate grid.
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
                const r_k = if (comptime isDual(U)) r[k].val else r[k];

                if (grid.at(mid * stride, k) <= r_k) low = mid + 1 else high = mid;
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
