const std = @import("std");

const Allocator = std.mem.Allocator;

const HartreeFockOptions = @import("hartree_fock.zig").Options;
const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

const ao2mo_pp = @import("integral_transform.zig").ao2mo_pp;
const ao2mo_pppp = @import("integral_transform.zig").ao2mo_pppp;
const ao2so_coef = @import("integral_transform.zig").ao2so_coef;
const ao2so_pp = @import("integral_transform.zig").ao2so_pp;
const ao2so_pppp = @import("integral_transform.zig").ao2so_pppp;
const eighSlice = @import("linear_algebra.zig").eighSlice;
const hartree_fock_run = @import("hartree_fock.zig").run;
const nuclearRepulsionGradient = @import("hartree_fock.zig").nuclearRepulsionGradient;
const primType = @import("value.zig").primType;
const printf = @import("read_write.zig").printf;

const ScalarDual = @import("dual.zig").ScalarDual;
const Value = @import("value.zig").Value;

pub const Options = struct {
    hartree_fock: HartreeFockOptions,

    excitations: []const u32 = &.{ 1, 2 },

    gradient: ?union(enum) {
        analytic: struct {
            state: u32 = 0,
        },
    } = null,
};

pub fn Result(comptime T: type) type {
    return struct {
        hartree_fock: HartreeFockResult(T),

        C: Matrix(T),

        energy: []T,

        gradient: []Matrix(T) = &.{},

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.hartree_fock.deinit(gpa);

            self.C.deinit(gpa);

            gpa.free(self.energy);

            for (0..self.gradient.len) |i| {
                self.gradient[i].deinit(gpa);
            }

            gpa.free(self.gradient);
        }
    };
}

pub fn generateDets(nel: usize, nsp: usize, excitations: []const u32, gpa: Allocator) !std.ArrayList([]const usize) {
    var dets = std.ArrayList([]const usize).empty;

    errdefer {
        for (0..dets.items.len) |i| gpa.free(dets.items[i]);

        dets.deinit(gpa);
    }

    {
        const ref_det = try gpa.alloc(usize, nel);

        for (0..nel) |i| {
            ref_det[i] = i;
        }

        try dets.append(gpa, ref_det);
    }

    for (0..excitations.len) |i| {
        const k = excitations[i];

        if (k == 0 or k > nel or k > nsp - nel) continue;

        var occ_combs = try generateCombinations(nel, k, 0, gpa);

        defer {
            for (0..occ_combs.items.len) |j| gpa.free(occ_combs.items[j]);

            occ_combs.deinit(gpa);
        }

        var vir_combs = try generateCombinations(nsp - nel, k, nel, gpa);

        defer {
            for (0..vir_combs.items.len) |j| gpa.free(vir_combs.items[j]);

            vir_combs.deinit(gpa);
        }

        for (0..occ_combs.items.len) |j| for (0..vir_combs.items.len) |l| {
            var list = try std.ArrayList(usize).initCapacity(gpa, nel);
            errdefer list.deinit(gpa);

            for (0..nel) |p| if (std.mem.indexOfScalar(usize, occ_combs.items[j], p) == null) {
                list.appendAssumeCapacity(p);
            };

            list.appendSliceAssumeCapacity(vir_combs.items[l]);

            std.mem.sort(usize, list.items, {}, std.sort.asc(usize));

            try dets.append(gpa, try list.toOwnedSlice(gpa));
        };
    }

    return dets;
}

pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator) !Result(T) {
    const generalized = opt.hartree_fock.generalized;

    var hf_opt = opt.hartree_fock;

    if (opt.gradient != null) {
        hf_opt.gradient = .{ .analytic = .{} };

        if (hf_opt.response == null) {
            hf_opt.response = .{};
        }
    }

    var timer = std.Io.Timestamp.now(io, .real);

    var hfres = try hartree_fock_run(T, io, hf_opt, log, gpa);
    errdefer hfres.deinit(gpa);

    if (log) {
        try printf(io, "\nHARTREE-FOCK CYCLE TIME: {f}\n\nCI EXCITATIONS CONSIDERED: [", .{timer.untilNow(io, .real)});

        for (0..opt.excitations.len) |i| {
            try printf(io, "{d}{s}", .{ opt.excitations[i], if (i == opt.excitations.len - 1) "]\n" else ", " });
        }

        try printf(io, "\nORBITAL TRANSFORMS TIME:", .{});
    }

    timer = std.Io.Timestamp.now(io, .real);

    var H_MS, var g_MS = try transformInts(T, hfres.C, hfres.ints.H.?, hfres.ints.g.?, generalized, gpa);

    defer {
        H_MS.deinit(gpa);
        g_MS.deinit(gpa);
    }

    if (log) {
        try printf(io, " {f}\nSLATER DETERMINANT TIME:", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    var dets = try generateDets(hfres.ints.sys.nel, 2 * hfres.ints.sys.nbf, opt.excitations, gpa);

    defer {
        for (0..dets.items.len) |i| gpa.free(dets.items[i]);

        dets.deinit(gpa);
    }

    if (log) {
        try printf(io, " {f}\n\nNUMBER OF DETERMINANTS: {d}\n\nCI MATRIX ASSEMBLY TIME:", .{ timer.untilNow(io, .real), dets.items.len });
    }

    timer = std.Io.Timestamp.now(io, .real);

    var H_CI = try Matrix(T).init(dets.items.len, dets.items.len, gpa);
    defer H_CI.deinit(gpa);

    for (0..dets.items.len) |i| for (i..dets.items.len) |j| {
        const val = slater(T, dets.items[i], dets.items[j], H_MS, g_MS);

        H_CI.ptr(i, j).* = val;
        H_CI.ptr(j, i).* = val;
    };

    if (log) {
        try printf(io, " {f}\nCI DIAGONALIZATION TIME:", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    var E, var C = try solveEigenvalueProblem(T, H_CI, try hfres.ints.sys.nrep(), gpa);

    if (log) {
        try printf(io, " {f}\n", .{timer.untilNow(io, .real)});
    }

    errdefer {
        E.deinit(gpa);
        C.deinit(gpa);
    }

    if (log) {
        try printf(io, "\nFINAL CI ENERGY: {d:.14} Eh\n", .{E.at(0)});
    }

    var grad = try gpa.alloc(Matrix(T), if (opt.gradient) |_| 1 else 0);
    errdefer if (opt.gradient) |_| gpa.free(grad);

    if (opt.gradient) |gradopt| {
        grad[0] = try gradient(T, hfres, C, dets, gradopt.analytic.state, gpa);
    }

    errdefer {
        if (opt.gradient) |_| grad[0].deinit(gpa);
    }

    if (log and opt.gradient != null) {
        try printf(io, "\nCI STATE {d} NUCLEAR ENERGY GRADIENT\n", .{opt.gradient.?.analytic.state});

        for (0..grad[0].nrow()) |j| for (0..grad[0].ncol()) |k| {
            try printf(io, "{d:20.14}{s}", .{ grad[0].at(j, k), if (k == 2) "\n" else " " });
        };
    }

    return Result(T){ .hartree_fock = hfres, .energy = E.data, .C = C, .gradient = grad };
}

pub fn slater(comptime T: type, A: []const usize, B: []const usize, H_MS: Matrix(T), g_MS: Tensor(T, 4)) T {
    var diff_A: [2]usize = undefined;
    var diff_B: [2]usize = undefined;

    var diff_a_len: usize = 0;
    var diff_b_len: usize = 0;

    var idx_a: usize = 0;
    var idx_b: usize = 0;

    while (idx_a < A.len and idx_b < B.len) {
        switch (std.math.order(A[idx_a], B[idx_b])) {
            .lt => {
                if (diff_a_len == 2) {
                    return Value(T).fromFloat(0).val;
                }

                diff_A[diff_a_len] = A[idx_a];

                diff_a_len, idx_a = .{ diff_a_len + 1, idx_a + 1 };
            },
            .gt => {
                if (diff_b_len == 2) {
                    return Value(T).fromFloat(0).val;
                }

                diff_B[diff_b_len] = B[idx_b];

                diff_b_len, idx_b = .{ diff_b_len + 1, idx_b + 1 };
            },
            .eq => {
                idx_a += 1;
                idx_b += 1;
            },
        }
    }

    while (idx_a < A.len) : (idx_a += 1) {
        if (diff_a_len == 2) {
            return Value(T).fromFloat(0).val;
        }

        diff_A[diff_a_len] = A[idx_a];

        diff_a_len += 1;
    }

    while (idx_b < B.len) : (idx_b += 1) {
        if (diff_b_len == 2) {
            return Value(T).fromFloat(0).val;
        }

        diff_B[diff_b_len] = B[idx_b];

        diff_b_len += 1;
    }

    const s_A = signOfExcitations(A, diff_A[0..diff_a_len]);
    const s_B = signOfExcitations(B, diff_B[0..diff_b_len]);

    if (diff_a_len == 0) {
        var sum = Value(T).fromFloat(0);

        for (0..A.len) |i| {
            sum = sum.add(Value(T).init(H_MS.at(A[i], A[i])));
        }

        for (0..A.len) |i| for (i + 1..A.len) |j| {
            const idx_i = A[i];
            const idx_j = A[j];

            const term1 = Value(T).init(g_MS.at(.{ idx_i, idx_j, idx_i, idx_j }));
            const term2 = Value(T).init(g_MS.at(.{ idx_i, idx_j, idx_j, idx_i }));

            sum = sum.add(term1.sub(term2));
        };

        return sum.val;
    }

    if (diff_a_len == 1) {
        const m = diff_A[0];
        const p = diff_B[0];

        var term = Value(T).init(H_MS.at(m, p));

        for (0..A.len) |i| if (A[i] != m) {
            const term1 = Value(T).init(g_MS.at(.{ m, A[i], p, A[i] }));
            const term2 = Value(T).init(g_MS.at(.{ m, A[i], A[i], p }));

            term = term.add(term1.sub(term2));
        };

        return term.muls(@as(primType(T), @floatFromInt(s_A * s_B))).val;
    }

    if (diff_a_len == 2) {
        const m = diff_A[0];
        const n = diff_A[1];
        const p = diff_B[0];
        const q = diff_B[1];

        const term1 = Value(T).init(g_MS.at(.{ m, n, p, q }));
        const term2 = Value(T).init(g_MS.at(.{ m, n, q, p }));

        return term1.sub(term2).muls(@as(primType(T), @floatFromInt(s_A * s_B))).val;
    }

    return Value(T).fromFloat(0).val;
}

fn generateCombinations(n: usize, k: usize, offset: usize, gpa: Allocator) !std.ArrayList([]const usize) {
    var results: std.ArrayList([]const usize) = .empty;

    errdefer {
        for (0..results.items.len) |i| gpa.free(results.items[i]);

        results.deinit(gpa);
    }

    if (k > n) return results;

    if (k == 0) {
        try results.append(gpa, try gpa.alloc(usize, 0));

        return results;
    }

    var current = try gpa.alloc(usize, k);
    defer gpa.free(current);

    for (0..k) |i| {
        current[i] = i;
    }

    while (true) {
        const copy = try gpa.alloc(usize, k);
        errdefer gpa.free(copy);

        for (0..k) |i| {
            copy[i] = current[i] + offset;
        }

        try results.append(gpa, copy);

        var i: usize = k - 1;

        while (current[i] == n - k + i) : (i -= 1) {
            if (i == 0) return results;
        }

        current[i] += 1;

        for (i + 1..k) |j| {
            current[j] = current[j - 1] + 1;
        }
    }

    return results;
}

fn gradient(comptime T: type, hfres: HartreeFockResult(T), C: Matrix(T), dets: std.ArrayList([]const usize), state: usize, gpa: Allocator) !Matrix(T) {
    const generalized = hfres.ints.sys.nbf != hfres.C.shape[0];

    var grad = try nuclearRepulsionGradient(T, hfres.ints.sys, gpa);
    errdefer grad.deinit(gpa);

    var H_d = try Matrix(ScalarDual(T)).init(hfres.ints.H.?.nrow(), hfres.ints.H.?.ncol(), gpa);
    defer H_d.deinit(gpa);

    var g_d = try Tensor(ScalarDual(T), 4).init(hfres.ints.g.?.shape, gpa);
    defer g_d.deinit(gpa);

    var C_d = try Matrix(ScalarDual(T)).init(hfres.C.nrow(), hfres.C.ncol(), gpa);
    defer C_d.deinit(gpa);

    for (0..3 * hfres.ints.sys.atoms.len) |i| {
        for (0..H_d.nrow()) |mu| for (0..H_d.ncol()) |nu| {
            H_d.ptr(mu, nu).* = ScalarDual(T).init(hfres.ints.H.?.at(mu, nu), hfres.ints.dH.?.at(.{ i, mu, nu }));
        };

        for (0..g_d.shape[0]) |mu| for (0..g_d.shape[1]) |lambda| for (0..g_d.shape[2]) |nu| for (0..g_d.shape[3]) |sigma| {
            const val = hfres.ints.g.?.at(.{ mu, lambda, nu, sigma });

            g_d.ptr(.{ mu, lambda, nu, sigma }).* = ScalarDual(T).init(val, hfres.ints.dg.?.at(.{ i, mu, lambda, nu, sigma }));
        };

        for (0..C_d.nrow()) |mu| for (0..C_d.ncol()) |p| {
            C_d.ptr(mu, p).* = ScalarDual(T).init(hfres.C.at(mu, p), hfres.dC.?.at(.{ i, mu, p }));
        };

        var H_MS_d, var g_MS_d = try transformInts(ScalarDual(T), C_d, H_d, g_d, generalized, gpa);

        defer {
            H_MS_d.deinit(gpa);
            g_MS_d.deinit(gpa);
        }

        var H_CI_d = try Matrix(ScalarDual(T)).init(dets.items.len, dets.items.len, gpa);
        defer H_CI_d.deinit(gpa);

        for (0..dets.items.len) |a| for (a..dets.items.len) |b| {
            const val = slater(ScalarDual(T), dets.items[a], dets.items[b], H_MS_d, g_MS_d);

            H_CI_d.ptr(a, b).* = val;
            H_CI_d.ptr(b, a).* = val;
        };

        for (0..dets.items.len) |a| for (0..dets.items.len) |b| {
            grad.ptr(i / 3, i % 3).* += C.at(a, state) * C.at(b, state) * H_CI_d.at(a, b).der;
        };
    }

    return grad;
}

fn signOfExcitations(S: []const usize, U: []const usize) i2 {
    var exp: usize = 0;

    for (0..U.len) |j| {
        exp += std.mem.indexOfScalar(usize, S, U[j]).? - j;
    }

    return if (exp % 2 == 0) 1 else -1;
}

fn solveEigenvalueProblem(comptime T: type, H_CI: Matrix(T), VN: T, gpa: Allocator) !struct { Vector(T), Matrix(T) } {
    var E = try Vector(T).init(H_CI.shape[0], gpa);
    errdefer E.deinit(gpa);

    var C = try Matrix(T).initZero(H_CI.shape[0], H_CI.shape[0], gpa);
    errdefer C.deinit(gpa);

    try eighSlice(T, E.data, C.data, H_CI.data);

    for (0..E.length()) |i| {
        E.ptr(i).* += VN;
    }

    return .{ E, C };
}

fn transformInts(comptime T: type, C: Matrix(T), H: Matrix(T), g: Tensor(T, 4), generalized: bool, gpa: Allocator) !struct { Matrix(T), Tensor(T, 4) } {
    const nsp = if (generalized) C.shape[0] else 2 * C.shape[0];

    var C_SO = if (generalized) C else try Matrix(T).initZero(nsp, nsp, gpa);
    defer if (!generalized) C_SO.deinit(gpa);

    if (!generalized) {
        ao2so_coef(T, &C_SO, C);
    }

    var H_SO = if (generalized) H else try Matrix(T).initZero(nsp, nsp, gpa);
    defer if (!generalized) H_SO.deinit(gpa);

    if (!generalized) {
        ao2so_pp(T, &H_SO, H);
    }

    var g_SO = if (generalized) g else try Tensor(T, 4).initZero(.{ nsp, nsp, nsp, nsp }, gpa);
    defer if (!generalized) g_SO.deinit(gpa);

    if (!generalized) {
        ao2so_pppp(T, &g_SO, g);
    }

    var H_MS = try Matrix(T).init(nsp, nsp, gpa);
    errdefer H_MS.deinit(gpa);

    ao2mo_pp(T, &H_MS, H_SO, C_SO);

    var g_MS = try Tensor(T, 4).init(.{ nsp, nsp, nsp, nsp }, gpa);
    errdefer g_MS.deinit(gpa);

    try ao2mo_pppp(T, &g_MS, g_SO, C_SO, gpa);

    return .{ H_MS, g_MS };
}
