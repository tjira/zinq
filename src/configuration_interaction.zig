const std = @import("std");

const Allocator = std.mem.Allocator;

const HartreeFockOptions = @import("hartree_fock.zig").Options;
const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

const ao2mo_pp = @import("integral_transform.zig").ao2mo_pp;
const ao2mo_pppp = @import("integral_transform.zig").ao2mo_pppp;
const hartree_fock_run = @import("hartree_fock.zig").run;
const printf = @import("read_write.zig").printf;
const eighSlice = @import("linear_algebra.zig").eighSlice;

const ScalarDual = @import("dual.zig").ScalarDual;
const Value = @import("value.zig").Value;
const Integrals = @import("molecular_integrals.zig").Integrals;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;

// OPTIONS =============================================================================================================

pub const Options = struct {
    hartree_fock: HartreeFockOptions,

    excitations: []const u32 = &.{1},
};

// CONFIGURATION INTERACTION FUNCTIONS =================================================================================

pub fn slater(comptime T: type, A: []const usize, B: []const usize, H_MO: Matrix(T), g_MO: Tensor(T, 4), gpa: Allocator) !T {
    var diff_A = std.ArrayList(usize).empty;
    defer diff_A.deinit(gpa);

    var diff_B = std.ArrayList(usize).empty;
    defer diff_B.deinit(gpa);

    var common = std.ArrayList(usize).empty;
    defer common.deinit(gpa);

    var idx_a: usize = 0;
    var idx_b: usize = 0;

    while (idx_a < A.len and idx_b < B.len) switch (std.math.order(A[idx_a], B[idx_b])) {
        .lt => {
            try diff_A.append(gpa, A[idx_a]);

            idx_a += 1;
        },
        .gt => {
            try diff_B.append(gpa, B[idx_b]);

            idx_b += 1;
        },
        .eq => {
            try common.append(gpa, A[idx_a]);

            idx_a += 1;
            idx_b += 1;
        },
    };

    try diff_A.appendSlice(gpa, A[idx_a..]);
    try diff_B.appendSlice(gpa, B[idx_b..]);

    if (diff_A.items.len > 2) {
        return 0;
    }

    const s_A = signOfExcitations(A, diff_A.items);
    const s_B = signOfExcitations(B, diff_B.items);

    if (diff_A.items.len == 0) {
        var sum: T = 0;

        for (0..A.len) |i| {
            sum += H_MO.at(A[i], A[i]);
        }

        for (0..A.len) |i| for (i + 1..A.len) |j| {
            const idx_i = A[i];
            const idx_j = A[j];

            sum += g_MO.at(.{ idx_i, idx_j, idx_i, idx_j }) - g_MO.at(.{ idx_i, idx_j, idx_j, idx_i });
        };

        return sum;
    }

    if (diff_A.items.len == 1) {
        const m = diff_A.items[0];
        const p = diff_B.items[0];

        var term = H_MO.at(m, p);

        for (0..common.items.len) |i| {
            term += g_MO.at(.{ m, common.items[i], p, common.items[i] }) - g_MO.at(.{ m, common.items[i], common.items[i], p });
        }

        return @as(T, @floatFromInt(s_A * s_B)) * term;
    }

    if (diff_A.items.len == 2) {
        const m = diff_A.items[0];
        const n = diff_A.items[1];
        const p = diff_B.items[0];
        const q = diff_B.items[1];

        const term = g_MO.at(.{ m, n, p, q }) - g_MO.at(.{ m, n, q, p });

        return @as(T, @floatFromInt(s_A * s_B)) * term;
    }

    return 0;
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

fn signOfExcitations(S: []const usize, U: []const usize) i2 {
    var exp: usize = 0;

    for (0..U.len) |j| {
        exp += std.mem.indexOfScalar(usize, S, U[j]).? - j;
    }

    return if (exp % 2 == 0) 1 else -1;
}

// RESULT STRUCT =======================================================================================================

pub fn Result(comptime T: type) type {
    return struct {
        hartree_fock: HartreeFockResult(T),

        energy: []T,

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.hartree_fock.deinit(gpa);

            gpa.free(self.energy);
        }
    };
}

// RUN =================================================================================================================

pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator) !Result(T) {
    const generalized = opt.hartree_fock.generalized;

    const hf_opt = opt.hartree_fock;

    var hfres = try hartree_fock_run(T, io, hf_opt, log, gpa);
    errdefer hfres.deinit(gpa);

    const nbf = 1 * hfres.ints.sys.nbf;
    const nsp = 2 * hfres.ints.sys.nbf;

    var C_SO = try Matrix(T).initZero(nsp, nsp, gpa);
    defer C_SO.deinit(gpa);

    var H_SO = try Matrix(T).initZero(nsp, nsp, gpa);
    defer H_SO.deinit(gpa);

    var g_SO = try Tensor(T, 4).initZero(.{ nsp, nsp, nsp, nsp }, gpa);
    defer g_SO.deinit(gpa);

    if (generalized) {
        @memcpy(C_SO.data, hfres.C.data);

        for (0..nsp) |i| for (0..nsp) |j| {
            H_SO.ptr(i, j).* = hfres.ints.K.?.at(i, j) + hfres.ints.V.?.at(i, j);
        };

        for (0..nsp) |i| for (0..nsp) |j| for (0..nsp) |k| for (0..nsp) |l| {
            g_SO.ptr(.{ i, j, k, l }).* = hfres.ints.g.?.at(.{ i, j, k, l });
        };
    }

    if (!generalized) {
        for (0..nbf) |i| {
            const col_a = 2 * i + 0;
            const col_b = 2 * i + 1;

            for (0..nbf) |mu| {
                const indmu0, const indmu1 = .{ mu, mu + nbf };

                C_SO.ptr(indmu0, col_a).* = hfres.C.at(mu, i);
                C_SO.ptr(indmu1, col_b).* = hfres.C.at(mu, i);
            }

            for (0..nbf) |j| {
                const val = hfres.ints.K.?.at(i, j) + hfres.ints.V.?.at(i, j);

                const indi0, const indi1 = .{ i, i + nbf };
                const indj0, const indj1 = .{ j, j + nbf };

                H_SO.ptr(indi0, indj0).* = val;
                H_SO.ptr(indi1, indj1).* = val;
            }
        }

        for (0..nbf) |i| for (0..nbf) |j| for (0..nbf) |k| for (0..nbf) |l| {
            g_SO.ptr(.{ i, j, k, l }).* = hfres.ints.g.?.at(.{ i, j, k, l });

            g_SO.ptr(.{ i + nbf, j, k + nbf, l }).* = hfres.ints.g.?.at(.{ i, j, k, l });
            g_SO.ptr(.{ i, j + nbf, k, l + nbf }).* = hfres.ints.g.?.at(.{ i, j, k, l });

            g_SO.ptr(.{ i + nbf, j + nbf, k + nbf, l + nbf }).* = hfres.ints.g.?.at(.{ i, j, k, l });
        };
    }

    var H_MO = try Matrix(T).initZero(nsp, nsp, gpa);
    defer H_MO.deinit(gpa);

    ao2mo_pp(T, &H_MO, H_SO, C_SO);

    var g_MO = try Tensor(T, 4).initZero(.{ nsp, nsp, nsp, nsp }, gpa);
    defer g_MO.deinit(gpa);

    try ao2mo_pppp(T, &g_MO, g_SO, C_SO, gpa);

    var dets = try generateDets(hfres.ints.sys.nel, nsp, opt.excitations, gpa);

    defer {
        for (0..dets.items.len) |i| gpa.free(dets.items[i]);

        dets.deinit(gpa);
    }

    if (log) {
        try printf(io, "\nNUMBER OF CONFIGURATIONS: {d}\n", .{dets.items.len});
    }

    var H_CI = try Matrix(T).initZero(dets.items.len, dets.items.len, gpa);
    defer H_CI.deinit(gpa);

    for (0..dets.items.len) |i| for (i..dets.items.len) |j| {
        const val = try slater(T, dets.items[i], dets.items[j], H_MO, g_MO, gpa);

        H_CI.ptr(i, j).* = val;
        H_CI.ptr(j, i).* = val;
    };

    const W = try gpa.alloc(T, dets.items.len);
    errdefer gpa.free(W);

    var U = try Matrix(T).initZero(dets.items.len, dets.items.len, gpa);
    defer U.deinit(gpa);

    try eighSlice(T, W, U.data, H_CI.data);

    const VN = try hfres.ints.sys.nrep();

    for (0..W.len) |i| {
        W[i] += VN;
    }

    if (log) {
        try printf(io, "\nFINAL CI ENERGY: {d:.14} Eh\n", .{W[0]});
    }

    return Result(T){ .hartree_fock = hfres, .energy = W };
}
