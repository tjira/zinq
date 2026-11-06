//! Generates and compares electronic Slater determinants for complete active space and truncated CI expansions.

const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;

const calculateTotalSpin = @import("spin_analysis.zig").calculateTotalSpin;
const mo2ms_coef = @import("integral_transform.zig").mo2ms_coef;
const printf = @import("read_write.zig").printf;

/// Holds determinant index, expansion coefficient, and contribution weight for printing.
const DetContrib = struct { d_idx: usize, coeff: f64, weight: f64 };

/// Evaluates the total spin squared expectation value <S^2> of a Slater determinant with orthonormal spatial orbitals.
pub fn calculateDeterminantSpin(det: []const usize) f64 {
    var n_a: usize = 0;
    var n_b: usize = 0;

    var n_pair: usize = 0;

    for (det) |k| {
        if (k % 2 == 0) {
            n_a += 1;

            if (std.mem.indexOfScalar(usize, det, k + 1) != null) {
                n_pair += 1;
            }
        } else {
            n_b += 1;
        }
    }

    const n_a_f = @as(f64, @floatFromInt(n_a));
    const n_b_f = @as(f64, @floatFromInt(n_b));
    const n_p_f = @as(f64, @floatFromInt(n_pair));

    const ms = 0.5 * (n_a_f - n_b_f);

    return ms * (ms + 1.0) + (n_b_f - n_p_f);
}

/// Evaluates the total spin squared expectation value <S^2> of a single Slater determinant from its density matrix.
pub fn calculateDeterminantTotalSpin(comptime T: type, det: []const usize, C_ms: Matrix(T), S_ao: Matrix(T), gpa: Allocator) !T {
    const nsp = C_ms.shape[0];

    var P = try Matrix(T).initZero(nsp, nsp, gpa);
    defer P.deinit(gpa);

    for (det) |k| for (0..nsp) |mu| {
        const c_mu = C_ms.at(mu, k);

        if (c_mu == 0) continue;

        for (0..nsp) |nu| {
            P.ptr(mu, nu).* += c_mu * C_ms.at(nu, k);
        }
    };

    return try calculateTotalSpin(T, P, S_ao, gpa);
}

/// Compares two Slater determinants and extracts orbital differences when excitation degree is at most two.
pub fn findDeterminantDifferences(A: []const usize, B: []const usize) ?struct { a: [2]usize, b: [2]usize, len: usize } {
    var diff_A: [2]usize = undefined;
    var diff_B: [2]usize = undefined;

    var diff_a_len: usize = 0;
    var diff_b_len: usize = 0;

    var idx_a: usize = 0;
    var idx_b: usize = 0;

    while (idx_a < A.len and idx_b < B.len) {
        switch (std.math.order(A[idx_a], B[idx_b])) {
            .lt => {
                if (diff_a_len == 2) return null;

                diff_A[diff_a_len] = A[idx_a];

                diff_a_len, idx_a = .{ diff_a_len + 1, idx_a + 1 };
            },
            .gt => {
                if (diff_b_len == 2) return null;

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
        if (diff_a_len == 2) return null;

        diff_A[diff_a_len] = A[idx_a];

        diff_a_len += 1;
    }

    while (idx_b < B.len) : (idx_b += 1) {
        if (diff_b_len == 2) return null;

        diff_B[diff_b_len] = B[idx_b];

        diff_b_len += 1;
    }

    if (diff_a_len != diff_b_len) return null;

    return .{ .a = diff_A, .b = diff_B, .len = diff_a_len };
}

/// Generates complete active space determinants preserving closed-shell core occupations.
pub fn generateCasDets(ncore: usize, nact: usize, nact_el: usize, mult: u32, gpa: Allocator) !std.ArrayList([]const usize) {
    var dets = std.ArrayList([]const usize).empty;

    errdefer {
        for (0..dets.items.len) |i| gpa.free(dets.items[i]);

        dets.deinit(gpa);
    }

    if ((nact_el + mult - 1) % 2 != 0) return error.InvalidInput;

    const n_alpha = (nact_el + mult - 1) / 2;

    const n_beta = nact_el - n_alpha;

    if (n_alpha > nact or n_beta > nact) return error.InvalidInput;

    var alpha_combs = try generateCombinations(nact, n_alpha, 0, gpa);

    defer {
        for (0..alpha_combs.items.len) |i| gpa.free(alpha_combs.items[i]);

        alpha_combs.deinit(gpa);
    }

    var beta_combs = try generateCombinations(nact, n_beta, 0, gpa);

    defer {
        for (0..beta_combs.items.len) |i| gpa.free(beta_combs.items[i]);

        beta_combs.deinit(gpa);
    }

    const nel = 2 * ncore + nact_el;

    for (alpha_combs.items) |ac| for (beta_combs.items) |bc| {
        const det = try gpa.alloc(usize, nel);
        errdefer gpa.free(det);

        for (0..2 * ncore) |k| {
            det[k] = k;
        }

        for (ac, 0..) |a, k| {
            det[2 * ncore + k] = 2 * (ncore + a);
        }

        for (bc, 0..) |b, k| {
            det[2 * ncore + n_alpha + k] = 2 * (ncore + b) + 1;
        }

        std.mem.sort(usize, det, {}, std.sort.asc(usize));

        try dets.append(gpa, det);
    };

    return dets;
}

/// Generates all k-combinations from a set of size n, representing orbital indices.
pub fn generateCombinations(n: usize, k: usize, offset: usize, gpa: Allocator) !std.ArrayList([]const usize) {
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

/// Generates excited Slater determinants up to specified excitation levels from a Hartree-Fock reference.
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

/// Generates all unique upper-triangular state index pairs for nonadiabatic coupling evaluation.
pub fn generatePairs(n: usize, gpa: Allocator) ![][2]u32 {
    const pairs = try gpa.alloc([2]u32, n * (n - 1) / 2);

    var idx: usize = 0;

    for (0..n) |i| for (i + 1..n) |j| {
        pairs[idx] = .{ @intCast(i), @intCast(j) };

        idx += 1;
    };

    return pairs;
}

/// Prints leading electronic determinants, configuration interaction coefficients, and weights for each state.
pub fn printDeterminants(comptime T: type, io: std.Io, C_ci: Matrix(T), C_mo: Matrix(T), S_ao: Matrix(T), dets: []const []const usize, nstates: usize, ncore: usize, norbs: usize, gen: bool, gpa: Allocator) !void {
    const threshold: f64 = 0.01;

    const nsp = if (gen) C_mo.shape[0] else 2 * C_mo.shape[0];

    var C_ms = if (gen) C_mo else try Matrix(T).initZero(nsp, nsp, gpa);
    defer if (!gen) C_ms.deinit(gpa);

    if (!gen) {
        mo2ms_coef(T, &C_ms, C_mo);
    }

    for (0..nstates) |s| {
        var list = std.ArrayList(DetContrib).empty;
        defer list.deinit(gpa);

        var max_w: f64 = 0;
        var max_c: f64 = 0;

        var max_i: usize = 0;

        for (0..C_ci.nrow()) |d| {
            const c_val: f64 = @floatCast(C_ci.at(d, s));

            const weight: f64 = c_val * c_val;

            if (weight > max_w) {
                max_w, max_c, max_i = .{ weight, c_val, d };
            }

            if (weight >= threshold) {
                try list.append(gpa, .{ .d_idx = d, .coeff = c_val, .weight = weight });
            }
        }

        if (list.items.len == 0 and C_ci.nrow() > 0) {
            try list.append(gpa, .{ .d_idx = max_i, .coeff = max_c, .weight = max_w });
        }

        std.mem.sort(DetContrib, list.items, {}, struct {
            fn desc(_: void, a: DetContrib, b: DetContrib) bool {
                return a.weight > b.weight;
            }
        }.desc);

        if (nstates == 1) {
            try printf(io, "\nLEADING DETERMINANT CONTRIBUTIONS (|C|^2 >= 0.0100)\n", .{});
        }

        if (nstates > 1) {
            try printf(io, "\nSTATE {d:2} LEADING DETERMINANT CONTRIBUTIONS (|C|^2 >= 0.0100)\n", .{s});
        }

        try printf(io, "{s:6} {s:18} {s:14} {s:14} {s:8}\n", .{ "DET", "OCCUPATION", "COEF (C)", "WEIGHT (|C|^2)", "<S^2>" });

        for (list.items) |item| {
            const desc = try formatDeterminant(dets[item.d_idx], ncore, norbs, gen, gpa);
            defer gpa.free(desc);

            const wght = item.weight * 100.0;

            const s_val = try calculateDeterminantTotalSpin(T, dets[item.d_idx], C_ms, S_ao, gpa);

            try printf(io, "{d:6} {s:18} {d:14.6} {d:12.2} % {d:8.2}\n", .{ item.d_idx + 1, desc, item.coeff, wght, s_val });
        }
    }
}

/// Computes the permutation sign (+1 or -1) associated with mapping two determinants via orbital excitations.
pub fn signOfExcitations(S: []const usize, U: []const usize) i2 {
    var exp: usize = 0;

    for (0..U.len) |j| {
        exp += std.mem.indexOfScalar(usize, S, U[j]).? - j;
    }

    return if (exp % 2 == 0) 1 else -1;
}

/// Formats a Slater determinant into an allocated occupation string representation over specified molecular orbitals.
fn formatDeterminant(det: []const usize, ncore: usize, norbs: usize, gen: bool, gpa: Allocator) ![]u8 {
    const len = 2 * norbs + 3;

    var buf = try gpa.alloc(u8, len);
    errdefer gpa.free(buf);

    var pos: usize = 0;

    buf[pos] = '|';

    pos += 1;

    for (0..norbs) |u| {
        const p = ncore + u;

        buf[pos] = ' ';

        pos += 1;

        if (gen) {
            const has_p = std.mem.indexOfScalar(usize, det, p) != null;

            buf[pos] = if (has_p) '1' else '0';
        }

        if (!gen) {
            const has_a = std.mem.indexOfScalar(usize, det, 2 * p + 0) != null;
            const has_b = std.mem.indexOfScalar(usize, det, 2 * p + 1) != null;

            buf[pos] = if (has_a and has_b) '2' else if (has_a) 'a' else if (has_b) 'b' else '0';
        }

        pos += 1;
    }

    buf[pos] = ' ';

    pos += 1;

    buf[pos] = '>';

    return buf;
}
