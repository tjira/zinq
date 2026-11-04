//! Generates and compares electronic Slater determinants for complete active space and truncated CI expansions.

const std = @import("std");

const Allocator = std.mem.Allocator;

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
pub fn generateCasDets(ncore: usize, nas: usize, nas_el: usize, mult: u32, gpa: Allocator) !std.ArrayList([]const usize) {
    var dets = std.ArrayList([]const usize).empty;

    errdefer {
        for (0..dets.items.len) |i| gpa.free(dets.items[i]);

        dets.deinit(gpa);
    }

    if ((nas_el + mult - 1) % 2 != 0) return error.InvalidInput;

    const n_alpha = (nas_el + mult - 1) / 2;

    const n_beta = nas_el - n_alpha;

    if (n_alpha > nas or n_beta > nas) return error.InvalidInput;

    var alpha_combs = try generateCombinations(nas, n_alpha, 0, gpa);

    defer {
        for (0..alpha_combs.items.len) |i| gpa.free(alpha_combs.items[i]);

        alpha_combs.deinit(gpa);
    }

    var beta_combs = try generateCombinations(nas, n_beta, 0, gpa);

    defer {
        for (0..beta_combs.items.len) |i| gpa.free(beta_combs.items[i]);

        beta_combs.deinit(gpa);
    }

    const nel = 2 * ncore + nas_el;

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

/// Computes the permutation sign (+1 or -1) associated with mapping two determinants via orbital excitations.
pub fn signOfExcitations(S: []const usize, U: []const usize) i2 {
    var exp: usize = 0;

    for (0..U.len) |j| {
        exp += std.mem.indexOfScalar(usize, S, U[j]).? - j;
    }

    return if (exp % 2 == 0) 1 else -1;
}