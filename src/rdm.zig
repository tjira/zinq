//! Constructs active-space one- and two-particle reduced density matrices from configuration coefficients.

const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;

const findDeterminantDifferences = @import("determinant.zig").findDeterminantDifferences;
const signOfExcitations = @import("determinant.zig").signOfExcitations;

/// Builds active-space 1-particle and 2-particle reduced density matrices from CAS-CI coefficients.
pub fn buildActRdm(comptime T: type, dets: []const []const usize, C_ci: Matrix(T), state: usize, ncore: usize, nact: usize, gpa: Allocator) !struct { Matrix(T), Tensor(T, 4) } {
    var D_act = try Matrix(T).initZero(nact, nact, gpa);
    errdefer D_act.deinit(gpa);

    var d_act = try Tensor(T, 4).initZero(.{ nact, nact, nact, nact }, gpa);
    errdefer d_act.deinit(gpa);

    for (0..dets.len) |i| for (0..dets.len) |j| {
        const c_prod = C_ci.at(i, state) * C_ci.at(j, state);

        if (@abs(c_prod) < 1e-15) continue;

        const diff = findDeterminantDifferences(dets[i], dets[j]) orelse continue;

        const sA = signOfExcitations(dets[i], diff.a[0..diff.len]);
        const sB = signOfExcitations(dets[j], diff.b[0..diff.len]);

        const phase: T = if (diff.len > 0) @as(T, @floatFromInt(sA * sB)) else 1;

        switch (diff.len) {
            0 => accumulateRdmDiagonal(T, &D_act, &d_act, dets[i], c_prod, ncore, nact),
            1 => accumulateRdmSingle(T, &D_act, &d_act, dets[i], diff.a[0], diff.b[0], phase, c_prod, ncore, nact),
            2 => accumulateRdmDouble(T, &d_act, diff.a, diff.b, phase, c_prod, ncore, nact),
            else => {},
        }
    };

    return .{ D_act, d_act };
}

/// Assembles state-averaged 1- and 2-particle reduced density matrices across active configurations.
pub fn buildStateAveragedActRdm(comptime T: type, dets: []const []const usize, C_ci: Matrix(T), weights: []const T, ncore: usize, nact: usize, gpa: Allocator) !struct { Matrix(T), Tensor(T, 4) } {
    var D_act = try Matrix(T).initZero(nact, nact, gpa);
    errdefer D_act.deinit(gpa);

    var d_act = try Tensor(T, 4).initZero(.{ nact, nact, nact, nact }, gpa);
    errdefer d_act.deinit(gpa);

    for (0..weights.len) |I| {
        const w_I = weights[I];

        if (w_I == 0) continue;

        var D_I, var d_I = try buildActRdm(T, dets, C_ci, I, ncore, nact, gpa);

        defer {
            D_I.deinit(gpa);
            d_I.deinit(gpa);
        }

        for (0..D_act.data.len) |idx| {
            D_act.data[idx] += w_I * D_I.data[idx];
        }

        for (0..d_act.data.len) |idx| {
            d_act.data[idx] += w_I * d_I.data[idx];
        }
    }

    return .{ D_act, d_act };
}

/// Accumulates diagonal 1-RDM and 2-RDM contributions between identical configurations.
fn accumulateRdmDiagonal(comptime T: type, D_act: *Matrix(T), d_act: *Tensor(T, 4), det: []const usize, c_prod: T, ncore: usize, nact: usize) void {
    for (0..nact) |u| {
        const pa = 2 * (ncore + u) + 0;
        const pb = 2 * (ncore + u) + 1;

        const occ_a = std.mem.indexOfScalar(usize, det, pa) != null;
        const occ_b = std.mem.indexOfScalar(usize, det, pb) != null;

        if (occ_a) D_act.ptr(u, u).* += c_prod;
        if (occ_b) D_act.ptr(u, u).* += c_prod;

        if (occ_a and occ_b) d_act.ptr(.{ u, u, u, u }).* += 2 * c_prod;

        for (0..nact) |w| if (u != w) {
            const qa = 2 * (ncore + w) + 0;
            const qb = 2 * (ncore + w) + 1;

            const w_a = std.mem.indexOfScalar(usize, det, qa) != null;
            const w_b = std.mem.indexOfScalar(usize, det, qb) != null;

            const nu: T = (if (occ_a) @as(T, 1) else @as(T, 0)) + (if (occ_b) @as(T, 1) else @as(T, 0));

            const nw: T = (if (w_a) @as(T, 1) else @as(T, 0)) + (if (w_b) @as(T, 1) else @as(T, 0));

            d_act.ptr(.{ u, u, w, w }).* += nu * nw * c_prod;

            var exch: T = 0;

            if (occ_a and w_a) exch += 1;
            if (occ_b and w_b) exch += 1;

            d_act.ptr(.{ u, w, w, u }).* -= exch * c_prod;
        };
    }
}

/// Accumulates 2-RDM contributions between configurations differing by two spin-orbitals.
fn accumulateRdmDouble(comptime T: type, d_act: *Tensor(T, 4), diff_a: [2]usize, diff_b: [2]usize, phase: T, c_prod: T, ncore: usize, nact: usize) void {
    const p, const q = .{ diff_a[0], diff_a[1] };
    const r, const s = .{ diff_b[0], diff_b[1] };

    if (p >= 2 * ncore and q >= 2 * ncore and r >= 2 * ncore and s >= 2 * ncore) {
        const u = (p / 2) - ncore;
        const v = (q / 2) - ncore;
        const w = (r / 2) - ncore;
        const x = (s / 2) - ncore;

        if (u < nact and v < nact and w < nact and x < nact) {
            if ((p % 2 == r % 2) and (q % 2 == s % 2)) {
                d_act.ptr(.{ u, w, v, x }).* += c_prod * phase;
                d_act.ptr(.{ v, x, u, w }).* += c_prod * phase;
            }

            if ((p % 2 == s % 2) and (q % 2 == r % 2)) {
                d_act.ptr(.{ u, x, v, w }).* -= c_prod * phase;
                d_act.ptr(.{ v, w, u, x }).* -= c_prod * phase;
            }
        }
    }
}

/// Accumulates 1-RDM and 2-RDM contributions between configurations differing by one spin-orbital.
fn accumulateRdmSingle(comptime T: type, D_act: *Matrix(T), d_act: *Tensor(T, 4), det: []const usize, p: usize, q: usize, phase: T, c_prod: T, ncore: usize, nact: usize) void {
    if (p >= 2 * ncore and q >= 2 * ncore) {
        const u = (p / 2) - ncore;
        const v = (q / 2) - ncore;

        if (u < nact and v < nact and (p % 2 == q % 2)) {
            D_act.ptr(u, v).* += c_prod * phase;

            for (0..nact) |w| {
                const ka = 2 * (ncore + w) + 0;
                const kb = 2 * (ncore + w) + 1;

                const w_a = std.mem.indexOfScalar(usize, det, ka) != null;
                const w_b = std.mem.indexOfScalar(usize, det, kb) != null;

                const occ_ka = w_a and ka != p and ka != q;
                const occ_kb = w_b and kb != p and kb != q;

                const nw: T = (if (occ_ka) @as(T, 1) else @as(T, 0)) + (if (occ_kb) @as(T, 1) else @as(T, 0));

                d_act.ptr(.{ u, v, w, w }).* += nw * c_prod * phase;
                d_act.ptr(.{ w, w, u, v }).* += nw * c_prod * phase;

                var exch: T = 0;

                if ((p % 2 == 0) and occ_ka) exch += 1;
                if ((p % 2 == 1) and occ_kb) exch += 1;

                d_act.ptr(.{ u, w, w, v }).* -= exch * c_prod * phase;
                d_act.ptr(.{ w, v, u, w }).* -= exch * c_prod * phase;
            }
        }
    }
}
