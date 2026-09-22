//! Calculates the total spin squared expectation value and spin contamination for open-shell wavefunctions.

const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;

const mm = @import("linear_algebra.zig").mm;
const printf = @import("read_write.zig").printf;

/// Evaluates the total spin squared expectation value <S^2> from the generalized density matrix.
pub fn calculateTotalSpin(comptime T: type, P: Matrix(T), S: Matrix(T), gpa: Allocator) !T {
    const nbf = if (S.shape[0] == P.shape[0]) S.shape[0] / 2 else S.shape[0];

    std.debug.assert(P.shape[0] == 2 * nbf);
    std.debug.assert(P.shape[1] == 2 * nbf);

    var P_a = try Matrix(T).init(nbf, nbf, gpa);
    defer P_a.deinit(gpa);

    var P_b = try Matrix(T).init(nbf, nbf, gpa);
    defer P_b.deinit(gpa);

    var S_spatial = try Matrix(T).init(nbf, nbf, gpa);
    defer S_spatial.deinit(gpa);

    for (0..nbf) |i| for (0..nbf) |j| {
        P_a.ptr(i, j).* = P.at(i + 0 * nbf, j + 0 * nbf);
        P_b.ptr(i, j).* = P.at(i + 1 * nbf, j + 1 * nbf);

        S_spatial.ptr(i, j).* = S.at(i, j);
    };

    var AS = try Matrix(T).init(nbf, nbf, gpa);
    defer AS.deinit(gpa);

    mm(T, &AS, P_a, S_spatial, 1, 0, false, false);

    var BS = try Matrix(T).init(nbf, nbf, gpa);
    defer BS.deinit(gpa);

    mm(T, &BS, P_b, S_spatial, 1, 0, false, false);

    var n_a: T = 0;
    var n_b: T = 0;

    var tr_psps: T = 0;

    for (0..nbf) |i| {
        n_a += AS.at(i, i);
        n_b += BS.at(i, i);

        for (0..nbf) |j| {
            tr_psps += AS.at(i, j) * BS.at(j, i);
        }
    }

    const sz = 0.5 * (n_a - n_b);

    return sz * (sz + 1) + n_b - tr_psps;
}

/// Formats and prints the calculated total spin expectation value and spin contamination to the output.
pub fn printTotalSpin(comptime T: type, io: std.Io, s2: T, multiplicity: u32, method_str: []const u8) !void {
    const s_exact = 0.5 * @as(T, @floatFromInt(multiplicity - 1));

    const s2_exact = s_exact * (s_exact + 1);

    const fmt = "\n{s} <S^2> (CALCULATED): {d:20.14}\n{s} <S^2> (EIGENVALUE): {d:20.14}\n{s} SPIN CONTAMINATION: {d:20.14}\n";

    try printf(io, fmt, .{ method_str, s2, method_str, s2_exact, method_str, s2 - s2_exact });
}
