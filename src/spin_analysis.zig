//! Calculates the total spin squared expectation value and spin contamination for open-shell wavefunctions.

const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;

const gemm = @import("linear_algebra.zig").gemm;
const printf = @import("read_write.zig").printf;

/// Holds spin-resolved blocks of a generalized single-particle density matrix.
pub fn SpinBlocks(comptime T: type) type {
    return struct {
        a: Matrix(T),
        b: Matrix(T),

        ab: Matrix(T),
        ba: Matrix(T),

        /// Extracts alpha, beta, and cross spin blocks from a generalized matrix.
        pub fn init(P: Matrix(T), nbf: usize, gpa: Allocator) !@This() {
            var a = try Matrix(T).init(nbf, nbf, gpa);
            errdefer a.deinit(gpa);

            var b = try Matrix(T).init(nbf, nbf, gpa);
            errdefer b.deinit(gpa);

            var ab = try Matrix(T).init(nbf, nbf, gpa);
            errdefer ab.deinit(gpa);

            var ba = try Matrix(T).init(nbf, nbf, gpa);
            errdefer ba.deinit(gpa);

            for (0..nbf) |i| for (0..nbf) |j| {
                a.ptr(i, j).* = P.at(i + 0 * nbf, j + 0 * nbf);
                b.ptr(i, j).* = P.at(i + 1 * nbf, j + 1 * nbf);

                ab.ptr(i, j).* = P.at(i + 0 * nbf, j + 1 * nbf);
                ba.ptr(i, j).* = P.at(i + 1 * nbf, j + 0 * nbf);
            };

            return .{ .a = a, .b = b, .ab = ab, .ba = ba };
        }

        /// Deallocates memory associated with the four spin-block matrices.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.a.deinit(gpa);
            self.b.deinit(gpa);

            self.ab.deinit(gpa);
            self.ba.deinit(gpa);
        }
    };
}

/// Evaluates the total spin squared expectation value <S^2> from the generalized density matrix.
pub fn calculateTotalSpin(comptime T: type, P: Matrix(T), S: Matrix(T), gpa: Allocator) !T {
    const nbf = if (S.shape[0] == P.shape[0]) S.shape[0] / 2 else S.shape[0];

    std.debug.assert(P.shape[0] == 2 * nbf);
    std.debug.assert(P.shape[1] == 2 * nbf);

    var blocks = try SpinBlocks(T).init(P, nbf, gpa);
    defer blocks.deinit(gpa);

    var S_spatial = try Matrix(T).init(nbf, nbf, gpa);
    defer S_spatial.deinit(gpa);

    for (0..nbf) |i| for (0..nbf) |j| {
        S_spatial.ptr(i, j).* = S.at(i, j);
    };

    var AS = try Matrix(T).init(nbf, nbf, gpa);
    defer AS.deinit(gpa);

    var BS = try Matrix(T).init(nbf, nbf, gpa);
    defer BS.deinit(gpa);

    var ABS = try Matrix(T).init(nbf, nbf, gpa);
    defer ABS.deinit(gpa);

    var BAS = try Matrix(T).init(nbf, nbf, gpa);
    defer BAS.deinit(gpa);

    gemm(T, &AS, blocks.a, S_spatial, 1, 0, false, false);
    gemm(T, &BS, blocks.b, S_spatial, 1, 0, false, false);

    gemm(T, &ABS, blocks.ab, S_spatial, 1, 0, false, false);
    gemm(T, &BAS, blocks.ba, S_spatial, 1, 0, false, false);

    var n_a: T = 0;
    var n_b: T = 0;

    var tr_abs: T = 0;
    var tr_bas: T = 0;

    var tr_asbs: T = 0;

    var tr_absbas: T = 0;

    for (0..nbf) |i| {
        n_a += AS.at(i, i);
        n_b += BS.at(i, i);

        tr_abs += ABS.at(i, i);
        tr_bas += BAS.at(i, i);

        for (0..nbf) |j| {
            tr_asbs += AS.at(i, j) * BS.at(j, i);

            tr_absbas += ABS.at(i, j) * BAS.at(j, i);
        }
    }

    const sx, const sz = .{ 0.5 * (tr_abs + tr_bas), 0.5 * (n_a - n_b) };

    return (sx * sx + sz * sz) + 0.5 * (n_a + n_b) - tr_asbs + tr_absbas;
}

/// Formats and prints the calculated total spin expectation value and spin contamination to the output.
pub fn printTotalSpin(comptime T: type, io: std.Io, s2: T, multiplicity: u32, method_str: []const u8) !void {
    const s_exact = 0.5 * @as(T, @floatFromInt(multiplicity - 1));

    const s2_exact = s_exact * (s_exact + 1);

    const fmt = "\n{s} <S^2> (CALCULATED): {d:20.14}\n{s} <S^2> (EIGENVALUE): {d:20.14}\n{s} SPIN CONTAMINATION: {d:20.14}\n";

    try printf(io, fmt, .{ method_str, s2, method_str, s2_exact, method_str, s2 - s2_exact });
}
