//! Implements Mulliken and Löwdin population analysis and Mayer and Wiberg bond order calculations.

const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const SpinBlocks = @import("spin_analysis.zig").SpinBlocks;
const Vector = @import("tensor.zig").Vector;

const dot = @import("linear_algebra.zig").dot;
const gemm = @import("linear_algebra.zig").gemm;
const getSymbol = @import("constant.zig").getSymbol;
const printf = @import("read_write.zig").printf;
const syevd = @import("linear_algebra.zig").syevd;

/// Computes Löwdin net atomic charges using symmetric orthogonalization of the overlap matrix.
pub fn lowdin(comptime T: type, sys: MolecularSystem(T), P: Matrix(T), S: Matrix(T), gpa: Allocator) !Vector(T) {
    const n = S.nrow();

    var Shalf = try matrixSquareRoot(T, S, gpa);
    defer Shalf.deinit(gpa);

    var SP = try Matrix(T).init(n, n, gpa);
    defer SP.deinit(gpa);

    gemm(T, &SP, Shalf, P, 1.0, 0.0, false, false);

    var net_populations = try gpa.alloc(T, sys.atoms.len);
    defer gpa.free(net_populations);

    @memset(net_populations, 0);

    for (0..n) |mu| {
        net_populations[@intCast(sys.bf2at[mu % sys.nbf])] += dot(T, SP.row(mu), Shalf.row(mu));
    }

    var charges = try Vector(T).init(sys.atoms.len, gpa);
    errdefer charges.deinit(gpa);

    for (0..sys.atoms.len) |A| {
        charges.data[A] = @as(T, @floatFromInt(sys.atoms[A])) - net_populations[A];
    }

    return charges;
}

/// Computes Mayer bond orders between all pairs of atoms from the density and overlap matrices.
pub fn mayer(comptime T: type, sys: MolecularSystem(T), P: Matrix(T), S: Matrix(T), gpa: Allocator) !Matrix(T) {
    const nbf = sys.nbf;

    const is_gen = (P.shape[0] == 2 * nbf);

    var bo = try Matrix(T).initZero(sys.atoms.len, sys.atoms.len, gpa);
    errdefer bo.deinit(gpa);

    var S_spatial = try Matrix(T).init(nbf, nbf, gpa);
    defer S_spatial.deinit(gpa);

    for (0..nbf) |mu| for (0..nbf) |nu| {
        S_spatial.ptr(mu, nu).* = S.at(mu, nu);
    };

    if (!is_gen) {
        var PS = try Matrix(T).init(nbf, nbf, gpa);
        defer PS.deinit(gpa);

        gemm(T, &PS, P, S_spatial, 1, 0, false, false);

        for (0..nbf) |mu| {
            const A: usize = @intCast(sys.bf2at[mu]);

            for (0..nbf) |nu| {
                const B: usize = @intCast(sys.bf2at[nu]);

                if (A != B) {
                    bo.ptr(A, B).* += PS.at(mu, nu) * PS.at(nu, mu);
                }
            }
        }
    }

    if (is_gen) {
        var blocks = try SpinBlocks(T).init(P, nbf, gpa);
        defer blocks.deinit(gpa);

        var PS_a = try Matrix(T).init(nbf, nbf, gpa);
        defer PS_a.deinit(gpa);

        var PS_b = try Matrix(T).init(nbf, nbf, gpa);
        defer PS_b.deinit(gpa);

        var PS_ab = try Matrix(T).init(nbf, nbf, gpa);
        defer PS_ab.deinit(gpa);

        gemm(T, &PS_a, blocks.a, S_spatial, 1, 0, false, false);
        gemm(T, &PS_b, blocks.b, S_spatial, 1, 0, false, false);

        gemm(T, &PS_ab, blocks.ab, S_spatial, 1, 0, false, false);

        for (0..nbf) |mu| {
            const A: usize = @intCast(sys.bf2at[mu]);

            for (0..nbf) |nu| {
                const B: usize = @intCast(sys.bf2at[nu]);

                if (A != B) {
                    const term_a = PS_a.at(mu, nu) * PS_a.at(nu, mu);
                    const term_b = PS_b.at(mu, nu) * PS_b.at(nu, mu);

                    const term_ab = PS_ab.at(mu, nu) * PS_ab.at(nu, mu);

                    bo.ptr(A, B).* += 2 * (term_a + term_b + 2 * term_ab);
                }
            }
        }
    }

    return bo;
}

/// Computes Mulliken net atomic charges by partitioning the electronic density matrix using the overlap matrix.
pub fn mulliken(comptime T: type, sys: MolecularSystem(T), P: Matrix(T), S: Matrix(T), gpa: Allocator) !Vector(T) {
    var net_populations = try gpa.alloc(T, sys.atoms.len);
    defer gpa.free(net_populations);

    @memset(net_populations, 0);

    for (0..P.shape[0]) |mu| {
        net_populations[@intCast(sys.bf2at[mu % sys.nbf])] += dot(T, P.row(mu), S.row(mu));
    }

    var charges = try Vector(T).init(sys.atoms.len, gpa);
    errdefer charges.deinit(gpa);

    for (0..sys.atoms.len) |A| {
        charges.data[A] = @as(T, @floatFromInt(sys.atoms[A])) - net_populations[A];
    }

    return charges;
}

/// Formats and prints the calculated Löwdin atomic charges to the output.
pub fn printLowdinCharges(comptime T: type, io: std.Io, sys: MolecularSystem(T), charges: Vector(T), method_str: []const u8) !void {
    try printf(io, "\n{s} LÖWDIN POPULATION ANALYSIS\n", .{method_str});

    for (0..sys.atoms.len) |A| {
        const sym = try getSymbol(sys.atoms[A]);

        try printf(io, "{s:4} {d:20.14}\n", .{ sym, charges.data[A] });
    }
}

/// Formats and prints the calculated Mayer bond orders to the output.
pub fn printMayerBondOrders(comptime T: type, io: std.Io, sys: MolecularSystem(T), bo: Matrix(T), method_str: []const u8) !void {
    try printf(io, "\n{s} MAYER BOND ORDERS\n", .{method_str});

    for (0..sys.atoms.len) |A| for (A + 1..sys.atoms.len) |B| {
        const sym_A = try getSymbol(sys.atoms[A]);
        const sym_B = try getSymbol(sys.atoms[B]);

        try printf(io, "{s:2}{d:<2} - {s:2}{d:<2} {d:20.14}\n", .{ sym_A, A + 1, sym_B, B + 1, bo.at(A, B) });
    };
}

/// Formats and prints the calculated Mulliken atomic charges to the output.
pub fn printMullikenCharges(comptime T: type, io: std.Io, sys: MolecularSystem(T), charges: Vector(T), method_str: []const u8) !void {
    try printf(io, "\n{s} MULLIKEN POPULATION ANALYSIS\n", .{method_str});

    for (0..sys.atoms.len) |A| {
        const sym = try getSymbol(sys.atoms[A]);

        try printf(io, "{s:4} {d:20.14}\n", .{ sym, charges.data[A] });
    }
}

/// Formats and prints the calculated Wiberg bond indices to the output.
pub fn printWibergBondOrders(comptime T: type, io: std.Io, sys: MolecularSystem(T), bo: Matrix(T), method_str: []const u8) !void {
    try printf(io, "\n{s} WIBERG BOND ORDERS\n", .{method_str});

    for (0..sys.atoms.len) |A| for (A + 1..sys.atoms.len) |B| {
        const sym_A = try getSymbol(sys.atoms[A]);
        const sym_B = try getSymbol(sys.atoms[B]);

        try printf(io, "{s:2}{d:<2} - {s:2}{d:<2} {d:20.14}\n", .{ sym_A, A + 1, sym_B, B + 1, bo.at(A, B) });
    };
}

/// Computes Wiberg bond indices between all pairs of atoms in the symmetrically orthogonalized Löwdin basis.
pub fn wiberg(comptime T: type, sys: MolecularSystem(T), P: Matrix(T), S: Matrix(T), gpa: Allocator) !Matrix(T) {
    const nbf = sys.nbf;

    const is_gen = (P.shape[0] == 2 * nbf);

    var S_spatial = try Matrix(T).init(nbf, nbf, gpa);
    defer S_spatial.deinit(gpa);

    for (0..nbf) |mu| for (0..nbf) |nu| {
        S_spatial.ptr(mu, nu).* = S.at(mu, nu);
    };

    var Shalf = try matrixSquareRoot(T, S_spatial, gpa);
    defer Shalf.deinit(gpa);

    var bo = try Matrix(T).initZero(sys.atoms.len, sys.atoms.len, gpa);
    errdefer bo.deinit(gpa);

    if (!is_gen) {
        var SP = try Matrix(T).init(nbf, nbf, gpa);
        defer SP.deinit(gpa);

        gemm(T, &SP, Shalf, P, 1, 0, false, false);

        var P_ortho = try Matrix(T).init(nbf, nbf, gpa);
        defer P_ortho.deinit(gpa);

        gemm(T, &P_ortho, SP, Shalf, 1, 0, false, false);

        for (0..nbf) |mu| {
            const A: usize = @intCast(sys.bf2at[mu]);

            for (0..nbf) |nu| {
                const B: usize = @intCast(sys.bf2at[nu]);

                if (A != B) {
                    const p_munu = P_ortho.at(mu, nu);

                    bo.ptr(A, B).* += p_munu * p_munu;
                }
            }
        }
    }

    if (is_gen) {
        var blocks = try SpinBlocks(T).init(P, nbf, gpa);
        defer blocks.deinit(gpa);

        var SP_a = try Matrix(T).init(nbf, nbf, gpa);
        defer SP_a.deinit(gpa);

        var SP_b = try Matrix(T).init(nbf, nbf, gpa);
        defer SP_b.deinit(gpa);

        var SP_ab = try Matrix(T).init(nbf, nbf, gpa);
        defer SP_ab.deinit(gpa);

        gemm(T, &SP_a, Shalf, blocks.a, 1, 0, false, false);
        gemm(T, &SP_b, Shalf, blocks.b, 1, 0, false, false);

        gemm(T, &SP_ab, Shalf, blocks.ab, 1, 0, false, false);

        var P_ortho_a = try Matrix(T).init(nbf, nbf, gpa);
        defer P_ortho_a.deinit(gpa);

        var P_ortho_b = try Matrix(T).init(nbf, nbf, gpa);
        defer P_ortho_b.deinit(gpa);

        var P_ortho_ab = try Matrix(T).init(nbf, nbf, gpa);
        defer P_ortho_ab.deinit(gpa);

        gemm(T, &P_ortho_a, SP_a, Shalf, 1, 0, false, false);
        gemm(T, &P_ortho_b, SP_b, Shalf, 1, 0, false, false);

        gemm(T, &P_ortho_ab, SP_ab, Shalf, 1, 0, false, false);

        for (0..nbf) |mu| {
            const A: usize = @intCast(sys.bf2at[mu]);

            for (0..nbf) |nu| {
                const B: usize = @intCast(sys.bf2at[nu]);

                if (A != B) {
                    const pa_munu = P_ortho_a.at(mu, nu);
                    const pb_munu = P_ortho_b.at(mu, nu);
                    const pab_munu = P_ortho_ab.at(mu, nu);

                    bo.ptr(A, B).* += 2.0 * (pa_munu * pa_munu + pb_munu * pb_munu + 2.0 * pab_munu * pab_munu);
                }
            }
        }
    }

    return bo;
}

/// Computes the symmetric square root of a positive semi-definite matrix using spectral decomposition.
fn matrixSquareRoot(comptime T: type, S: Matrix(T), gpa: Allocator) !Matrix(T) {
    const n = S.nrow();

    var eigvals = try Vector(T).init(n, gpa);
    defer eigvals.deinit(gpa);

    var U = try Matrix(T).init(n, n, gpa);
    defer U.deinit(gpa);

    try syevd(T, &eigvals, &U, S);

    var U_scaled = try Matrix(T).init(n, n, gpa);
    defer U_scaled.deinit(gpa);

    for (0..n) |i| for (0..n) |j| {
        U_scaled.ptr(i, j).* = U.at(i, j) * @sqrt(@max(0, eigvals.at(j)));
    };

    var Shalf = try Matrix(T).init(n, n, gpa);
    errdefer Shalf.deinit(gpa);

    gemm(T, &Shalf, U_scaled, U, 1, 0, false, true);

    return Shalf;
}
