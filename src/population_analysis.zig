//! Implements Mulliken and Löwdin population analysis and Mayer and Wiberg bond order calculations.

const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Vector = @import("tensor.zig").Vector;

const dot = @import("linear_algebra.zig").dot;
const eigh = @import("linear_algebra.zig").eigh;
const getSymbol = @import("constant.zig").getSymbol;
const mm = @import("linear_algebra.zig").mm;
const printf = @import("read_write.zig").printf;

/// Computes Löwdin net atomic charges using symmetric orthogonalization of the overlap matrix.
pub fn lowdin(comptime T: type, sys: MolecularSystem(T), P: Matrix(T), S: Matrix(T), gpa: Allocator) !Vector(T) {
    const n = S.nrow();

    var eigvals = try Vector(T).init(n, gpa);
    defer eigvals.deinit(gpa);

    var U = try Matrix(T).init(n, n, gpa);
    defer U.deinit(gpa);

    try eigh(T, &eigvals, &U, S);

    var U_scaled = try Matrix(T).init(n, n, gpa);
    defer U_scaled.deinit(gpa);

    for (0..n) |i| for (0..n) |j| {
        U_scaled.ptr(i, j).* = U.at(i, j) * @sqrt(@max(0, eigvals.at(j)));
    };

    var Shalf = try Matrix(T).init(n, n, gpa);
    defer Shalf.deinit(gpa);

    mm(T, &Shalf, U_scaled, U, 1.0, 0.0, false, true);

    var SP = try Matrix(T).init(n, n, gpa);
    defer SP.deinit(gpa);

    mm(T, &SP, Shalf, P, 1.0, 0.0, false, false);

    var net_populations = try gpa.alloc(T, sys.atoms.len);
    defer gpa.free(net_populations);

    @memset(net_populations, 0);

    for (0..n) |u| {
        net_populations[@intCast(sys.bf2at[u % sys.nbf])] += dot(T, SP.row(u), Shalf.row(u));
    }

    var charges = try Vector(T).init(sys.atoms.len, gpa);
    errdefer charges.deinit(gpa);

    for (0..sys.atoms.len) |i| {
        charges.data[i] = @as(T, @floatFromInt(sys.atoms[i])) - net_populations[i];
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

    for (0..nbf) |i| for (0..nbf) |j| {
        S_spatial.ptr(i, j).* = S.at(i, j);
    };

    if (!is_gen) {
        var PS = try Matrix(T).init(nbf, nbf, gpa);
        defer PS.deinit(gpa);

        mm(T, &PS, P, S_spatial, 1, 0, false, false);

        for (0..nbf) |u| {
            const at_u: usize = @intCast(sys.bf2at[u]);

            for (0..nbf) |v| {
                const at_v: usize = @intCast(sys.bf2at[v]);

                if (at_u != at_v) {
                    bo.ptr(at_u, at_v).* += PS.at(u, v) * PS.at(v, u);
                }
            }
        }
    }

    if (is_gen) {
        var P_a = try Matrix(T).init(nbf, nbf, gpa);
        defer P_a.deinit(gpa);

        var P_b = try Matrix(T).init(nbf, nbf, gpa);
        defer P_b.deinit(gpa);

        var P_ab = try Matrix(T).init(nbf, nbf, gpa);
        defer P_ab.deinit(gpa);

        for (0..nbf) |i| for (0..nbf) |j| {
            P_a.ptr(i, j).* = P.at(i + 0 * nbf, j + 0 * nbf);
            P_b.ptr(i, j).* = P.at(i + 1 * nbf, j + 1 * nbf);

            P_ab.ptr(i, j).* = P.at(i + 0 * nbf, j + 1 * nbf);
        };

        var PS_a = try Matrix(T).init(nbf, nbf, gpa);
        defer PS_a.deinit(gpa);

        var PS_b = try Matrix(T).init(nbf, nbf, gpa);
        defer PS_b.deinit(gpa);

        var PS_ab = try Matrix(T).init(nbf, nbf, gpa);
        defer PS_ab.deinit(gpa);

        mm(T, &PS_a, P_a, S_spatial, 1, 0, false, false);
        mm(T, &PS_b, P_b, S_spatial, 1, 0, false, false);

        mm(T, &PS_ab, P_ab, S_spatial, 1, 0, false, false);

        for (0..nbf) |u| {
            const at_u: usize = @intCast(sys.bf2at[u]);

            for (0..nbf) |v| {
                const at_v: usize = @intCast(sys.bf2at[v]);

                if (at_u != at_v) {
                    const term_a = PS_a.at(u, v) * PS_a.at(v, u);
                    const term_b = PS_b.at(u, v) * PS_b.at(v, u);

                    const term_ab = PS_ab.at(u, v) * PS_ab.at(v, u);

                    bo.ptr(at_u, at_v).* += 2 * (term_a + term_b + 2 * term_ab);
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

    for (0..P.shape[0]) |u| {
        net_populations[@intCast(sys.bf2at[u % sys.nbf])] += dot(T, P.row(u), S.row(u));
    }

    var charges = try Vector(T).init(sys.atoms.len, gpa);
    errdefer charges.deinit(gpa);

    for (0..sys.atoms.len) |i| {
        charges.data[i] = @as(T, @floatFromInt(sys.atoms[i])) - net_populations[i];
    }

    return charges;
}

/// Formats and prints the calculated Löwdin atomic charges to the output.
pub fn printLowdinCharges(comptime T: type, io: std.Io, sys: MolecularSystem(T), charges: Vector(T), method_str: []const u8) !void {
    try printf(io, "\n{s} LÖWDIN POPULATION ANALYSIS\n", .{method_str});

    for (0..sys.atoms.len) |i| {
        const sym = try getSymbol(sys.atoms[i]);

        try printf(io, "{s:4} {d:20.14}\n", .{ sym, charges.data[i] });
    }
}

/// Formats and prints the calculated Mayer bond orders to the output.
pub fn printMayerBondOrders(comptime T: type, io: std.Io, sys: MolecularSystem(T), bo: Matrix(T), method_str: []const u8) !void {
    try printf(io, "\n{s} MAYER BOND ORDERS\n", .{method_str});

    for (0..sys.atoms.len) |i| for (i + 1..sys.atoms.len) |j| {
        const sym_i = try getSymbol(sys.atoms[i]);
        const sym_j = try getSymbol(sys.atoms[j]);

        try printf(io, "{s:2}{d:<2} - {s:2}{d:<2} {d:20.14}\n", .{ sym_i, i + 1, sym_j, j + 1, bo.at(i, j) });
    };
}

/// Formats and prints the calculated Mulliken atomic charges to the output.
pub fn printMullikenCharges(comptime T: type, io: std.Io, sys: MolecularSystem(T), charges: Vector(T), method_str: []const u8) !void {
    try printf(io, "\n{s} MULLIKEN POPULATION ANALYSIS\n", .{method_str});

    for (0..sys.atoms.len) |i| {
        const sym = try getSymbol(sys.atoms[i]);

        try printf(io, "{s:4} {d:20.14}\n", .{ sym, charges.data[i] });
    }
}

/// Formats and prints the calculated Wiberg bond indices to the output.
pub fn printWibergBondOrders(comptime T: type, io: std.Io, sys: MolecularSystem(T), bo: Matrix(T), method_str: []const u8) !void {
    try printf(io, "\n{s} WIBERG BOND ORDERS\n", .{method_str});

    for (0..sys.atoms.len) |i| for (i + 1..sys.atoms.len) |j| {
        const sym_i = try getSymbol(sys.atoms[i]);
        const sym_j = try getSymbol(sys.atoms[j]);

        try printf(io, "{s:2}{d:<2} - {s:2}{d:<2} {d:20.14}\n", .{ sym_i, i + 1, sym_j, j + 1, bo.at(i, j) });
    };
}

/// Computes Wiberg bond indices between all pairs of atoms in the symmetrically orthogonalized Löwdin basis.
pub fn wiberg(comptime T: type, sys: MolecularSystem(T), P: Matrix(T), S: Matrix(T), gpa: Allocator) !Matrix(T) {
    const nbf = sys.nbf;

    const is_gen = (P.shape[0] == 2 * nbf);

    var S_spatial = try Matrix(T).init(nbf, nbf, gpa);
    defer S_spatial.deinit(gpa);

    for (0..nbf) |i| for (0..nbf) |j| {
        S_spatial.ptr(i, j).* = S.at(i, j);
    };

    var eigvals = try Vector(T).init(nbf, gpa);
    defer eigvals.deinit(gpa);

    var U = try Matrix(T).init(nbf, nbf, gpa);
    defer U.deinit(gpa);

    try eigh(T, &eigvals, &U, S_spatial);

    var U_scaled = try Matrix(T).init(nbf, nbf, gpa);
    defer U_scaled.deinit(gpa);

    for (0..nbf) |i| for (0..nbf) |j| {
        U_scaled.ptr(i, j).* = U.at(i, j) * @sqrt(@max(0, eigvals.at(j)));
    };

    var Shalf = try Matrix(T).init(nbf, nbf, gpa);
    defer Shalf.deinit(gpa);

    mm(T, &Shalf, U_scaled, U, 1, 0, false, true);

    var bo = try Matrix(T).initZero(sys.atoms.len, sys.atoms.len, gpa);
    errdefer bo.deinit(gpa);

    if (!is_gen) {
        var SP = try Matrix(T).init(nbf, nbf, gpa);
        defer SP.deinit(gpa);

        mm(T, &SP, Shalf, P, 1, 0, false, false);

        var P_ortho = try Matrix(T).init(nbf, nbf, gpa);
        defer P_ortho.deinit(gpa);

        mm(T, &P_ortho, SP, Shalf, 1, 0, false, false);

        for (0..nbf) |u| {
            const at_u: usize = @intCast(sys.bf2at[u]);

            for (0..nbf) |v| {
                const at_v: usize = @intCast(sys.bf2at[v]);

                if (at_u != at_v) {
                    const p_uv = P_ortho.at(u, v);

                    bo.ptr(at_u, at_v).* += p_uv * p_uv;
                }
            }
        }
    }

    if (is_gen) {
        var P_a = try Matrix(T).init(nbf, nbf, gpa);
        defer P_a.deinit(gpa);

        var P_b = try Matrix(T).init(nbf, nbf, gpa);
        defer P_b.deinit(gpa);

        var P_ab = try Matrix(T).init(nbf, nbf, gpa);
        defer P_ab.deinit(gpa);

        for (0..nbf) |i| for (0..nbf) |j| {
            P_a.ptr(i, j).* = P.at(i + 0 * nbf, j + 0 * nbf);
            P_b.ptr(i, j).* = P.at(i + 1 * nbf, j + 1 * nbf);

            P_ab.ptr(i, j).* = P.at(i + 0 * nbf, j + 1 * nbf);
        };

        var SP_a = try Matrix(T).init(nbf, nbf, gpa);
        defer SP_a.deinit(gpa);

        var SP_b = try Matrix(T).init(nbf, nbf, gpa);
        defer SP_b.deinit(gpa);

        var SP_ab = try Matrix(T).init(nbf, nbf, gpa);
        defer SP_ab.deinit(gpa);

        mm(T, &SP_a, Shalf, P_a, 1, 0, false, false);
        mm(T, &SP_b, Shalf, P_b, 1, 0, false, false);

        mm(T, &SP_ab, Shalf, P_ab, 1, 0, false, false);

        var P_ortho_a = try Matrix(T).init(nbf, nbf, gpa);
        defer P_ortho_a.deinit(gpa);

        var P_ortho_b = try Matrix(T).init(nbf, nbf, gpa);
        defer P_ortho_b.deinit(gpa);

        var P_ortho_ab = try Matrix(T).init(nbf, nbf, gpa);
        defer P_ortho_ab.deinit(gpa);

        mm(T, &P_ortho_a, SP_a, Shalf, 1, 0, false, false);
        mm(T, &P_ortho_b, SP_b, Shalf, 1, 0, false, false);

        mm(T, &P_ortho_ab, SP_ab, Shalf, 1, 0, false, false);

        for (0..nbf) |u| {
            const at_u: usize = @intCast(sys.bf2at[u]);

            for (0..nbf) |v| {
                const at_v: usize = @intCast(sys.bf2at[v]);

                if (at_u != at_v) {
                    const pa_uv = P_ortho_a.at(u, v);
                    const pb_uv = P_ortho_b.at(u, v);
                    const pab_uv = P_ortho_ab.at(u, v);

                    bo.ptr(at_u, at_v).* += 2.0 * (pa_uv * pa_uv + pb_uv * pb_uv + 2.0 * pab_uv * pab_uv);
                }
            }
        }
    }

    return bo;
}
