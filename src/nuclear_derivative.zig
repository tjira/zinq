//! Numerical nuclear gradients and Hessians via finite differences for molecular geometry optimization.

const std = @import("std");

const libint = @import("cimport.zig").libint;

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Vector = @import("tensor.zig").Vector;

const ao2so_coef = @import("integral_transform.zig").ao2so_coef;
const determinant = @import("linear_algebra.zig").determinant;
const dot = @import("linear_algebra.zig").dot;
const exportIfBuiltin = @import("molecular_integrals.zig").exportIfBuiltin;
const generateDets = @import("configuration_interaction.zig").generateDets;
const getSymbol = @import("constant.zig").getSymbol;
const mm = @import("linear_algebra.zig").mm;
const printf = @import("read_write.zig").printf;

const A2BOHR = @import("constant.zig").A2BOHR;

/// Computes the nuclear gradient of the electronic energy using central finite differences.
pub fn calculateNumericalGradient(comptime T: type, io: std.Io, runFn: anytype, opt: anytype, sys: *MolecularSystem(T), log: bool, gpa: Allocator) anyerror!Matrix(T) {
    const generalized = if (@hasField(@TypeOf(opt), "hartree_fock")) opt.hartree_fock.generalized else opt.generalized;

    var grad = try Matrix(T).initZero(sys.atoms.len, 3, gpa);
    errdefer grad.deinit(gpa);

    const h = @as(T, @floatCast(opt.gradient.?.numeric.step));

    const state = if (@hasField(@TypeOf(opt.gradient.?.numeric), "state")) (opt.gradient.?.numeric.state orelse 0) else 0;

    if (log) try std.Io.File.stdout().writeStreamingAll(io, "\n");

    const nbf = if (generalized) 2 * sys.nbf else sys.nbf;

    var Pg = try Matrix(T).initZero(nbf, nbf, gpa);
    defer Pg.deinit(gpa);

    _ = try getE(T, io, runFn, opt, sys, .{}, null, &Pg, state, gpa);

    for (0..sys.atoms.len) |i| for (0..3) |j| {
        const d1, const d2 = .{ .{ i * 3 + j, h }, .{ i * 3 + j, -h } };

        const ep = try getE(T, io, runFn, opt, sys, d1, Pg, null, state, gpa);
        const em = try getE(T, io, runFn, opt, sys, d2, Pg, null, state, gpa);

        grad.ptr(i, j).* = (ep - em) / (2 * h);

        if (log) {
            const fmt = "DERIVATIVE OF ATOM {d:02} ({s:2}) IN DIRECTION {s}: {d:20.14}\n";

            const sym = try getSymbol(sys.atoms[i]);

            const dirs = &[_][]const u8{ "'x'", "'y'", "'z'" };

            try printf(io, fmt, .{ i, sym, dirs[j], grad.at(i, j) });
        }
    };

    return grad;
}

/// Computes the nuclear Hessian matrix of the electronic energy using finite differences of energies.
pub fn calculateNumericalHessian(comptime T: type, io: std.Io, runFn: anytype, opt: anytype, sys: *MolecularSystem(T), log: bool, gpa: Allocator) anyerror!Matrix(T) {
    const generalized = if (@hasField(@TypeOf(opt), "hartree_fock")) opt.hartree_fock.generalized else opt.generalized;

    var hess = try Matrix(T).initZero(3 * sys.atoms.len, 3 * sys.atoms.len, gpa);
    errdefer hess.deinit(gpa);

    const h = @as(T, @floatCast(opt.hessian.?.numeric.step));

    const state = if (@hasField(@TypeOf(opt.hessian.?.numeric), "state")) opt.hessian.?.numeric.state else 0;

    if (log) try std.Io.File.stdout().writeStreamingAll(io, "\n");

    const nbf = if (generalized) 2 * sys.nbf else sys.nbf;

    var Pg = try Matrix(T).initZero(nbf, nbf, gpa);
    defer Pg.deinit(gpa);

    const er = try getE(T, io, runFn, opt, sys, .{}, null, &Pg, state, gpa);

    for (0..3 * sys.atoms.len) |a| {
        const d1, const d2 = .{ .{ a, h }, .{ a, -h } };

        const ep = try getE(T, io, runFn, opt, sys, d1, Pg, null, state, gpa);
        const em = try getE(T, io, runFn, opt, sys, d2, Pg, null, state, gpa);

        hess.ptr(a, a).* = (ep - 2 * er + em) / (h * h);

        if (log) {
            const fmt = "SECOND DERIVATIVE OF ATOM {d:02} ({s:2}) DIR {s} AND ATOM {d:02} ({s:2}) DIR {s}: {d:20.14}\n";

            const sym = try getSymbol(sys.atoms[a / 3]);

            const dirs = &[_][]const u8{ "'x'", "'y'", "'z'" };

            try printf(io, fmt, .{ a / 3, sym, dirs[a % 3], a / 3, sym, dirs[a % 3], hess.at(a, a) });
        }

        for (a + 1..3 * sys.atoms.len) |b| {
            const epp = try getE(T, io, runFn, opt, sys, .{ a, h, b, h }, Pg, null, state, gpa);

            const epm = try getE(T, io, runFn, opt, sys, .{ a, h, b, -h }, Pg, null, state, gpa);
            const emp = try getE(T, io, runFn, opt, sys, .{ a, -h, b, h }, Pg, null, state, gpa);

            const emm = try getE(T, io, runFn, opt, sys, .{ a, -h, b, -h }, Pg, null, state, gpa);

            const val = (epp - epm - emp + emm) / (4 * h * h);

            hess.ptr(a, b).* = val;
            hess.ptr(b, a).* = val;

            if (log) {
                const fmt = "SECOND DERIVATIVE OF ATOM {d:02} ({s:2}) DIR {s} AND " ++
                    "ATOM {d:02} ({s:2}) DIR {s}: {d:20.14}\n";

                const sym_a = try getSymbol(sys.atoms[a / 3]);
                const sym_b = try getSymbol(sys.atoms[b / 3]);

                const dirs = &[_][]const u8{ "'x'", "'y'", "'z'" };

                try printf(io, fmt, .{ a / 3, sym_a, dirs[a % 3], b / 3, sym_b, dirs[b % 3], val });
            }
        }
    }

    return hess;
}

/// Computes nonadiabatic coupling vectors between electronic states using central finite differences.
pub fn calculateNumericalNacv(comptime T: type, io: std.Io, runFn: anytype, opt: anytype, sys: *MolecularSystem(T), log: bool, gpa: Allocator) anyerror!Matrix(T) {
    const generalized = if (@hasField(@TypeOf(opt), "hartree_fock")) opt.hartree_fock.generalized else opt.generalized;

    var nacv = try Matrix(T).initZero(sys.atoms.len, 3, gpa);
    errdefer nacv.deinit(gpa);

    const h = @as(T, @floatCast(opt.nacv.?.numeric.step));

    const s_j = @as(usize, @intCast(opt.nacv.?.numeric.states[0]));
    const s_k = @as(usize, @intCast(opt.nacv.?.numeric.states[1]));

    if (log) try std.Io.File.stdout().writeStreamingAll(io, "\n");

    const nbf = if (generalized) 2 * sys.nbf else sys.nbf;

    var Pg = try Matrix(T).initZero(nbf, nbf, gpa);
    defer Pg.deinit(gpa);

    var modified_opt = opt;

    modified_opt.gradient, modified_opt.hessian, modified_opt.optimize = .{ null, null, null };

    if (comptime @hasField(@TypeOf(modified_opt), "nacv")) {
        modified_opt.nacv = null;
    }

    if (comptime @hasField(@TypeOf(modified_opt), "nstate")) {
        const required = @as(u32, @intCast(@max(s_j, s_k) + 1));

        if (modified_opt.nstate < required) {
            std.log.err("REQUESTED NACV STATE COUNT IS NOT AVAILABLE IN THE ELECTRONIC STRUCTURE CALCULATION", .{});

            return error.InvalidInput;
        }
    }

    if (comptime @hasField(@TypeOf(modified_opt), "frequency")) {
        modified_opt.frequency = null;
    }

    if (comptime @hasField(@TypeOf(modified_opt), "lowdin")) {
        modified_opt.lowdin = false;
    }

    if (comptime @hasField(@TypeOf(modified_opt), "mayer")) {
        modified_opt.mayer = false;
    }

    if (comptime @hasField(@TypeOf(modified_opt), "mulliken")) {
        modified_opt.mulliken = false;
    }

    if (comptime @hasField(@TypeOf(modified_opt), "wiberg")) {
        modified_opt.wiberg = false;
    }

    if (comptime @hasField(@TypeOf(modified_opt), "response")) {
        modified_opt.response = null;
    }

    if (comptime @hasField(@TypeOf(modified_opt), "write")) {
        modified_opt.write = .{};
    }

    if (comptime @hasField(@TypeOf(modified_opt), "hartree_fock")) {
        modified_opt.hartree_fock.frequency = null;

        modified_opt.hartree_fock.gradient = null;
        modified_opt.hartree_fock.response = null;

        modified_opt.hartree_fock.hessian = null;
        modified_opt.hartree_fock.lowdin = false;
        modified_opt.hartree_fock.wiberg = false;

        modified_opt.hartree_fock.mayer, modified_opt.hartree_fock.mulliken = .{ false, false };

        modified_opt.hartree_fock.write = .{};
    }

    var res_ref = try runFn(T, io, modified_opt, sys, null, false, gpa);
    defer res_ref.deinit(gpa);

    @memcpy(Pg.data, (if (@hasField(@TypeOf(res_ref), "hartree_fock")) res_ref.hartree_fock.P else res_ref.P).data);

    if (res_ref.C.ncol() <= s_j or res_ref.C.ncol() <= s_k) {
        std.log.err("STATE REQUESTED FOR NUMERICAL NACV IS NOT AVAILABLE IN THE ELECTRONIC STRUCTURE CALCULATION", .{});

        return error.InvalidInput;
    }

    var c_j_ref = try Vector(T).init(res_ref.C.nrow(), gpa);
    defer c_j_ref.deinit(gpa);

    var c_k_ref = try Vector(T).init(res_ref.C.nrow(), gpa);
    defer c_k_ref.deinit(gpa);

    for (0..res_ref.C.nrow()) |r| {
        c_j_ref.data[r] = res_ref.C.at(r, s_j);
        c_k_ref.data[r] = res_ref.C.at(r, s_k);
    }

    var dets = try generateDets(sys.nel, 2 * sys.nbf, opt.excitations, gpa);

    defer {
        for (0..dets.items.len) |i| gpa.free(dets.items[i]);

        dets.deinit(gpa);
    }

    var M = try Matrix(T).init(sys.nel, sys.nel, gpa);
    defer M.deinit(gpa);

    const ipiv = try gpa.alloc(i32, sys.nel);
    defer gpa.free(ipiv);

    var sys_m = try sys.clone(gpa);
    defer sys_m.deinit(gpa);

    var sys_p = try sys.clone(gpa);
    defer sys_p.deinit(gpa);

    const nthreads = if (@hasField(@TypeOf(opt), "nthreads")) opt.nthreads else 1;

    const nsp = 2 * sys.nbf;

    for (0..sys.atoms.len) |i| for (0..3) |j| {
        const coord_idx = i * 3 + j;

        sys_m.coors[coord_idx] -= h;

        libint.libint_update_coords(sys_m.ptr, sys_m.coors.ptr);

        defer {
            sys_m.coors[coord_idx] += h;

            libint.libint_update_coords(sys_m.ptr, sys_m.coors.ptr);
        }

        sys_p.coors[coord_idx] += h;

        libint.libint_update_coords(sys_p.ptr, sys_p.coors.ptr);

        defer {
            sys_p.coors[coord_idx] -= h;

            libint.libint_update_coords(sys_p.ptr, sys_p.coors.ptr);
        }

        var res_m = try runFn(T, io, modified_opt, &sys_m, Pg, false, gpa);
        defer res_m.deinit(gpa);

        var res_p = try runFn(T, io, modified_opt, &sys_p, Pg, false, gpa);
        defer res_p.deinit(gpa);

        var c_j_m = try Vector(T).init(res_m.C.nrow(), gpa);
        defer c_j_m.deinit(gpa);

        var c_k_m = try Vector(T).init(res_m.C.nrow(), gpa);
        defer c_k_m.deinit(gpa);

        for (0..res_m.C.nrow()) |r| {
            c_j_m.data[r] = res_m.C.at(r, s_j);
            c_k_m.data[r] = res_m.C.at(r, s_k);
        }

        var c_j_p = try Vector(T).init(res_p.C.nrow(), gpa);
        defer c_j_p.deinit(gpa);

        var c_k_p = try Vector(T).init(res_p.C.nrow(), gpa);
        defer c_k_p.deinit(gpa);

        for (0..res_p.C.nrow()) |r| {
            c_j_p.data[r] = res_p.C.at(r, s_j);
            c_k_p.data[r] = res_p.C.at(r, s_k);
        }

        if (dot(T, c_j_m, c_j_ref) < 0) c_j_m.muls(-1);
        if (dot(T, c_k_m, c_k_ref) < 0) c_k_m.muls(-1);
        if (dot(T, c_j_p, c_j_ref) < 0) c_j_p.muls(-1);
        if (dot(T, c_k_p, c_k_ref) < 0) c_k_p.muls(-1);

        var S_SO = try sys_m.overlapBetweenSpin(sys_p, nthreads, gpa);
        defer S_SO.deinit(gpa);

        const C_m = if (@hasField(@TypeOf(res_m), "hartree_fock")) res_m.hartree_fock.C else res_m.C;
        const C_p = if (@hasField(@TypeOf(res_p), "hartree_fock")) res_p.hartree_fock.C else res_p.C;

        var C_m_SO = if (generalized) C_m else try Matrix(T).initZero(nsp, nsp, gpa);
        defer if (!generalized) C_m_SO.deinit(gpa);

        if (!generalized) ao2so_coef(T, &C_m_SO, C_m);

        var C_p_SO = if (generalized) C_p else try Matrix(T).initZero(nsp, nsp, gpa);
        defer if (!generalized) C_p_SO.deinit(gpa);

        if (!generalized) ao2so_coef(T, &C_p_SO, C_p);

        var tmp = try Matrix(T).init(nsp, nsp, gpa);
        defer tmp.deinit(gpa);

        mm(T, &tmp, S_SO, C_p_SO, 1, 0, false, false);

        var S_MO = try Matrix(T).init(nsp, nsp, gpa);
        defer S_MO.deinit(gpa);

        mm(T, &S_MO, C_m_SO, tmp, 1, 0, true, false);

        var s_jk: T = 0;
        var s_kj: T = 0;

        for (0..dets.items.len) |a| for (0..dets.items.len) |b| {
            const prod_jk = c_j_m.data[a] * c_k_p.data[b];
            const prod_kj = c_k_m.data[a] * c_j_p.data[b];

            if (@abs(prod_jk) < 1e-15 and @abs(prod_kj) < 1e-15) continue;

            for (0..sys.nel) |r| for (0..sys.nel) |c| {
                M.ptr(r, c).* = S_MO.at(dets.items[a][r], dets.items[b][c]);
            };

            const det = try determinant(T, &M, ipiv);

            s_jk += prod_jk * det;
            s_kj += prod_kj * det;
        };

        nacv.ptr(i, j).* = (s_jk - s_kj) / (4 * h);

        if (log) {
            const fmt = "NACV OF ATOM {d:02} ({s:2}) IN DIRECTION {s}: {d:20.14}\n";

            const dirs = &[_][]const u8{ "'x'", "'y'", "'z'" };

            try printf(io, fmt, .{ i, try getSymbol(sys.atoms[i]), dirs[j], nacv.at(i, j) });
        }
    };

    return nacv;
}

/// Evaluates the electronic energy for a perturbed molecular geometry.
fn getE(comptime T: type, io: std.Io, runFn: anytype, opt: anytype, sys: *MolecularSystem(T), perts: anytype, Pg: ?Matrix(T), P: ?*Matrix(T), state: usize, gpa: Allocator) !T {
    std.debug.assert(@typeInfo(@TypeOf(perts)).@"struct".field_names.len % 2 == 0);

    var modified_opt = opt;

    modified_opt.gradient, modified_opt.hessian, modified_opt.optimize = .{ null, null, null };

    if (comptime @hasField(@TypeOf(modified_opt), "nacv")) {
        modified_opt.nacv = null;
    }

    if (comptime @hasField(@TypeOf(modified_opt), "frequency")) {
        modified_opt.frequency = null;
    }

    if (comptime @hasField(@TypeOf(modified_opt), "lowdin")) {
        modified_opt.lowdin = false;
    }

    if (comptime @hasField(@TypeOf(modified_opt), "mayer")) {
        modified_opt.mayer = false;
    }

    if (comptime @hasField(@TypeOf(modified_opt), "mulliken")) {
        modified_opt.mulliken = false;
    }

    if (comptime @hasField(@TypeOf(modified_opt), "wiberg")) {
        modified_opt.wiberg = false;
    }

    if (comptime @hasField(@TypeOf(modified_opt), "response")) {
        modified_opt.response = null;
    }

    if (comptime @hasField(@TypeOf(modified_opt), "write")) {
        modified_opt.write = .{};
    }

    if (comptime @hasField(@TypeOf(modified_opt), "hartree_fock")) {
        modified_opt.hartree_fock.frequency = null;

        modified_opt.hartree_fock.gradient = null;
        modified_opt.hartree_fock.response = null;

        modified_opt.hartree_fock.hessian = null;
        modified_opt.hartree_fock.lowdin = false;
        modified_opt.hartree_fock.wiberg = false;

        modified_opt.hartree_fock.mayer, modified_opt.hartree_fock.mulliken = .{ false, false };

        modified_opt.hartree_fock.write = .{};
    }

    var original_coors: [@typeInfo(@TypeOf(perts)).@"struct".field_names.len / 2]T = undefined;

    inline for (0..@typeInfo(@TypeOf(perts)).@"struct".field_names.len / 2) |idx| {
        original_coors[idx] = sys.coors[perts[idx * 2]];
        sys.coors[perts[idx * 2]] += perts[idx * 2 + 1];
    }

    libint.libint_update_coords(sys.ptr, sys.coors.ptr);

    defer {
        inline for (0..@typeInfo(@TypeOf(perts)).@"struct".field_names.len / 2) |idx| {
            sys.coors[perts[idx * 2 + 0]] = original_coors[idx];
        }

        libint.libint_update_coords(sys.ptr, sys.coors.ptr);
    }

    var res = try runFn(T, io, modified_opt, sys, Pg, false, gpa);
    defer res.deinit(gpa);

    if (P) |out| {
        @memcpy(out.data, (if (@hasField(@TypeOf(res), "hartree_fock")) res.hartree_fock.P else res.P).data);
    }

    if (res.energy.len <= state) {
        const msg = "STATE REQUESTED FOR NUMERICAL DERIVATIVE IS NOT AVAILABLE IN THE ELECTRONIC STRUCTURE CALCULATION";

        std.log.err(msg, .{});

        return error.InvalidInput;
    }

    return res.energy[state];
}
