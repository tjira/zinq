//! Numerical nuclear gradients and Hessians via finite differences for molecular geometry optimization.

const std = @import("std");

const libint = @import("cimport.zig").libint;

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Vector = @import("tensor.zig").Vector;

const ao2so_coef = @import("integral_transform.zig").ao2so_coef;
const det = @import("linear_algebra.zig").det;
const dot = @import("linear_algebra.zig").dot;
const gemm = @import("linear_algebra.zig").gemm;
const generateCasDets = @import("configuration_interaction.zig").generateCasDets;
const generateDets = @import("configuration_interaction.zig").generateDets;
const getSymbol = @import("constant.zig").getSymbol;
const printf = @import("read_write.zig").printf;

/// Execution context encapsulating molecular system geometry, IO streams, memory allocator, and logging settings.
pub fn DerivativeContext(comptime T: type) type {
    return struct { gpa: Allocator, io: std.Io, log: bool = false, sys: *MolecularSystem(T), Pg: ?Matrix(T) = null };
}

/// Computes nuclear gradients of electronic states using central finite differences.
pub fn calculateNumericalGradient(comptime T: type, ctx: DerivativeContext(T), runFn: anytype, opt: anytype, states: []const u32) anyerror![]Matrix(T) {
    const generalized = if (@hasField(@TypeOf(opt), "hartree_fock")) opt.hartree_fock.generalized else opt.generalized;

    const io, const sys, const gpa, const log = .{ ctx.io, ctx.sys, ctx.gpa, ctx.log };

    var grad = try gpa.alloc(Matrix(T), states.len);
    errdefer gpa.free(grad);

    var grad_done: usize = 0;

    errdefer for (0..grad_done) |k| grad[k].deinit(gpa);

    for (0..states.len) |k| {
        grad[k] = try Matrix(T).initZero(sys.atoms.len, 3, gpa);

        grad_done += 1;
    }

    const h = @as(T, @floatCast(opt.gradient.?.numeric.step));

    if (log) try printf(io, "\n", .{});

    const nbf = if (generalized) 2 * sys.nbf else sys.nbf;

    var Pg = try Matrix(T).initZero(nbf, nbf, gpa);
    defer Pg.deinit(gpa);

    var max_state: usize = 0;

    for (states) |s| {
        if (s > max_state) max_state = s;
    }

    const ep = try gpa.alloc(T, max_state + 1);
    defer gpa.free(ep);

    const em = try gpa.alloc(T, max_state + 1);
    defer gpa.free(em);

    try getEnergies(T, ctx, runFn, opt, .{}, ctx.Pg, &Pg, ep);

    for (0..sys.atoms.len) |i| for (0..3) |j| {
        const d1, const d2 = .{ .{ i * 3 + j, h }, .{ i * 3 + j, -h } };

        try getEnergies(T, ctx, runFn, opt, d1, Pg, null, ep);
        try getEnergies(T, ctx, runFn, opt, d2, Pg, null, em);

        for (states, 0..) |s, s_idx| {
            grad[s_idx].ptr(i, j).* = (ep[s] - em[s]) / (2 * h);

            if (log) {
                const sym = try getSymbol(sys.atoms[i]);

                const dirs = &[_][]const u8{ "'x'", "'y'", "'z'" };

                if (comptime @hasField(@TypeOf(opt), "nstate")) {
                    const fmt = "NUCLEAR GRADIENT STATE {d} FOR ATOM {d:02} ({s:2}) IN DIRECTION {s}: {d:20.14}\n";

                    try printf(io, fmt, .{ s, i, sym, dirs[j], grad[s_idx].at(i, j) });
                }

                if (!comptime @hasField(@TypeOf(opt), "nstate")) {
                    const fmt = "NUCLEAR GRADIENT FOR ATOM {d:02} ({s:2}) IN DIRECTION {s}: {d:20.14}\n";

                    try printf(io, fmt, .{ i, sym, dirs[j], grad[s_idx].at(i, j) });
                }
            }
        }
    };

    return grad;
}

/// Computes nuclear Hessians of electronic states using finite differences of energies.
pub fn calculateNumericalHessian(comptime T: type, ctx: DerivativeContext(T), runFn: anytype, opt: anytype, states: []const u32) anyerror![]Matrix(T) {
    const generalized = if (@hasField(@TypeOf(opt), "hartree_fock")) opt.hartree_fock.generalized else opt.generalized;

    const io, const sys, const gpa, const log = .{ ctx.io, ctx.sys, ctx.gpa, ctx.log };

    var hess = try gpa.alloc(Matrix(T), states.len);
    errdefer gpa.free(hess);

    var hess_done: usize = 0;

    errdefer for (0..hess_done) |k| hess[k].deinit(gpa);

    for (0..states.len) |k| {
        hess[k] = try Matrix(T).initZero(3 * sys.atoms.len, 3 * sys.atoms.len, gpa);

        hess_done += 1;
    }

    const h = @as(T, @floatCast(opt.hessian.?.numeric.step));

    if (log) try printf(io, "\n", .{});

    const nbf = if (generalized) 2 * sys.nbf else sys.nbf;

    var Pg = try Matrix(T).initZero(nbf, nbf, gpa);
    defer Pg.deinit(gpa);

    var max_state: usize = 0;

    for (states) |s| {
        if (s > max_state) max_state = s;
    }

    const er = try gpa.alloc(T, max_state + 1);
    defer gpa.free(er);

    const ep = try gpa.alloc(T, max_state + 1);
    defer gpa.free(ep);

    const em = try gpa.alloc(T, max_state + 1);
    defer gpa.free(em);

    const epp = try gpa.alloc(T, max_state + 1);
    defer gpa.free(epp);

    const epm = try gpa.alloc(T, max_state + 1);
    defer gpa.free(epm);

    const emp = try gpa.alloc(T, max_state + 1);
    defer gpa.free(emp);

    const emm = try gpa.alloc(T, max_state + 1);
    defer gpa.free(emm);

    try getEnergies(T, ctx, runFn, opt, .{}, ctx.Pg, &Pg, er);

    for (0..3 * sys.atoms.len) |a| {
        const d1, const d2 = .{ .{ a, h }, .{ a, -h } };

        try getEnergies(T, ctx, runFn, opt, d1, Pg, null, ep);
        try getEnergies(T, ctx, runFn, opt, d2, Pg, null, em);

        for (states, 0..) |s, s_idx| {
            hess[s_idx].ptr(a, a).* = (ep[s] - 2 * er[s] + em[s]) / (h * h);

            if (log) {
                const at = a / 3;

                const sym, const dirs = .{ try getSymbol(sys.atoms[at]), &[_][]const u8{ "'x'", "'y'", "'z'" } };

                if (comptime @hasField(@TypeOf(opt), "nstate")) {
                    const fmt = "NUCLEAR HESSIAN STATE {d} FOR ATOM {d:02} ({s:2}) DIR {s} AND " ++
                        "ATOM {d:02} ({s:2}) DIR {s}: {d:20.14}\n";

                    try printf(io, fmt, .{ s, at, sym, dirs[a % 3], at, sym, dirs[a % 3], hess[s_idx].at(a, a) });
                }

                if (!comptime @hasField(@TypeOf(opt), "nstate")) {
                    const fmt = "NUCLEAR HESSIAN FOR ATOM {d:02} ({s:2}) DIR {s} AND " ++
                        "ATOM {d:02} ({s:2}) DIR {s}: {d:20.14}\n";

                    try printf(io, fmt, .{ at, sym, dirs[a % 3], at, sym, dirs[a % 3], hess[s_idx].at(a, a) });
                }
            }
        }

        for (a + 1..3 * sys.atoms.len) |b| {
            try getEnergies(T, ctx, runFn, opt, .{ a, h, b, h }, Pg, null, epp);
            try getEnergies(T, ctx, runFn, opt, .{ a, h, b, -h }, Pg, null, epm);
            try getEnergies(T, ctx, runFn, opt, .{ a, -h, b, h }, Pg, null, emp);
            try getEnergies(T, ctx, runFn, opt, .{ a, -h, b, -h }, Pg, null, emm);

            for (states, 0..) |s, s_idx| {
                const val = (epp[s] - epm[s] - emp[s] + emm[s]) / (4 * h * h);

                hess[s_idx].ptr(a, b).* = val;
                hess[s_idx].ptr(b, a).* = val;

                if (log) {
                    const sym_a = try getSymbol(sys.atoms[a / 3]);
                    const sym_b = try getSymbol(sys.atoms[b / 3]);

                    const dirs = &[_][]const u8{ "'x'", "'y'", "'z'" };

                    if (comptime @hasField(@TypeOf(opt), "nstate")) {
                        const fmt = "NUCLEAR HESSIAN STATE {d} FOR ATOM {d:02} ({s:2}) DIR {s} AND " ++
                            "ATOM {d:02} ({s:2}) DIR {s}: {d:20.14}\n";

                        try printf(io, fmt, .{ s, a / 3, sym_a, dirs[a % 3], b / 3, sym_b, dirs[b % 3], val });
                    }

                    if (!comptime @hasField(@TypeOf(opt), "nstate")) {
                        const fmt = "NUCLEAR HESSIAN FOR ATOM {d:02} ({s:2}) DIR {s} AND " ++
                            "ATOM {d:02} ({s:2}) DIR {s}: {d:20.14}\n";

                        try printf(io, fmt, .{ a / 3, sym_a, dirs[a % 3], b / 3, sym_b, dirs[b % 3], val });
                    }
                }
            }
        }
    }

    return hess;
}

/// Computes nonadiabatic coupling vectors between state pairs using central finite differences.
pub fn calculateNumericalNacv(comptime T: type, ctx: DerivativeContext(T), runFn: anytype, opt: anytype, pairs: []const [2]u32) anyerror![]Matrix(T) {
    const generalized = if (@hasField(@TypeOf(opt), "hartree_fock")) opt.hartree_fock.generalized else opt.generalized;

    const io, const sys, const gpa, const log = .{ ctx.io, ctx.sys, ctx.gpa, ctx.log };

    var nacvs = try gpa.alloc(Matrix(T), pairs.len);
    errdefer gpa.free(nacvs);

    var nacvs_done: usize = 0;

    errdefer for (0..nacvs_done) |k| nacvs[k].deinit(gpa);

    for (0..pairs.len) |k| {
        nacvs[k] = try Matrix(T).initZero(sys.atoms.len, 3, gpa);

        nacvs_done += 1;
    }

    const h = @as(T, @floatCast(opt.nacv.?.numeric.step));

    var max_state: usize = 0;

    for (pairs) |pair| {
        if (pair[0] > max_state) max_state = pair[0];
        if (pair[1] > max_state) max_state = pair[1];
    }

    if (log) try printf(io, "\n", .{});

    const nbf = if (generalized) 2 * sys.nbf else sys.nbf;

    var Pg = try Matrix(T).initZero(nbf, nbf, gpa);
    defer Pg.deinit(gpa);

    const modified_opt = stripSecondaryOptions(opt, .{});

    if (comptime @hasField(@TypeOf(modified_opt), "nstate")) {
        const required = @as(u32, @intCast(max_state + 1));

        if (modified_opt.nstate < required) {
            std.log.err("REQUESTED NACV STATE COUNT IS NOT AVAILABLE IN THE ELECTRONIC STRUCTURE CALCULATION", .{});

            return error.InvalidInput;
        }
    }

    var res_ref = try runFn(T, io, modified_opt, sys, ctx.Pg, false, gpa);
    defer res_ref.deinit(gpa);

    @memcpy(Pg.data, (if (@hasField(@TypeOf(res_ref), "hartree_fock")) res_ref.hartree_fock.P else res_ref.P).data);

    const C_ci_ref = res_ref.C_ci;

    if (C_ci_ref.ncol() <= max_state) {
        std.log.err("STATE REQUESTED FOR NUMERICAL NACV IS NOT AVAILABLE IN THE ELECTRONIC STRUCTURE CALCULATION", .{});

        return error.InvalidInput;
    }

    var c_ref = try gpa.alloc(Vector(T), max_state + 1);

    var c_ref_done: usize = 0;

    defer {
        for (0..c_ref_done) |k| c_ref[k].deinit(gpa);

        gpa.free(c_ref);
    }

    for (0..max_state + 1) |k| {
        c_ref[k] = try Vector(T).init(C_ci_ref.nrow(), gpa);

        c_ref_done += 1;

        for (0..C_ci_ref.nrow()) |r| {
            c_ref[k].data[r] = C_ci_ref.at(r, k);
        }
    }

    var dets: std.ArrayList([]const usize) = undefined;

    if (comptime @hasField(@TypeOf(opt), "active_electrons")) {
        const ncore, const mult = .{ (sys.nel - opt.active_electrons) / 2, opt.hartree_fock.multiplicity };

        dets = try generateCasDets(ncore, opt.active_orbitals, opt.active_electrons, mult, gpa);
    }

    if (comptime @hasField(@TypeOf(opt), "space")) {
        dets = switch (opt.space) {
            .cas => |c| blk: {
                const ncore, const mult = .{ (sys.nel - c.active_electrons) / 2, opt.hartree_fock.multiplicity };

                break :blk try generateCasDets(ncore, c.active_orbitals, c.active_electrons, mult, gpa);
            },
            .excitations => |ex| try generateDets(sys.nel, 2 * sys.nbf, ex, gpa),
        };
    }

    if (!comptime @hasField(@TypeOf(opt), "active_electrons") and !@hasField(@TypeOf(opt), "space")) {
        std.log.err("NACV CALCULATION REQUIRES EITHER ACTIVE SPACE OR EXCITATION SPACE SPECIFICATION", .{});

        return error.InvalidInput;
    }

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

    const C_ref = if (comptime @hasField(@TypeOf(res_ref), "C")) res_ref.C else res_ref.hartree_fock.C;

    var C_ref_SO = if (generalized) C_ref else try Matrix(T).initZero(nsp, nsp, gpa);
    defer if (!generalized) C_ref_SO.deinit(gpa);

    if (!generalized) ao2so_coef(T, &C_ref_SO, C_ref);

    var c_m = try gpa.alloc(Vector(T), max_state + 1);

    var c_m_done: usize = 0;

    defer {
        for (0..c_m_done) |k| c_m[k].deinit(gpa);

        gpa.free(c_m);
    }

    var c_p = try gpa.alloc(Vector(T), max_state + 1);

    var c_p_done: usize = 0;

    defer {
        for (0..c_p_done) |k| c_p[k].deinit(gpa);

        gpa.free(c_p);
    }

    for (0..max_state + 1) |k| {
        c_m[k] = try Vector(T).init(C_ci_ref.nrow(), gpa);

        c_m_done += 1;

        c_p[k] = try Vector(T).init(C_ci_ref.nrow(), gpa);

        c_p_done += 1;
    }

    var s_jk = try gpa.alloc(T, pairs.len);
    defer gpa.free(s_jk);

    var s_kj = try gpa.alloc(T, pairs.len);
    defer gpa.free(s_kj);

    var prods_jk = try gpa.alloc(T, pairs.len);
    defer gpa.free(prods_jk);

    var prods_kj = try gpa.alloc(T, pairs.len);
    defer gpa.free(prods_kj);

    var tmp = try Matrix(T).init(nsp, nsp, gpa);
    defer tmp.deinit(gpa);

    var S_MO = try Matrix(T).init(nsp, nsp, gpa);
    defer S_MO.deinit(gpa);

    var C_m_buf = if (!generalized) try Matrix(T).initZero(nsp, nsp, gpa) else undefined;
    defer if (!generalized) C_m_buf.deinit(gpa);

    var C_p_buf = if (!generalized) try Matrix(T).initZero(nsp, nsp, gpa) else undefined;
    defer if (!generalized) C_p_buf.deinit(gpa);

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

        for (0..max_state + 1) |k| {
            const C_ci_m = res_m.C_ci;
            const C_ci_p = res_p.C_ci;

            for (0..C_ci_m.nrow()) |r| {
                c_m[k].data[r] = C_ci_m.at(r, k);
            }

            for (0..C_ci_p.nrow()) |r| {
                c_p[k].data[r] = C_ci_p.at(r, k);
            }
        }

        var S_SO = try sys_m.overlapBetweenSpin(sys_p, nthreads, gpa);
        defer S_SO.deinit(gpa);

        const C_m = if (comptime @hasField(@TypeOf(res_m), "C")) res_m.C else res_m.hartree_fock.C;
        const C_p = if (comptime @hasField(@TypeOf(res_p), "C")) res_p.C else res_p.hartree_fock.C;

        var C_m_SO = if (generalized) C_m else C_m_buf;

        if (!generalized) ao2so_coef(T, &C_m_SO, C_m);

        var C_p_SO = if (generalized) C_p else C_p_buf;

        if (!generalized) ao2so_coef(T, &C_p_SO, C_p);

        for (0..nsp) |p| {
            var d_m: T = 0;
            var d_p: T = 0;

            for (0..nsp) |mu| {
                d_m += C_m_SO.at(mu, p) * C_ref_SO.at(mu, p);
                d_p += C_p_SO.at(mu, p) * C_ref_SO.at(mu, p);
            }

            if (d_m < 0) {
                for (0..nsp) |mu| {
                    C_m_SO.ptr(mu, p).* = -C_m_SO.at(mu, p);
                }

                for (0..dets.items.len) |a| {
                    if (std.mem.indexOfScalar(usize, dets.items[a], p) != null) {
                        for (0..max_state + 1) |k| {
                            c_m[k].data[a] = -c_m[k].data[a];
                        }
                    }
                }
            }

            if (d_p < 0) {
                for (0..nsp) |mu| {
                    C_p_SO.ptr(mu, p).* = -C_p_SO.at(mu, p);
                }

                for (0..dets.items.len) |a| {
                    if (std.mem.indexOfScalar(usize, dets.items[a], p) != null) {
                        for (0..max_state + 1) |k| {
                            c_p[k].data[a] = -c_p[k].data[a];
                        }
                    }
                }
            }
        }

        for (0..max_state + 1) |k| {
            if (dot(T, c_m[k], c_ref[k]) < 0) c_m[k].muls(-1);
            if (dot(T, c_p[k], c_ref[k]) < 0) c_p[k].muls(-1);
        }

        gemm(T, &tmp, S_SO, C_p_SO, 1, 0, false, false);

        gemm(T, &S_MO, C_m_SO, tmp, 1, 0, true, false);

        @memset(s_jk, 0);
        @memset(s_kj, 0);

        for (0..dets.items.len) |a| for (0..dets.items.len) |b| {
            var any_active = false;

            for (pairs, 0..) |pair, p_idx| {
                prods_jk[p_idx] = c_m[pair[0]].data[a] * c_p[pair[1]].data[b];
                prods_kj[p_idx] = c_m[pair[1]].data[a] * c_p[pair[0]].data[b];

                if (@abs(prods_jk[p_idx]) >= 1e-15 or @abs(prods_kj[p_idx]) >= 1e-15) {
                    any_active = true;
                }
            }

            if (!any_active) continue;

            for (0..sys.nel) |r| for (0..sys.nel) |c| {
                M.ptr(r, c).* = S_MO.at(dets.items[a][r], dets.items[b][c]);
            };

            const d = try det(T, &M, ipiv);

            for (0..pairs.len) |p_idx| {
                s_jk[p_idx] += prods_jk[p_idx] * d;
                s_kj[p_idx] += prods_kj[p_idx] * d;
            }
        };

        for (pairs, 0..) |pair, p_idx| {
            nacvs[p_idx].ptr(i, j).* = (s_jk[p_idx] - s_kj[p_idx]) / (4 * h);

            if (log) {
                const dirs = &[_][]const u8{ "'x'", "'y'", "'z'" };

                const sym = try getSymbol(sys.atoms[i]);

                const fmt = "NONADIABATIC COUPLING STATES {d}/{d} FOR ATOM {d:02} ({s:2}) " ++
                    "IN DIRECTION {s}: {d:20.14}\n";

                try printf(io, fmt, .{ pair[0], pair[1], i, sym, dirs[j], nacvs[p_idx].at(i, j) });
            }
        }
    };

    return nacvs;
}

/// Strips secondary calculation blocks and property requests from electronic structure options.
pub fn stripSecondaryOptions(opt: anytype, preserve: anytype) @TypeOf(opt) {
    if (@typeInfo(@TypeOf(opt)) == .@"union") {
        return switch (opt) {
            inline else => |m, tag| @unionInit(@TypeOf(opt), @tagName(tag), stripSecondaryOptions(m, preserve)),
        };
    }

    var clean = opt;

    const preserve_grad = if (comptime @hasField(@TypeOf(preserve), "gradient")) preserve.gradient else false;

    const preserve_hess = if (comptime @hasField(@TypeOf(preserve), "hessian")) preserve.hessian else false;

    const preserve_nacv = if (comptime @hasField(@TypeOf(preserve), "nacv")) preserve.nacv else false;

    if (comptime @hasField(@TypeOf(clean), "optimize")) {
        clean.optimize = null;
    }

    if (comptime @hasField(@TypeOf(clean), "frequency")) {
        clean.frequency = null;
    }

    if (comptime @hasField(@TypeOf(clean), "response")) {
        clean.response = null;
    }

    if (comptime @hasField(@TypeOf(clean), "write")) {
        clean.write = .{};
    }

    if (!preserve_grad and comptime @hasField(@TypeOf(clean), "gradient")) {
        clean.gradient = null;
    }

    if (!preserve_hess and comptime @hasField(@TypeOf(clean), "hessian")) {
        clean.hessian = null;
    }

    if (!preserve_nacv and comptime @hasField(@TypeOf(clean), "nacv")) {
        clean.nacv = null;
    }

    inline for (.{ "lowdin", "mayer", "mulliken", "wiberg" }) |pop| {
        if (comptime @hasField(@TypeOf(clean), pop)) @field(clean, pop) = false;
    }

    if (comptime @hasField(@TypeOf(clean), "hartree_fock")) {
        clean.hartree_fock = stripSecondaryOptions(clean.hartree_fock, .{});
    }

    return clean;
}

/// Evaluates electronic state energies for a perturbed molecular geometry.
fn getEnergies(comptime T: type, ctx: DerivativeContext(T), runFn: anytype, opt: anytype, perts: anytype, Pg: ?Matrix(T), P: ?*Matrix(T), out: []T) !void {
    std.debug.assert(@typeInfo(@TypeOf(perts)).@"struct".field_names.len % 2 == 0);

    const modified_opt, const sys = .{ stripSecondaryOptions(opt, .{}), ctx.sys };

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

    var res = try runFn(T, ctx.io, modified_opt, sys, Pg, false, ctx.gpa);
    defer res.deinit(ctx.gpa);

    if (P) |out_P| {
        @memcpy(out_P.data, (if (@hasField(@TypeOf(res), "hartree_fock")) res.hartree_fock.P else res.P).data);
    }

    if (res.energy.len <= out.len - 1) {
        const msg = "STATE REQUESTED FOR NUMERICAL DERIVATIVE IS NOT AVAILABLE IN THE ELECTRONIC STRUCTURE CALCULATION";

        std.log.err(msg, .{});

        return error.InvalidInput;
    }

    @memcpy(out, res.energy[0..out.len]);
}
