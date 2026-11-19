//! Numerical nuclear gradients and Hessians via finite differences for molecular geometry optimization.

const std = @import("std");

const libint = @import("cimport.zig").libint;

const Allocator = std.mem.Allocator;

const electronic_structure = @import("electronic_structure.zig");

const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;

const calculateStateSpin = @import("determinant.zig").calculateStateSpin;
const det = @import("linear_algebra.zig").det;
const gemm = @import("linear_algebra.zig").gemm;
const getSymbol = @import("constants.zig").getSymbol;
const mo2ms_coef = @import("integral_transform.zig").mo2ms_coef;
const printf = @import("read_write.zig").printf;

/// Execution context encapsulating molecular system geometry, IO streams, memory allocator, and logging settings.
pub fn DerivativeContext(comptime T: type) type {
    return struct { gpa: Allocator, io: std.Io, log: bool = false, sys: *MolecularSystem(T), Pg: ?Matrix(T) = null, C_mo: ?Matrix(T) = null };
}

/// Computes nuclear gradients of electronic states using central finite differences.
pub fn calculateNumericalGradient(comptime T: type, ctx: DerivativeContext(T), runFn: anytype, opt: anytype, states: []const u32) anyerror![]Matrix(T) {
    const generalized = electronic_structure.isGeneralized(opt);

    var grad = try ctx.gpa.alloc(Matrix(T), states.len);
    errdefer ctx.gpa.free(grad);

    var grad_done: usize = 0;

    errdefer for (0..grad_done) |k| grad[k].deinit(ctx.gpa);

    for (0..states.len) |k| {
        grad[k] = try Matrix(T).initZero(ctx.sys.atoms.len, 3, ctx.gpa);

        grad_done += 1;
    }

    const h = @as(T, @floatCast(opt.gradient.?.numeric.step));

    if (ctx.log) try printf(ctx.io, "\n", .{});

    const nbf = if (generalized) 2 * ctx.sys.nbf else ctx.sys.nbf;

    var Pg = try Matrix(T).initZero(nbf, nbf, ctx.gpa);
    defer Pg.deinit(ctx.gpa);

    var max_state: usize = 0;

    for (states) |s| {
        if (s > max_state) max_state = s;
    }

    const n = max_state + 1;

    const energies = try ctx.gpa.alloc(T, 2 * n);
    defer ctx.gpa.free(energies);

    const ep, const em = .{ energies[0..n], energies[n .. 2 * n] };

    var res_ref = try runFn(T, ctx.io, stripSecondaryOptions(opt, .{}), ctx.sys, ctx.Pg, ctx.C_mo, false, ctx.gpa);
    defer res_ref.deinit(ctx.gpa);

    if (res_ref.energy.len <= max_state) {
        std.log.err("STATE REQUESTED FOR NUMERICAL DERIVATIVE IS NOT AVAILABLE IN THE ELECTRONIC STRUCTURE CALCULATION", .{});

        return error.InvalidInput;
    }

    @memcpy(Pg.data, electronic_structure.getDensity(T, res_ref).data);

    for (0..ctx.sys.atoms.len) |i| for (0..3) |j| {
        try getEnergies(T, ctx, runFn, opt, .{ i * 3 + j, 0 + h }, Pg, ep);
        try getEnergies(T, ctx, runFn, opt, .{ i * 3 + j, 0 - h }, Pg, em);

        for (states, 0..) |s, s_idx| {
            grad[s_idx].ptr(i, j).* = (ep[s] - em[s]) / (2 * h);

            if (ctx.log) {
                try logGradient(T, ctx.io, opt, ctx.sys.atoms, s, i, j, grad[s_idx].at(i, j));
            }
        }
    };

    return grad;
}

/// Computes nuclear Hessians of electronic states using finite differences of energies.
pub fn calculateNumericalHessian(comptime T: type, ctx: DerivativeContext(T), runFn: anytype, opt: anytype, states: []const u32) anyerror![]Matrix(T) {
    const generalized = electronic_structure.isGeneralized(opt);

    var hess = try ctx.gpa.alloc(Matrix(T), states.len);
    errdefer ctx.gpa.free(hess);

    var hess_done: usize = 0;

    errdefer for (0..hess_done) |k| hess[k].deinit(ctx.gpa);

    for (0..states.len) |k| {
        hess[k] = try Matrix(T).initZero(3 * ctx.sys.atoms.len, 3 * ctx.sys.atoms.len, ctx.gpa);

        hess_done += 1;
    }

    const h = @as(T, @floatCast(opt.hessian.?.numeric.step));

    if (ctx.log) try printf(ctx.io, "\n", .{});

    const nbf = if (generalized) 2 * ctx.sys.nbf else ctx.sys.nbf;

    var Pg = try Matrix(T).initZero(nbf, nbf, ctx.gpa);
    defer Pg.deinit(ctx.gpa);

    var max_state: usize = 0;

    for (states) |s| {
        if (s > max_state) max_state = s;
    }

    const n = max_state + 1;

    const energies = try ctx.gpa.alloc(T, 7 * n);
    defer ctx.gpa.free(energies);

    const er, const ep, const em = .{ energies[0..n], energies[n .. 2 * n], energies[2 * n .. 3 * n] };

    const epp, const epm = .{ energies[3 * n .. 4 * n], energies[4 * n .. 5 * n] };
    const emp, const emm = .{ energies[5 * n .. 6 * n], energies[6 * n .. 7 * n] };

    var res_ref = try runFn(T, ctx.io, stripSecondaryOptions(opt, .{}), ctx.sys, ctx.Pg, ctx.C_mo, false, ctx.gpa);
    defer res_ref.deinit(ctx.gpa);

    if (res_ref.energy.len <= max_state) {
        std.log.err("STATE REQUESTED FOR NUMERICAL DERIVATIVE IS NOT AVAILABLE IN THE ELECTRONIC STRUCTURE CALCULATION", .{});

        return error.InvalidInput;
    }

    @memcpy(Pg.data, electronic_structure.getDensity(T, res_ref).data);

    @memcpy(er, res_ref.energy[0..er.len]);

    for (0..3 * ctx.sys.atoms.len) |a| {
        try getEnergies(T, ctx, runFn, opt, .{ a, 0 + h }, Pg, ep);
        try getEnergies(T, ctx, runFn, opt, .{ a, 0 - h }, Pg, em);

        for (states, 0..) |s, s_idx| {
            hess[s_idx].ptr(a, a).* = (ep[s] - 2 * er[s] + em[s]) / (h * h);

            if (ctx.log) {
                try logHessian(T, ctx.io, opt, ctx.sys.atoms, s, a, a, hess[s_idx].at(a, a));
            }
        }

        for (a + 1..3 * ctx.sys.atoms.len) |b| {
            try getEnergies(T, ctx, runFn, opt, .{ a, h, b, h }, Pg, epp);

            try getEnergies(T, ctx, runFn, opt, .{ a, h, b, -h }, Pg, epm);
            try getEnergies(T, ctx, runFn, opt, .{ a, -h, b, h }, Pg, emp);

            try getEnergies(T, ctx, runFn, opt, .{ a, -h, b, -h }, Pg, emm);

            for (states, 0..) |s, s_idx| {
                const val = (epp[s] - epm[s] - emp[s] + emm[s]) / (4 * h * h);

                hess[s_idx].ptr(a, b).* = val;
                hess[s_idx].ptr(b, a).* = val;

                if (ctx.log) {
                    try logHessian(T, ctx.io, opt, ctx.sys.atoms, s, a, b, val);
                }
            }
        }
    }

    return hess;
}

/// Computes nonadiabatic coupling vectors between state pairs using central finite differences.
pub fn calculateNumericalNacv(comptime T: type, ctx: DerivativeContext(T), runFn: anytype, opt: anytype, pairs: []const [2]u32) anyerror![]Matrix(T) {
    const generalized = electronic_structure.isGeneralized(opt);

    var nacvs = try ctx.gpa.alloc(Matrix(T), pairs.len);
    errdefer ctx.gpa.free(nacvs);

    var nacvs_done: usize = 0;

    errdefer for (0..nacvs_done) |k| nacvs[k].deinit(ctx.gpa);

    for (0..pairs.len) |k| {
        nacvs[k] = try Matrix(T).initZero(ctx.sys.atoms.len, 3, ctx.gpa);

        nacvs_done += 1;
    }

    const h = @as(T, @floatCast(opt.nacv.?.numeric.step));

    var max_state: usize = 0;

    for (pairs) |pair| {
        if (pair[0] > max_state) max_state = pair[0];
        if (pair[1] > max_state) max_state = pair[1];
    }

    if (ctx.log) try printf(ctx.io, "\n", .{});

    const nbf = if (generalized) 2 * ctx.sys.nbf else ctx.sys.nbf;

    var Pg = try Matrix(T).initZero(nbf, nbf, ctx.gpa);
    defer Pg.deinit(ctx.gpa);

    const modified_opt = stripSecondaryOptions(opt, .{});

    if (comptime electronic_structure.isMultiState(@TypeOf(modified_opt))) {
        const required = @as(u32, @intCast(max_state + 1));

        if (modified_opt.nstate < required) {
            std.log.err("REQUESTED NACV STATE COUNT IS NOT AVAILABLE IN THE ELECTRONIC STRUCTURE CALCULATION", .{});

            return error.InvalidInput;
        }
    }

    var res_ref = try runFn(T, ctx.io, modified_opt, ctx.sys, ctx.Pg, ctx.C_mo, false, ctx.gpa);
    defer res_ref.deinit(ctx.gpa);

    @memcpy(Pg.data, electronic_structure.getDensity(T, res_ref).data);

    const C_ci_ref = res_ref.C_ci;

    if (C_ci_ref.ncol() <= max_state) {
        std.log.err("STATE REQUESTED FOR NUMERICAL NACV IS NOT AVAILABLE IN THE ELECTRONIC STRUCTURE CALCULATION", .{});

        return error.InvalidInput;
    }

    var dets = try electronic_structure.generateDeterminants(ctx.sys.nel, ctx.sys.nbf, opt, ctx.gpa);

    defer {
        for (0..dets.items.len) |i| ctx.gpa.free(dets.items[i]);

        dets.deinit(ctx.gpa);
    }

    for (pairs) |pair| {
        const s2_0 = calculateStateSpin(T, C_ci_ref, pair[0], dets.items);
        const s2_1 = calculateStateSpin(T, C_ci_ref, pair[1], dets.items);

        if (@abs(s2_0 - s2_1) > 0.5) {
            const params = .{ s2_0, s2_1 };

            std.log.err("NONADIABATIC COUPLING REQUIRES STATES WITH THE SAME SPIN MULTIPLICITY (<S^2> = {d:.2} vs {d:.2})", params);

            return error.InvalidInput;
        }
    }

    var sys_m = try ctx.sys.clone(ctx.gpa);
    defer sys_m.deinit(ctx.gpa);

    var sys_p = try ctx.sys.clone(ctx.gpa);
    defer sys_p.deinit(ctx.gpa);

    const nthreads = if (@hasField(@TypeOf(opt), "nthreads")) opt.nthreads else 1;

    const nsp = 2 * ctx.sys.nbf;

    const C_ref = electronic_structure.getOrbitals(T, res_ref);

    var C_ref_ms = if (generalized) C_ref else try Matrix(T).initZero(nsp, nsp, ctx.gpa);
    defer if (!generalized) C_ref_ms.deinit(ctx.gpa);

    if (!generalized) mo2ms_coef(T, &C_ref_ms, C_ref);

    var S_ms = try Matrix(T).init(nsp, nsp, ctx.gpa);
    defer S_ms.deinit(ctx.gpa);

    var tmp = try Matrix(T).init(nsp, nsp, ctx.gpa);
    defer tmp.deinit(ctx.gpa);

    var C_m_buf = if (!generalized) try Matrix(T).initZero(nsp, nsp, ctx.gpa) else undefined;
    defer if (!generalized) C_m_buf.deinit(ctx.gpa);

    var C_p_buf = if (!generalized) try Matrix(T).initZero(nsp, nsp, ctx.gpa) else undefined;
    defer if (!generalized) C_p_buf.deinit(ctx.gpa);

    const couplings = try ctx.gpa.alloc(T, pairs.len);
    defer ctx.gpa.free(couplings);

    for (0..ctx.sys.atoms.len) |i| for (0..3) |j| {
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

        var res_m = try runFn(T, ctx.io, modified_opt, &sys_m, Pg, C_ref, false, ctx.gpa);
        defer res_m.deinit(ctx.gpa);

        var res_p = try runFn(T, ctx.io, modified_opt, &sys_p, Pg, C_ref, false, ctx.gpa);
        defer res_p.deinit(ctx.gpa);

        var S_as = try sys_m.overlapBetweenSpin(sys_p, nthreads, ctx.gpa);
        defer S_as.deinit(ctx.gpa);

        const C_m = electronic_structure.getOrbitals(T, res_m);
        const C_p = electronic_structure.getOrbitals(T, res_p);

        var C_m_ms = if (generalized) C_m else C_m_buf;

        if (!generalized) mo2ms_coef(T, &C_m_ms, C_m);

        var C_p_ms = if (generalized) C_p else C_p_buf;

        if (!generalized) mo2ms_coef(T, &C_p_ms, C_p);

        alignDisplacedPhases(T, .{ &C_m_ms, &C_p_ms }, C_ref_ms, .{ &res_m.C_ci, &res_p.C_ci }, C_ci_ref, dets.items);

        gemm(T, &tmp, S_as, C_p_ms, 1, 0, false, false);

        gemm(T, &S_ms, C_m_ms, tmp, 1, 0, true, false);

        try accumulateNacv(T, couplings, S_ms, .{ res_m.C_ci, res_p.C_ci }, dets.items, pairs, ctx.gpa);

        for (0..pairs.len) |p_idx| {
            nacvs[p_idx].ptr(i, j).* = couplings[p_idx] / (4 * h);
        }

        if (ctx.log) {
            const sym = try getSymbol(ctx.sys.atoms[i]);

            const fmt = "NONADIABATIC COUPLING STATES {d}/{d} FOR ATOM {d:02} ({s:2}) IN DIRECTION {s}: {d:20.14}\n";

            const dirs = [_][]const u8{ "'x'", "'y'", "'z'" };

            for (pairs, 0..) |pair, p_idx| {
                try printf(ctx.io, fmt, .{ pair[0], pair[1], i, sym, dirs[j], nacvs[p_idx].at(i, j) });
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

    inline for (.{ "frequency", "optimize", "response" }) |field| {
        if (comptime @hasField(@TypeOf(clean), field)) @field(clean, field) = null;
    }

    inline for (.{ "gradient", "hessian", "nacv" }) |field| {
        const keep = if (comptime @hasField(@TypeOf(preserve), field)) @field(preserve, field) else false;

        if (!keep and comptime @hasField(@TypeOf(clean), field)) @field(clean, field) = null;
    }

    inline for (.{ "lowdin", "mayer", "mulliken", "wiberg" }) |pop| {
        if (comptime @hasField(@TypeOf(clean), pop)) @field(clean, pop) = false;
    }

    if (comptime @hasField(@TypeOf(clean), "write")) clean.write = .{};

    if (comptime @hasField(@TypeOf(clean), "hartree_fock")) {
        clean.hartree_fock = stripSecondaryOptions(clean.hartree_fock, .{});
    }

    return clean;
}

/// Accumulates nonadiabatic coupling vector elements by contracting determinant overlaps with CI coefficients.
fn accumulateNacv(comptime T: type, couplings: []T, S_MS: Matrix(T), C_ci: [2]Matrix(T), dets: []const []const usize, pairs: []const [2]u32, gpa: Allocator) !void {
    const nel = dets[0].len;

    var M = try Matrix(T).init(nel, nel, gpa);
    defer M.deinit(gpa);

    const ipiv = try gpa.alloc(i32, nel);
    defer gpa.free(ipiv);

    @memset(couplings, 0);

    for (0..dets.len) |a| for (0..dets.len) |b| {
        var any_active = false;

        for (pairs) |pair| {
            const C_0011 = C_ci[0].at(a, pair[0]) * C_ci[1].at(b, pair[1]);
            const C_0110 = C_ci[0].at(a, pair[1]) * C_ci[1].at(b, pair[0]);

            const diff = C_0011 - C_0110;

            if (@abs(diff) >= 1e-15) {
                any_active = true;

                break;
            }
        }

        if (!any_active) continue;

        for (0..nel) |r| for (0..nel) |c| {
            M.ptr(r, c).* = S_MS.at(dets[a][r], dets[b][c]);
        };

        const d = try det(T, &M, ipiv);

        for (pairs, 0..) |pair, p_idx| {
            const C_0011 = C_ci[0].at(a, pair[0]) * C_ci[1].at(b, pair[1]);
            const C_0110 = C_ci[0].at(a, pair[1]) * C_ci[1].at(b, pair[0]);

            const diff = C_0011 - C_0110;

            couplings[p_idx] += diff * d;
        }
    };
}

/// Aligns spin-orbital and CI state phase signs across displaced geometries.
fn alignDisplacedPhases(comptime T: type, C_ms: [2]*Matrix(T), C_ref: Matrix(T), C_ci: [2]*Matrix(T), C_ci_ref: Matrix(T), dets: []const []const usize) void {
    const nsp = C_ref.nrow();

    for (0..2) |idx| {
        for (0..nsp) |p| {
            var d: T = 0;

            for (0..nsp) |mu| {
                d += C_ms[idx].at(mu, p) * C_ref.at(mu, p);
            }

            if (d < 0) {
                for (0..nsp) |mu| {
                    C_ms[idx].ptr(mu, p).* = -C_ms[idx].at(mu, p);
                }

                for (0..dets.len) |a| {
                    if (std.mem.indexOfScalar(usize, dets[a], p) != null) {
                        for (0..C_ci[idx].ncol()) |k| {
                            C_ci[idx].ptr(a, k).* = -C_ci[idx].at(a, k);
                        }
                    }
                }
            }
        }

        for (0..C_ci_ref.ncol()) |k| {
            var state_dot: T = 0;

            for (0..C_ci_ref.nrow()) |r| {
                state_dot += C_ci[idx].at(r, k) * C_ci_ref.at(r, k);
            }

            if (state_dot < 0) {
                for (0..C_ci[idx].nrow()) |r| {
                    C_ci[idx].ptr(r, k).* = -C_ci[idx].at(r, k);
                }
            }
        }
    }
}

/// Evaluates electronic state energies for a perturbed molecular geometry.
fn getEnergies(comptime T: type, ctx: DerivativeContext(T), runFn: anytype, opt: anytype, perts: anytype, Pg: ?Matrix(T), out: []T) !void {
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

    var res = try runFn(T, ctx.io, modified_opt, sys, Pg, ctx.C_mo, false, ctx.gpa);
    defer res.deinit(ctx.gpa);

    if (res.energy.len <= out.len - 1) {
        const msg = "STATE REQUESTED FOR NUMERICAL DERIVATIVE IS NOT AVAILABLE IN THE ELECTRONIC STRUCTURE CALCULATION";

        std.log.err(msg, .{});

        return error.InvalidInput;
    }

    @memcpy(out, res.energy[0..out.len]);
}

/// Prints numerical nuclear gradient components for single-state or multi-state calculations.
fn logGradient(comptime T: type, io: std.Io, opt: anytype, atoms: []const i32, s: u32, i: usize, j: usize, val: T) !void {
    const sym = try getSymbol(atoms[i]);

    const dirs = &[_][]const u8{ "'x'", "'y'", "'z'" };

    if (comptime electronic_structure.isMultiState(@TypeOf(opt))) {
        const fmt = "NUCLEAR GRADIENT STATE {d} FOR ATOM {d:02} ({s:2}) IN DIRECTION {s}: {d:20.14}\n";

        try printf(io, fmt, .{ s, i, sym, dirs[j], val });
    }

    if (!comptime electronic_structure.isMultiState(@TypeOf(opt))) {
        const fmt = "NUCLEAR GRADIENT FOR ATOM {d:02} ({s:2}) IN DIRECTION {s}: {d:20.14}\n";

        try printf(io, fmt, .{ i, sym, dirs[j], val });
    }
}

/// Prints numerical nuclear Hessian components for single-state or multi-state calculations.
fn logHessian(comptime T: type, io: std.Io, opt: anytype, atoms: []const i32, s: u32, a: usize, b: usize, val: T) !void {
    const at_a, const at_b = .{ a / 3, b / 3 };

    const sym_a, const sym_b = .{ try getSymbol(atoms[at_a]), try getSymbol(atoms[at_b]) };

    const dirs = &[_][]const u8{ "'x'", "'y'", "'z'" };

    if (comptime electronic_structure.isMultiState(@TypeOf(opt))) {
        const fmt = "NUCLEAR HESSIAN STATE {d} FOR ATOM {d:02} ({s:2}) DIR {s} AND ATOM {d:02} ({s:2}) DIR {s}: {d:20.14}\n";

        try printf(io, fmt, .{ s, at_a, sym_a, dirs[a % 3], at_b, sym_b, dirs[b % 3], val });
    }

    if (!comptime electronic_structure.isMultiState(@TypeOf(opt))) {
        const fmt = "NUCLEAR HESSIAN FOR ATOM {d:02} ({s:2}) DIR {s} AND ATOM {d:02} ({s:2}) DIR {s}: {d:20.14}\n";

        try printf(io, fmt, .{ at_a, sym_a, dirs[a % 3], at_b, sym_b, dirs[b % 3], val });
    }
}
