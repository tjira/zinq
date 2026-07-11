const std = @import("std");

const libint = @import("cimport.zig").libint;

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;

const exportIfBuiltin = @import("molecular_integrals.zig").exportIfBuiltin;
const getSymbol = @import("constant.zig").getSymbol;
const printf = @import("read_write.zig").printf;

const A2BOHR = @import("constant.zig").A2BOHR;

pub fn calculateNumericalGradient(comptime T: type, io: std.Io, runFn: anytype, opt: anytype, sys: *MolecularSystem(T), log: bool, gpa: Allocator) anyerror!Matrix(T) {
    var grad = try Matrix(T).initZero(sys.atoms.len, 3, gpa);
    errdefer grad.deinit(gpa);

    const h = @as(T, @floatCast(opt.gradient.?.numeric.step));

    const state = if (@hasField(@TypeOf(opt.gradient.?.numeric), "state")) opt.gradient.?.numeric.state else 0;

    if (log) try std.Io.File.stdout().writeStreamingAll(io, "\n");

    for (0..sys.atoms.len) |i| for (0..3) |j| {
        const d1, const d2 = .{ .{ i * 3 + j, h }, .{ i * 3 + j, -h } };

        const ep = try getE(T, io, runFn, opt, sys, d1, state, gpa);
        const em = try getE(T, io, runFn, opt, sys, d2, state, gpa);

        grad.ptr(i, j).* = (ep - em) / (2 * h);

        if (log) {
            const sym = try getSymbol(sys.atoms[i]);

            try printf(io, "DERIVATIVE OF ATOM {d:02} ({s:2}) IN DIRECTION {d}: {d:20.14}\n", .{ i, sym, j, grad.at(i, j) });
        }
    };

    return grad;
}

pub fn calculateNumericalHessian(comptime T: type, io: std.Io, runFn: anytype, opt: anytype, sys: *MolecularSystem(T), log: bool, gpa: Allocator) anyerror!Matrix(T) {
    var hess = try Matrix(T).initZero(3 * sys.atoms.len, 3 * sys.atoms.len, gpa);
    errdefer hess.deinit(gpa);

    const h = @as(T, @floatCast(opt.hessian.?.numeric.step));

    const state = if (@hasField(@TypeOf(opt.hessian.?.numeric), "state")) opt.hessian.?.numeric.state else 0;

    if (log) try std.Io.File.stdout().writeStreamingAll(io, "\n");

    const er = try getE(T, io, runFn, opt, sys, .{}, state, gpa);

    for (0..3 * sys.atoms.len) |a| {
        const d1, const d2 = .{ .{ a, h }, .{ a, -h } };

        const ep = try getE(T, io, runFn, opt, sys, d1, state, gpa);
        const em = try getE(T, io, runFn, opt, sys, d2, state, gpa);

        hess.ptr(a, a).* = (ep - 2 * er + em) / (h * h);

        if (log) {
            const sym = try getSymbol(sys.atoms[a / 3]);

            const params = .{ a / 3, sym, a % 3, a % 3, hess.at(a, a) };

            try printf(io, "SECOND DERIVATIVE OF ATOM {d:02} ({s:2}) IN DIRECTION {d} AND {d}: {d:20.14}\n", params);
        }

        for (a + 1..3 * sys.atoms.len) |b| {
            const epp = try getE(T, io, runFn, opt, sys, .{ a, h, b, h }, state, gpa);

            const epm = try getE(T, io, runFn, opt, sys, .{ a, h, b, -h }, state, gpa);
            const emp = try getE(T, io, runFn, opt, sys, .{ a, -h, b, h }, state, gpa);

            const emm = try getE(T, io, runFn, opt, sys, .{ a, -h, b, -h }, state, gpa);

            const val = (epp - epm - emp + emm) / (4 * h * h);

            hess.ptr(a, b).* = val;
            hess.ptr(b, a).* = val;

            if (log) {
                const sym = try getSymbol(sys.atoms[a / 3]);

                const params = .{ a / 3, sym, a % 3, b % 3, val };

                try printf(io, "SECOND DERIVATIVE OF ATOM {d:02} ({s:2}) IN DIRECTION {d} AND {d}: {d:20.14}\n", params);
            }
        }
    }

    return hess;
}

fn getE(comptime T: type, io: std.Io, runFn: anytype, opt: anytype, sys: *MolecularSystem(T), perts: anytype, state: usize, gpa: Allocator) !T {
    std.debug.assert(@typeInfo(@TypeOf(perts)).@"struct".fields.len % 2 == 0);

    var modified_opt = opt;

    modified_opt.gradient, modified_opt.hessian = .{ null, null };

    if (comptime @hasField(@TypeOf(modified_opt), "hartree_fock")) {
        modified_opt.hartree_fock.gradient, modified_opt.hartree_fock.hessian = .{ null, null };
    }

    if (comptime @hasField(@TypeOf(modified_opt), "mulliken")) {
        modified_opt.mulliken = false;
    }

    if (comptime @hasField(@TypeOf(modified_opt), "hartree_fock")) {
        modified_opt.hartree_fock.mulliken = false;
    }

    var original_coors: [@typeInfo(@TypeOf(perts)).@"struct".fields.len / 2]T = undefined;

    inline for (0..@typeInfo(@TypeOf(perts)).@"struct".fields.len / 2) |idx| {
        original_coors[idx] = sys.coors[perts[idx * 2]];
        sys.coors[perts[idx * 2]] += perts[idx * 2 + 1];
    }

    libint.libint_update_coords(sys.ptr, sys.coors.ptr);

    defer {
        inline for (0..@typeInfo(@TypeOf(perts)).@"struct".fields.len / 2) |idx| {
            sys.coors[perts[idx * 2 + 0]] = original_coors[idx];
        }

        libint.libint_update_coords(sys.ptr, sys.coors.ptr);
    }

    var res = try runFn(T, io, modified_opt, sys, false, gpa);
    defer res.deinit(gpa);

    if (res.energy.len <= state) {
        std.log.err("STATE REQUESTED FOR NUMERICAL DERIVATIVE IS NOT AVAILABLE IN THE ELECTRONIC STRUCTURE CALCULATION", .{});

        return error.InvalidInput;
    }

    return res.energy[state];
}
