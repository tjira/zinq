const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;

const printf = @import("read_write.zig").printf;

const A2BOHR = @import("constant.zig").A2BOHR;

pub fn XyzData(comptime T: type) type {
    return struct {
        sym: [][]const u8,
        coords: Matrix(T),

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.coords.deinit(gpa);

            for (self.sym) |sym| {
                gpa.free(sym);
            }

            gpa.free(self.sym);
        }
    };
}

pub fn calculateNumericalGradient(comptime T: type, io: std.Io, runFn: anytype, opt: anytype, log: bool, gpa: Allocator) anyerror!Matrix(T) {
    const sys_path = if (@hasField(@TypeOf(opt), "system")) opt.system else opt.hartree_fock.system;

    var xyz = try readXyz(T, io, sys_path, gpa);
    defer xyz.deinit(gpa);

    var grad = try Matrix(T).initZero(xyz.sym.len, 3, gpa);
    errdefer grad.deinit(gpa);

    const h = @as(T, @floatCast(opt.gradient.?.numeric.step));

    const state = if (@hasField(@TypeOf(opt.gradient.?.numeric), "state")) opt.gradient.?.numeric.state else 0;

    if (log) try std.Io.File.stdout().writeStreamingAll(io, "\n");

    for (0..xyz.sym.len) |i| for (0..3) |j| {
        const d1, const d2 = .{ .{ i * 3 + j, h }, .{ i * 3 + j, -h } };

        const ep = try getE(T, io, runFn, opt, &xyz, d1, state, gpa);
        const em = try getE(T, io, runFn, opt, &xyz, d2, state, gpa);

        grad.ptr(i, j).* = (ep - em) / (2 * h);

        if (log) {
            try printf(io, "DERIVATIVE OF ATOM {d:02} ({s:2}) IN DIRECTION {d}: {d:20.14}\n", .{ i, xyz.sym[i], j, grad.at(i, j) });
        }
    };

    return grad;
}

pub fn calculateNumericalHessian(comptime T: type, io: std.Io, runFn: anytype, opt: anytype, log: bool, gpa: Allocator) anyerror!Matrix(T) {
    const sys_path = if (@hasField(@TypeOf(opt), "system")) opt.system else opt.hartree_fock.system;

    var xyz = try readXyz(T, io, sys_path, gpa);
    defer xyz.deinit(gpa);

    var hess = try Matrix(T).initZero(3 * xyz.sym.len, 3 * xyz.sym.len, gpa);
    errdefer hess.deinit(gpa);

    const h = @as(T, @floatCast(opt.hessian.?.numeric.step));

    const state = if (@hasField(@TypeOf(opt.hessian.?.numeric), "state")) opt.hessian.?.numeric.state else 0;

    if (log) try std.Io.File.stdout().writeStreamingAll(io, "\n");

    const er = try getE(T, io, runFn, opt, &xyz, .{}, state, gpa);

    for (0..3 * xyz.sym.len) |a| {
        const d1, const d2 = .{ .{ a, h }, .{ a, -h } };

        const ep = try getE(T, io, runFn, opt, &xyz, d1, state, gpa);
        const em = try getE(T, io, runFn, opt, &xyz, d2, state, gpa);

        hess.ptr(a, a).* = (ep - 2 * er + em) / (h * h);

        if (log) {
            const params = .{ a / 3, xyz.sym[a / 3], a % 3, a % 3, hess.at(a, a) };

            try printf(io, "SECOND DERIVATIVE OF ATOM {d:02} ({s:2}) IN DIRECTION {d} AND {d}: {d:20.14}\n", params);
        }

        for (a + 1..3 * xyz.sym.len) |b| {
            const epp = try getE(T, io, runFn, opt, &xyz, .{ a, h, b, h }, state, gpa);

            const epm = try getE(T, io, runFn, opt, &xyz, .{ a, h, b, -h }, state, gpa);
            const emp = try getE(T, io, runFn, opt, &xyz, .{ a, -h, b, h }, state, gpa);

            const emm = try getE(T, io, runFn, opt, &xyz, .{ a, -h, b, -h }, state, gpa);

            const val = (epp - epm - emp + emm) / (4 * h * h);

            hess.ptr(a, b).* = val;
            hess.ptr(b, a).* = val;

            if (log) {
                const params = .{ a / 3, xyz.sym[a / 3], a % 3, b % 3, val };

                try printf(io, "SECOND DERIVATIVE OF ATOM {d:02} ({s:2}) IN DIRECTION {d} AND {d}: {d:20.14}\n", params);
            }
        }
    }

    return hess;
}

pub fn readXyz(comptime T: type, io: std.Io, path: []const u8, gpa: Allocator) !XyzData(T) {
    const content = try std.Io.Dir.cwd().readFileAlloc(io, path, gpa, .unlimited);
    defer gpa.free(content);

    var line_it = std.mem.splitScalar(u8, content, '\n');

    const first_line = line_it.next() orelse return error.InvalidFormat;

    const natoms = try std.fmt.parseInt(usize, std.mem.trim(u8, first_line, " \r\t"), 10);

    _ = line_it.next() orelse return error.InvalidFormat;

    var sym, var i: usize = .{ try gpa.alloc([]const u8, natoms), 0 };
    errdefer gpa.free(sym);

    errdefer for (0..i) |j| {
        gpa.free(sym[j]);
    };

    var coords = try Matrix(T).init(natoms, 3, gpa);
    errdefer coords.deinit(gpa);

    while (i < natoms) : (i += 1) {
        var tok_it = std.mem.tokenizeAny(u8, line_it.next() orelse return error.InvalidFormat, " \t\r");

        sym[i] = try gpa.dupe(u8, tok_it.next() orelse return error.InvalidFormat);

        const x_str = tok_it.next() orelse return error.InvalidFormat;
        const y_str = tok_it.next() orelse return error.InvalidFormat;
        const z_str = tok_it.next() orelse return error.InvalidFormat;

        coords.ptr(i, 0).* = try std.fmt.parseFloat(T, x_str) * A2BOHR;
        coords.ptr(i, 1).* = try std.fmt.parseFloat(T, y_str) * A2BOHR;
        coords.ptr(i, 2).* = try std.fmt.parseFloat(T, z_str) * A2BOHR;
    }

    return XyzData(T){ .sym = sym, .coords = coords };
}

pub fn writeXyz(comptime T: type, io: std.Io, path: []const u8, sym: [][]const u8, coords: Matrix(T)) !void {
    var file = try std.Io.Dir.cwd().createFile(io, path, .{});
    defer file.close(io);

    var buffer: [65536]u8 = undefined;
    var writer = file.writer(io, &buffer);

    try writer.interface.print("{d}\n\n", .{sym.len});

    for (0..sym.len) |i| {
        const x = coords.at(i, 0) / A2BOHR;
        const y = coords.at(i, 1) / A2BOHR;
        const z = coords.at(i, 2) / A2BOHR;

        try writer.interface.print("{s} {d:12.8} {d:12.8} {d:12.8}\n", .{ sym[i], x, y, z });
    }

    try writer.interface.flush();
}

fn getE(comptime T: type, io: std.Io, runFn: anytype, opt: anytype, xyz: *XyzData(T), perts: anytype, state: usize, gpa: Allocator) !T {
    std.debug.assert(@typeInfo(@TypeOf(perts)).@"struct".fields.len % 2 == 0);

    var modified_opt, const fname = .{ opt, "molecule.xyz" };

    if (@hasField(@TypeOf(opt), "system")) {
        modified_opt.system = fname;
    }

    if (!@hasField(@TypeOf(opt), "system")) {
        modified_opt.hartree_fock.system = fname;
    }

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

    inline for (0..@typeInfo(@TypeOf(perts)).@"struct".fields.len / 2) |idx| {
        const i_val = perts[idx * 2 + 0];
        const v_val = perts[idx * 2 + 1];

        const i = i_val / 3;
        const j = i_val % 3;

        xyz.coords.ptr(i, j).* += v_val;
    }

    try writeXyz(T, io, fname, xyz.sym, xyz.coords);
    defer std.Io.Dir.cwd().deleteFile(io, fname) catch {};

    var res = try runFn(T, io, modified_opt, false, gpa);
    defer res.deinit(gpa);

    inline for (0..@typeInfo(@TypeOf(perts)).@"struct".fields.len / 2) |idx| {
        const i_val = perts[idx * 2 + 0];
        const v_val = perts[idx * 2 + 1];

        const i = i_val / 3;
        const j = i_val % 3;

        xyz.coords.ptr(i, j).* -= v_val;
    }

    return res.energy[state];
}
