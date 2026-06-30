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

    const h, var modified_opt = .{ @as(T, @floatCast(opt.gradient.?.numeric.step)), opt };

    if (@hasField(@TypeOf(opt), "system")) {
        modified_opt.system = "molecule.xyz";
    }

    if (!@hasField(@TypeOf(opt), "system")) {
        modified_opt.hartree_fock.system = "molecule.xyz";
    }

    modified_opt.gradient = null;

    const state = if (@hasField(@TypeOf(opt.gradient.?.numeric), "state")) opt.gradient.?.numeric.state else 0;

    if (log) try std.Io.File.stdout().writeStreamingAll(io, "\n");

    for (0..xyz.sym.len) |i| for (0..3) |j| {
        const original_value = xyz.coords.at(i, j);

        xyz.coords.ptr(i, j).* = original_value + h;

        try writeXyz(T, io, "molecule.xyz", xyz.sym, xyz.coords);

        var res_plus = try runFn(T, io, modified_opt, false, gpa);
        defer res_plus.deinit(gpa);

        const E_plus = res_plus.energy[state];

        xyz.coords.ptr(i, j).* = original_value - h;

        try writeXyz(T, io, "molecule.xyz", xyz.sym, xyz.coords);

        var res_minus = try runFn(T, io, modified_opt, false, gpa);
        defer res_minus.deinit(gpa);

        const E_minus = res_minus.energy[state];

        xyz.coords.ptr(i, j).* = original_value;

        grad.ptr(i, j).* = (E_plus - E_minus) / (2 * h);

        if (log) {
            try printf(io, "DERIVATIVE OF ATOM {d:02} ({s:2}) IN DIRECTION {d}: {d:20.14}\n", .{ i, xyz.sym[i], j, grad.at(i, j) });
        }
    };

    std.Io.Dir.cwd().deleteFile(io, "molecule.xyz") catch {};

    return grad;
}

pub fn readXyz(comptime T: type, io: std.Io, path: []const u8, gpa: Allocator) !XyzData(T) {
    const content = try std.Io.Dir.cwd().readFileAlloc(io, path, gpa, .unlimited);
    defer gpa.free(content);

    var line_it = std.mem.splitScalar(u8, content, '\n');

    const first_line = line_it.next() orelse return error.InvalidFormat;

    const natoms = try std.fmt.parseInt(usize, std.mem.trim(u8, first_line, " \r\t"), 10);

    _ = line_it.next() orelse return error.InvalidFormat;

    var sym, var i: usize = .{try gpa.alloc([]const u8, natoms), 0};
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
