const std = @import("std");
const zinq = @import("zinq");

const BasisSet = zinq.basis_set.BasisSet;
const RealMatrix = zinq.real_matrix.RealMatrix;
const RealTensor4 = zinq.real_tensor_four.RealTensor4;

const exportRealMatrix = zinq.device_write.exportRealMatrix;
const exportRealTensorFour = zinq.device_write.exportRealTensorFour;
const overlap = zinq.molecular_integrals.overlap;
const nuclear = zinq.molecular_integrals.nuclear;
const kinetic = zinq.molecular_integrals.kinetic;
const coulomb = zinq.molecular_integrals.coulomb;
const print = zinq.device_write.print;

pub fn help(io: std.Io) !void {
    try print(io, "USAGE: zinq-molint [-b BASIS] -m [MOLECULE] [-h]\n", .{});
}

pub fn parse(io: std.Io, system: *[]const u8, basis: *[]const u8, overlap_path: *?[]const u8, kinetic_path: *?[]const u8, nuclear_path: *?[]const u8, coulomb_path: *?[]const u8, args: []const []const u8, h: *bool) !void {
    var i: usize = 0;

    while (i < args.len) : (i += 1) {
        const arg = args[i];

        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            h.* = true;
            try help(io);
            return;
        }

        if (std.mem.eql(u8, arg, "-b") or std.mem.eql(u8, arg, "--basis")) {
            i += 1;
            basis.* = if (i < args.len) args[i] else return error.InvalidArgument;
            continue;
        }
        if (std.mem.eql(u8, arg, "-j") or std.mem.eql(u8, arg, "--coulomb")) {
            i += 1;
            coulomb_path.* = if (i < args.len) args[i] else return error.InvalidArgument;
            continue;
        }
        if (std.mem.eql(u8, arg, "-m") or std.mem.eql(u8, arg, "--molecule")) {
            i += 1;
            system.* = if (i < args.len) args[i] else return error.InvalidArgument;
            continue;
        }
        if (std.mem.eql(u8, arg, "-s") or std.mem.eql(u8, arg, "--overlap")) {
            i += 1;
            overlap_path.* = if (i < args.len) args[i] else return error.InvalidArgument;
            continue;
        }
        if (std.mem.eql(u8, arg, "-t") or std.mem.eql(u8, arg, "--kinetic")) {
            i += 1;
            kinetic_path.* = if (i < args.len) args[i] else return error.InvalidArgument;
            continue;
        }
        if (std.mem.eql(u8, arg, "-v") or std.mem.eql(u8, arg, "--nuclear")) {
            i += 1;
            nuclear_path.* = if (i < args.len) args[i] else return error.InvalidArgument;
            continue;
        }
    }
}

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    const io = init.io;

    var system_path: []const u8 = "molecule.xyz";
    var basis_name: []const u8 = "sto-3g";

    var overlap_path: ?[]const u8 = null;
    var kinetic_path: ?[]const u8 = null;
    var nuclear_path: ?[]const u8 = null;
    var coulomb_path: ?[]const u8 = null;

    var timer_total = std.Io.Timestamp.now(io, .real);
    var h = false;

    const args = try init.minimal.args.toSlice(init.arena.allocator());

    try parse(io, &system_path, &basis_name, &overlap_path, &kinetic_path, &nuclear_path, &coulomb_path, args[1..], &h);

    if (h) return;

    try print(io, "INTEGRAL OVER ATOMIC BASIS FUNCTIONS - SYSTEM: {s}, BASIS: {s}\n", .{ system_path, basis_name });

    {
        const system = try zinq.classical_particle.read(f64, io, system_path, 0, 0, allocator);
        defer system.deinit(allocator);

        var basis = try BasisSet(f64).init(system, basis_name, allocator);
        defer basis.deinit(allocator);

        var timer_integrals = std.Io.Timestamp.now(io, .real);
        try print(io, "\nCOMPUTING THE INTEGRALS: ", .{});

        var S: ?RealMatrix(f64) = null;
        var T: ?RealMatrix(f64) = null;
        var V: ?RealMatrix(f64) = null;
        var J: ?RealTensor4(f64) = null;

        if (coulomb_path) |_| J = try coulomb(f64, system, basis, 1, allocator);
        if (kinetic_path) |_| T = try kinetic(f64, system, basis, 1, allocator);
        if (nuclear_path) |_| V = try nuclear(f64, system, basis, 1, allocator);
        if (overlap_path) |_| S = try overlap(f64, system, basis, 1, allocator);

        try print(io, "{f}\n", .{timer_integrals.untilNow(io, .real)});

        if (overlap_path) |path| if (path.len > 0) {
            try exportRealMatrix(f64, io, path, S.?);
        };
        if (kinetic_path) |path| if (path.len > 0) {
            try exportRealMatrix(f64, io, path, T.?);
        };
        if (nuclear_path) |path| if (path.len > 0) {
            try exportRealMatrix(f64, io, path, V.?);
        };

        if (coulomb_path) |path| if (path.len > 0) {
            try exportRealTensorFour(f64, io, path, J.?);
        };

        if (S) |s| s.deinit(allocator);
        if (T) |t| t.deinit(allocator);
        if (V) |v| v.deinit(allocator);
        if (J) |j| j.deinit(allocator);
    }

    try print(io, "\nTOTAL EXECUTION TIME: {f}\n", .{timer_total.untilNow(io, .real)});
}
