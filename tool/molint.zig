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

pub fn help() !void {
    try print("USAGE: zinq-molint [-b BASIS] -m [MOLECULE] [-h]\n", .{});
}

pub fn parse(system: *[]const u8, basis: *[]const u8, overlap_path: *?[]const u8, kinetic_path: *?[]const u8, nuclear_path: *?[]const u8, coulomb_path: *?[]const u8, allocator: std.mem.Allocator, h: *bool) !std.process.ArgIterator {
    var argc: usize = 0;
    var argv = try std.process.argsWithAllocator(allocator);
    _ = argv.next();

    while (argv.next()) |arg| {
        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            h.* = true;
            try help();
            return argv;
        }

        if (std.mem.eql(u8, arg, "-b") or std.mem.eql(u8, arg, "--basis")) {
            basis.* = argv.next() orelse return error.InvalidArgument;
            argc += 1;
        }
        if (std.mem.eql(u8, arg, "-j") or std.mem.eql(u8, arg, "--coulomb")) {
            coulomb_path.* = argv.next() orelse return error.InvalidArgument;
            argc += 1;
        }
        if (std.mem.eql(u8, arg, "-m") or std.mem.eql(u8, arg, "--molecule")) {
            system.* = argv.next() orelse return error.InvalidArgument;
            argc += 1;
        }
        if (std.mem.eql(u8, arg, "-s") or std.mem.eql(u8, arg, "--overlap")) {
            overlap_path.* = argv.next() orelse return error.InvalidArgument;
            argc += 1;
        }
        if (std.mem.eql(u8, arg, "-t") or std.mem.eql(u8, arg, "--kinetic")) {
            kinetic_path.* = argv.next() orelse return error.InvalidArgument;
            argc += 1;
        }
        if (std.mem.eql(u8, arg, "-v") or std.mem.eql(u8, arg, "--nuclear")) {
            nuclear_path.* = argv.next() orelse return error.InvalidArgument;
            argc += 1;
        }

        argc += 1;
    }

    return argv;
}

pub fn main() !void {
    var system_path: []const u8 = "molecule.xyz";
    var basis_name: []const u8 = "sto-3g";

    var overlap_path: ?[]const u8 = null;
    var kinetic_path: ?[]const u8 = null;
    var nuclear_path: ?[]const u8 = null;
    var coulomb_path: ?[]const u8 = null;

    var timer_total = try std.time.Timer.start();
    var h = false;

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();

    defer {
        if (gpa.deinit() == .leak) std.log.err("MEMORY LEAK DETECTED IN THE ALLOCATOR\n", .{});
    }

    var argv = try parse(&system_path, &basis_name, &overlap_path, &kinetic_path, &nuclear_path, &coulomb_path, allocator, &h);
    defer argv.deinit();
    if (h) return;

    try print("INTEGRAL OVER ATOMIC BASIS FUNCTIONS - SYSTEM: {s}, BASIS: {s}\n", .{ system_path, basis_name });

    {
        const system = try zinq.classical_particle.read(f64, system_path, 0, 0, allocator);
        defer system.deinit(allocator);

        var basis = try BasisSet(f64).init(system, basis_name, allocator);
        defer basis.deinit(allocator);

        var timer_integrals = try std.time.Timer.start();
        try print("\nCOMPUTING THE INTEGRALS: ", .{});

        var S: ?RealMatrix(f64) = null;
        var T: ?RealMatrix(f64) = null;
        var V: ?RealMatrix(f64) = null;
        var J: ?RealTensor4(f64) = null;

        if (coulomb_path) |_| J = try coulomb(f64, system, basis, 1, allocator);
        if (kinetic_path) |_| T = try kinetic(f64, system, basis, 1, allocator);
        if (nuclear_path) |_| V = try nuclear(f64, system, basis, 1, allocator);
        if (overlap_path) |_| S = try overlap(f64, system, basis, 1, allocator);

        try print("{D}\n", .{timer_integrals.read()});

        if (overlap_path) |path| if (path.len > 0) {
            try exportRealMatrix(f64, path, S.?);
        };
        if (kinetic_path) |path| if (path.len > 0) {
            try exportRealMatrix(f64, path, T.?);
        };
        if (nuclear_path) |path| if (path.len > 0) {
            try exportRealMatrix(f64, path, V.?);
        };

        if (coulomb_path) |path| if (path.len > 0) {
            try exportRealTensorFour(f64, path, J.?);
        };

        if (S) |s| s.deinit(allocator);
        if (T) |t| t.deinit(allocator);
        if (V) |v| v.deinit(allocator);
        if (J) |j| j.deinit(allocator);
    }

    try print("\nTOTAL EXECUTION TIME: {D}\n", .{timer_total.read()});
}
