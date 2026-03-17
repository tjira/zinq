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

pub fn parse(system: *[]const u8, basis: *[]const u8, do_overlap: *bool, do_kinetic: *bool, do_nuclear: *bool, do_coulomb: *bool, do_export: *bool, allocator: std.mem.Allocator, h: *bool) !std.process.ArgIterator {
    var argc: usize = 0; var argv = try std.process.argsWithAllocator(allocator); _ = argv.next();

    while (argv.next()) |arg| {

        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {h.* = true; try help(); return argv;}

        if (std.mem.eql(u8, arg, "-b") or std.mem.eql(u8, arg, "--basis"   )) {basis.* = argv.next() orelse return error.InvalidArgument; argc += 1;}
        if (std.mem.eql(u8, arg, "-m") or std.mem.eql(u8, arg, "--molecule")) {system.* = argv.next() orelse return error.InvalidArgument; argc += 1;}

        if (std.mem.eql(u8, arg, "--overlap")) {do_overlap.* = true;}
        if (std.mem.eql(u8, arg, "--kinetic")) {do_kinetic.* = true;}
        if (std.mem.eql(u8, arg, "--nuclear")) {do_nuclear.* = true;}
        if (std.mem.eql(u8, arg, "--coulomb")) {do_coulomb.* = true;}
        if (std.mem.eql(u8, arg, "--export" )) {do_export.* = true;}

        argc += 1;
    }

    return argv;
}

pub fn main() !void {
    var system_path: []const u8 = "molecule.xyz"; var basis_name: []const u8 = "sto-3g"; var do_export = false;

    var do_overlap = false; var do_kinetic = false; var do_nuclear = false; var do_coulomb = false;

    var timer_total = try std.time.Timer.start(); var h = false;

    var gpa = std.heap.GeneralPurposeAllocator(.{}){}; const allocator = gpa.allocator();

    defer {
        if (gpa.deinit() == .leak) std.log.err("MEMORY LEAK DETECTED IN THE ALLOCATOR\n", .{});
    }

    var argv = try parse(&system_path, &basis_name, &do_overlap, &do_kinetic, &do_nuclear, &do_coulomb, &do_export, allocator, &h); defer argv.deinit(); if (h) return;

    try print("INTEGRAL OVER ATOMIC BASIS FUNCTIONS - SYSTEM: {s}, BASIS: {s}\n", .{system_path, basis_name});

    {
        const system = try zinq.classical_particle.read(f64, system_path, 0, 0, allocator); defer system.deinit(allocator);
        var basis = try BasisSet(f64).init(system, basis_name, allocator); defer basis.deinit(allocator);

        var timer_integrals = try std.time.Timer.start(); try print("\nCOMPUTING THE INTEGRALS: ", .{});

        var S: ?RealMatrix(f64) = null; var T: ?RealMatrix(f64) = null; var V: ?RealMatrix(f64) = null; var J: ?RealTensor4(f64) = null;

        if (do_overlap) S = try overlap(f64, basis, 1, allocator);
        if (do_kinetic) T = try kinetic(f64, basis, 1, allocator);
        if (do_coulomb) J = try coulomb(f64, basis, 1, allocator);

        if (do_nuclear) V = try nuclear(f64, system, basis, 1, allocator);

        try print("{D}\n", .{timer_integrals.read()});

        if (do_export) if (S) |s| try exportRealMatrix(f64, "S_AO.mat", s);
        if (do_export) if (T) |t| try exportRealMatrix(f64, "T_AO.mat", t);
        if (do_export) if (V) |v| try exportRealMatrix(f64, "V_AO.mat", v);

        if (do_export) if (J) |j| try exportRealTensorFour(f64, "J_AO.mat", j);

        if (S) |s| s.deinit(allocator);
        if (T) |t| t.deinit(allocator);
        if (V) |v| v.deinit(allocator);
        if (J) |j| j.deinit(allocator);
    }

    try print("\nTOTAL EXECUTION TIME: {D}\n", .{timer_total.read()});
}
