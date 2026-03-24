const std = @import("std");
const zinq = @import("zinq");

const classical_particle = zinq.classical_particle;
const device_write = zinq.device_write;
const global_variables = zinq.global_variables;
const system_alignment = zinq.system_alignment;
const trajectory_thermodynamics = zinq.trajectory_thermodynamics;

const RealMatrix = zinq.real_matrix.RealMatrix;

const exportRealMatrix = zinq.device_write.exportRealMatrix;
const print = zinq.device_write.print;
const printRealMatrix = zinq.device_write.printRealMatrix;
const printClassicalParticleAsMolecule = zinq.device_write.printClassicalParticleAsMolecule;
const readTrajectoryFromXYZ = zinq.system_alignment.readTrajectoryFromXYZ;
const sreEntropy = zinq.trajectory_thermodynamics.sreEntropy;

const Na = global_variables.Na;
const Eh = global_variables.Eh;

pub fn help() !void {
    try print("USAGE: zinq-sre [-s TRAJECTORY_XYZ] [-t TEMPERATURE] [-h]\n", .{});
}

pub fn parse(trajectory: *[]const u8, temperature: *f64, time_step: *f64, allocator: std.mem.Allocator, h: *bool) !std.process.ArgIterator {
    var argc: usize = 0; var argv = try std.process.argsWithAllocator(allocator); _ = argv.next();

    while (argv.next()) |arg| {

        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {h.* = true; try help(); return argv;}

        if (std.mem.eql(u8, arg, "-d") or std.mem.eql(u8, arg, "--dt"         )) {time_step.* = try std.fmt.parseFloat(f64, argv.next() orelse return error.InvalidArgument); argc += 1;}
        if (std.mem.eql(u8, arg, "-s") or std.mem.eql(u8, arg, "--trajectory" )) {trajectory.* = argv.next() orelse return error.InvalidArgument; argc += 1;}
        if (std.mem.eql(u8, arg, "-t") or std.mem.eql(u8, arg, "--temperature")) {temperature.* = try std.fmt.parseFloat(f64, argv.next() orelse return error.InvalidArgument); argc += 1;}

        argc += 1;
    }

    return argv;
}

pub fn main() !void {
    var trajectory_path: []const u8 = "trajectory.xyz"; var temperature: f64 = 300; var time_step: f64 = 10;

    var timer_total = try std.time.Timer.start(); var h = false;

    var gpa = std.heap.GeneralPurposeAllocator(.{}){}; const allocator = gpa.allocator();

    defer {
        if (gpa.deinit() == .leak) std.log.err("MEMORY LEAK DETECTED IN THE ALLOCATOR\n", .{});
    }

    var argv = try parse(&trajectory_path, &temperature, &time_step, allocator, &h); defer argv.deinit(); if (h) return;

    try print("CALCULATING THE SRE ENTROPY OF A TRAJECTORY - TRAJECTORY: {s}, TEMPERATURE: {d:.3}, TIME STEP: {d:.3}\n", .{trajectory_path, temperature, time_step});

    {
        const traj = try readTrajectoryFromXYZ(f64, trajectory_path, allocator); defer traj.positions.deinit(allocator); defer traj.masses.deinit(allocator);

        var timer_entropy = try std.time.Timer.start(); try print("\nCALCULATING THE SRE ENTROPY: ", .{});

        const entropy = try sreEntropy(f64, traj.positions, traj.masses, temperature, time_step, traj.positions.cols - 6, true, allocator);

        try print("{D}\n", .{timer_entropy.read()});

        try print("\nSRE ENTROPY: {d:.8} Eh/K = {d:.8} J/MOL/K\n", .{entropy, entropy * Na * Eh});
    }

    try print("\nTOTAL EXECUTION TIME: {D}\n", .{timer_total.read()});
}
