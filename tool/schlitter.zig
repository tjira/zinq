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
const schlitterEntropy = zinq.trajectory_thermodynamics.schlitterEntropy;

const Na = global_variables.Na;
const Eh = global_variables.Eh;

pub fn help() !void {
    try print("USAGE: zinq-schlitter [-s TRAJECTORY_XYZ] [-t TEMPERATURE] [-h]\n", .{});
}

pub fn parse(trajectory: *[]const u8, temperature: *f64, allocator: std.mem.Allocator, h: *bool) !std.process.ArgIterator {
    var argc: usize = 0;
    var argv = try std.process.argsWithAllocator(allocator);
    _ = argv.next();

    while (argv.next()) |arg| {
        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            h.* = true;
            try help();
            return argv;
        }

        if (std.mem.eql(u8, arg, "-s") or std.mem.eql(u8, arg, "--trajectory")) {
            trajectory.* = argv.next() orelse return error.InvalidArgument;
            argc += 1;
        }
        if (std.mem.eql(u8, arg, "-t") or std.mem.eql(u8, arg, "--temperature")) {
            temperature.* = try std.fmt.parseFloat(f64, argv.next() orelse return error.InvalidArgument);
            argc += 1;
        }

        argc += 1;
    }

    return argv;
}

pub fn main() !void {
    var trajectory_path: []const u8 = "trajectory.xyz";
    var temperature: f64 = 300;

    var timer_total = try std.time.Timer.start();
    var h = false;

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();

    defer {
        if (gpa.deinit() == .leak) std.log.err("MEMORY LEAK DETECTED IN THE ALLOCATOR\n", .{});
    }

    var argv = try parse(&trajectory_path, &temperature, allocator, &h);
    defer argv.deinit();
    if (h) return;

    try print("CALCULATING THE SCHLITTER ENTROPY OF A TRAJECTORY - TRAJECTORY: {s}, TEMPERATURE: {d:.3}\n", .{ trajectory_path, temperature });

    {
        const traj = try readTrajectoryFromXYZ(f64, trajectory_path, allocator);
        defer traj.positions.deinit(allocator);
        defer traj.masses.deinit(allocator);

        var timer_entropy = try std.time.Timer.start();
        try print("\nCALCULATING THE SCHLITTER ENTROPY: ", .{});

        const entropy = try schlitterEntropy(f64, traj.positions, traj.masses, temperature, true, allocator);

        try print("{D}\n", .{timer_entropy.read()});

        try print("\nSCHLITTER ENTROPY: {d:.8} Eh/K = {d:.8} J/MOL/K\n", .{ entropy, entropy * Na * Eh });
    }

    try print("\nTOTAL EXECUTION TIME: {D}\n", .{timer_total.read()});
}
