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

pub fn help(io: std.Io) !void {
    try print(io, "USAGE: zinq-schlitter [-s TRAJECTORY_XYZ] [-t TEMPERATURE] [-h]\n", .{});
}

pub fn parse(io: std.Io, trajectory: *[]const u8, temperature: *f64, args: []const []const u8, h: *bool) !void {
    var i: usize = 0;

    while (i < args.len) : (i += 1) {
        const arg = args[i];

        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            h.* = true;
            try help(io);
            return;
        }

        if (std.mem.eql(u8, arg, "-s") or std.mem.eql(u8, arg, "--trajectory")) {
            i += 1;
            trajectory.* = if (i < args.len) args[i] else return error.InvalidArgument;
            continue;
        }
        if (std.mem.eql(u8, arg, "-t") or std.mem.eql(u8, arg, "--temperature")) {
            i += 1;
            temperature.* = try std.fmt.parseFloat(f64, if (i < args.len) args[i] else return error.InvalidArgument);
            continue;
        }
    }
}

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    const io = init.io;

    var trajectory_path: []const u8 = "trajectory.xyz";
    var temperature: f64 = 300;

    var timer_total = std.Io.Timestamp.now(io, .real);
    var h = false;

    const args = try init.minimal.args.toSlice(init.arena.allocator());

    try parse(io, &trajectory_path, &temperature, args[1..], &h);

    if (h) return;

    try print(io, "CALCULATING THE SCHLITTER ENTROPY OF A TRAJECTORY - TRAJECTORY: {s}, TEMPERATURE: {d:.3}\n", .{ trajectory_path, temperature });

    {
        const traj = try readTrajectoryFromXYZ(f64, io, trajectory_path, allocator);
        defer traj.positions.deinit(allocator);
        defer traj.masses.deinit(allocator);

        var timer_entropy = std.Io.Timestamp.now(io, .real);
        try print(io, "\nCALCULATING THE SCHLITTER ENTROPY: ", .{});

        const entropy = try schlitterEntropy(f64, traj.positions, traj.masses, temperature, true, allocator);

        try print(io, "{f}\n", .{timer_entropy.untilNow(io, .real)});

        try print(io, "\nSCHLITTER ENTROPY: {d:.8} Eh/K = {d:.8} J/MOL/K\n", .{ entropy, entropy * Na * Eh });
    }

    try print(io, "\nTOTAL EXECUTION TIME: {f}\n", .{timer_total.untilNow(io, .real)});
}
