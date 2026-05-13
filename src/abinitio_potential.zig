//! File with ab initio potential struct and functions.

const std = @import("std");

const bias_potential = @import("bias_potential.zig");
const classical_particle = @import("classical_particle.zig");
const device_read = @import("device_read.zig");
const device_write = @import("device_write.zig");
const global_variables = @import("global_variables.zig");
const prcess = @import("process.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const BiasPotential = bias_potential.BiasPotential;
const ClassicalParticle = classical_particle.ClassicalParticle;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const executeCommand = prcess.executeCommand;
const readRealMatrix = device_read.readRealMatrix;

const WRITE_BUFFER_SIZE = global_variables.WRITE_BUFFER_SIZE;

/// Struct holding parameters for the multidimensional ab initio potential.
pub fn AbInitioPotential(comptime T: type) type {
    return struct {
        command: []const u8,
        states: usize = 1,

        /// Evaluate the adiabatic potential energy matrix at given system state and time.
        pub fn evaluateAdiabatic(self: @This(), io: std.Io, adiabatic_potential: *RealMatrix(T), dir: std.Io.Dir, allocator: std.mem.Allocator) !void {
            const dirname = try dir.realPathFileAlloc(io, ".", allocator);
            defer allocator.free(dirname);

            const path = try std.mem.concat(allocator, u8, &.{ dirname, "/ENERGY.mat" });
            defer allocator.free(path);

            const ENERGY = try readRealMatrix(T, io, path, allocator);
            defer ENERGY.deinit(allocator);

            if (ENERGY.cols != 1) {
                std.log.err("THE ENERHY FILE 'ENERGY.mat' MUST BE A COLUMN VECTOR", .{});

                return error.InvalidInput;
            }

            if (ENERGY.rows != self.states) {
                std.log.err("THE NUMBER OF ROWS IN THE ENERGY FILE 'ENERGY.mat' MUST BE EQUAL TO THE NUMBER OF STATES SPECIFIED IN THE AB INITIO POTENTIAL", .{});

                return error.InvalidInput;
            }

            adiabatic_potential.fill(0);

            for (0..ENERGY.rows) |i| adiabatic_potential.ptr(i, i).* = ENERGY.at(i, 0);
        }

        /// Comptime adiabatic force evaluation. The evaluateAdiabatic function needs to be run first.
        pub fn forceAdiabatic(self: @This(), io: std.Io, i: usize, position: RealVector(T), _: T, state: usize, bias: ?BiasPotential(T), dir: std.Io.Dir, allocator: std.mem.Allocator) !T {
            const dirname = try dir.realPathFileAlloc(io, ".", allocator);
            defer allocator.free(dirname);

            const path_gradient = try std.mem.concat(allocator, u8, &.{ dirname, "/GRADIENT.mat" });
            defer allocator.free(path_gradient);

            const path_energy = try std.mem.concat(allocator, u8, &.{ dirname, "/ENERGY.mat" });
            defer allocator.free(path_energy);

            const ENERGY = try readRealMatrix(T, io, path_energy, allocator);
            defer ENERGY.deinit(allocator);

            const GRADIENT = try readRealMatrix(T, io, path_gradient, allocator);
            defer GRADIENT.deinit(allocator);

            if (ENERGY.cols != 1) {
                std.log.err("THE ENERGY FILE 'ENERGY.mat' MUST BE A COLUMN VECTOR", .{});

                return error.InvalidInput;
            }

            if (GRADIENT.cols != 1) {
                std.log.err("THE GRADIENT FILE 'GRADIENT.mat' MUST BE A COLUMN VECTOR", .{});

                return error.InvalidInput;
            }

            if (ENERGY.rows != self.states) {
                std.log.err("THE NUMBER OF ROWS IN THE ENERGY FILE 'ENERGY.mat' MUST BE EQUAL TO THE NUMBER OF STATES SPECIFIED IN THE AB INITIO POTENTIAL", .{});

                return error.InvalidInput;
            }

            var adiabatic = try RealMatrix(T).initZero(ENERGY.rows, ENERGY.rows, allocator);
            defer adiabatic.deinit(allocator);

            for (0..ENERGY.rows) |j| adiabatic.ptr(j, j).* = ENERGY.at(j, 0);

            const bias_force = if (bias) |bs| try bs.force(adiabatic, state, i) else 0;

            return -GRADIENT.at(state * position.len + i, 0) + if (bias) |bs| switch (bs.variable) {
                .potential_energy => bias_force * GRADIENT.at(state * position.len + i, 0),
                .potential_energy_difference => |value| blk: {
                    const s1 = value.states[0];
                    const s2 = value.states[1];

                    break :blk bias_force * (GRADIENT.at(s2 * position.len + i, 0) - GRADIENT.at(s1 * position.len + i, 0));
                },
            } else 0;
        }

        pub fn runElectronicStructureCalculation(self: @This(), io: std.Io, system: ClassicalParticle(T), dir: std.Io.Dir, allocator: std.mem.Allocator) !void {
            try system.writeCoordinatesToXYZ(io, "molecule.xyz", dir);

            try appendFileToAnother(io, "molecule.xyz", "trajectory.xyz", dir);

            const path = try dir.realPathFileAlloc(io, ".", allocator);
            defer allocator.free(path);

            const states_str = try std.fmt.allocPrint(allocator, "{d}", .{self.states});
            defer allocator.free(states_str);

            const command = try std.mem.concat(allocator, u8, &.{ self.command, " -k ", states_str, " -s ", path, "/molecule.xyz" });
            defer allocator.free(command);

            const command_output = try executeCommand(io, command, allocator);
            defer allocator.free(command_output);
        }

        /// Number of states getter.
        pub fn nstate(self: @This()) usize {
            return self.states;
        }
    };
}

/// Append the contents of one file to another, creating the target file if it does not exist.
pub fn appendFileToAnother(io: std.Io, source_path: []const u8, target_path: []const u8, dir: std.Io.Dir) !void {
    dir.access(io, target_path, .{}) catch {
        const file = try dir.createFile(io, target_path, .{});
        defer file.close(io);
    };

    const source_file = try dir.openFile(io, source_path, .{ .mode = .read_only });
    defer source_file.close(io);

    const target_file = try dir.openFile(io, target_path, .{ .mode = .write_only });
    defer target_file.close(io);

    var write_buffer: [WRITE_BUFFER_SIZE]u8 = undefined;
    var writer = target_file.writer(io, &write_buffer);
    try writer.seekTo(try target_file.length(io));
    var writer_interface = &writer.interface;

    var buffer: [4096]u8 = undefined;
    const buffers = [_][]u8{ &buffer };

    while (true) {
        const bytes_read = source_file.readStreaming(io, &buffers) catch |err| {
            if (err == error.EndOfStream) {
                break;
            }
            else return err;
        };
        if (bytes_read == 0) break;
        try writer_interface.writeAll(buffer[0..bytes_read]);
    }
}
