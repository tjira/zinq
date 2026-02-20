//! File with ab initio potential struct and functions.

const std = @import("std");

const bias_potential = @import("bias_potential.zig");
const classical_particle = @import("classical_particle.zig");
const device_read = @import("device_read.zig");
const device_write = @import("device_write.zig");
const error_handling = @import("error_handling.zig");
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
const throw = error_handling.throw;

/// Struct holding parameters for the multidimensional ab initio potential.
pub fn AbInitioPotential(comptime T: type) type {
    return struct {
        command: []const u8,
        states: usize = 1,
        ncoord: usize,

        /// Evaluate the adiabatic potential energy matrix at given system state and time.
        pub fn evaluateAdiabatic(_: @This(), adiabatic_potential: *RealMatrix(T), dir: std.fs.Dir, allocator: std.mem.Allocator) !void {
            const dirname = try dir.realpathAlloc(allocator, "."); defer allocator.free(dirname);
            const path = try std.mem.concat(allocator, u8, &.{dirname, "/ENERGY.mat"}); defer allocator.free(path);

            const ENERGY = try readRealMatrix(T, path, allocator); defer ENERGY.deinit(allocator);

            adiabatic_potential.fill(0);

            for (0..ENERGY.rows) |i| adiabatic_potential.ptr(i, i).* = ENERGY.at(i, 0);
        }

        /// Comptime adiabatic force evaluation. The evaluateAdiabatic function needs to be run first.
        pub fn forceAdiabatic(_: @This(), i: usize, position: RealVector(T), _: T, state: usize, bias: ?BiasPotential(T), dir: std.fs.Dir, allocator: std.mem.Allocator) !T {
            const dirname = try dir.realpathAlloc(allocator, "."); defer allocator.free(dirname);
            const path_gradient = try std.mem.concat(allocator, u8, &.{dirname, "/GRADIENT.mat"}); defer allocator.free(path_gradient);
            const path_energy = try std.mem.concat(allocator, u8, &.{dirname, "/ENERGY.mat"}); defer allocator.free(path_energy);

            const ENERGY = try readRealMatrix(T, path_energy, allocator); defer ENERGY.deinit(allocator);
            const GRADIENT = try readRealMatrix(T, path_gradient, allocator); defer GRADIENT.deinit(allocator);

            var adiabatic = try RealMatrix(T).initZero(ENERGY.rows, ENERGY.rows, allocator); defer adiabatic.deinit(allocator);

            for (0..ENERGY.rows) |j| adiabatic.ptr(j, j).* = ENERGY.at(j, 0);

            return -GRADIENT.at(state * position.len + i, 0) + if (bias) |bs| bs.force(adiabatic, state, i, GRADIENT.at(state * position.len + i, 0)) else 0;
        }

        pub fn runElectronicStructureCalculation(self: @This(), system: ClassicalParticle(T), dir: std.fs.Dir, allocator: std.mem.Allocator) !void {
            try system.writeCoordinatesToXYZ("molecule.xyz", dir);

            const path = try dir.realpathAlloc(allocator, "."); defer allocator.free(path);
            const states_str = try std.fmt.allocPrint(allocator, "{d}", .{self.states}); defer allocator.free(states_str);
            const command = try std.mem.concat(allocator, u8, &.{self.command, " -k ", states_str, " -s ", path, "/molecule.xyz"}); defer allocator.free(command);

            const command_output = try executeCommand(command, allocator); defer allocator.free(command_output);
        }

        /// Number of dimensions getter.
        pub fn ndim(self: @This()) usize {
            return self.ncoord;
        }

        /// Number of states getter.
        pub fn nstate(self: @This()) usize {
            return self.states;
        }
    };
}
