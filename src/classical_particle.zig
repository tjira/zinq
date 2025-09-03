//! General system struct that holds all information about the classical like coordinates, velocities, etc.

const std = @import("std");

const electronic_potential = @import("electronic_potential.zig");
const global_variables = @import("global_variables.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const ElectronicPotential = electronic_potential.ElectronicPotential;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const A2AU = global_variables.A2AU;
const AN2M = global_variables.AN2M;
const FINITE_DIFFERENCES_STEP = global_variables.FINITE_DIFFERENCES_STEP;
const SM2AN = global_variables.SM2AN;

/// Classical particle struct that holds all information about the classical-like coordinates, velocities and masses.
pub fn ClassicalParticle(comptime T: type) type {
    return struct {
        atoms: ?[]usize,
        acceleration: RealVector(T),
        charge: i32,
        masses: RealVector(T),
        ndim: usize,
        position: RealVector(T),
        velocity: RealVector(T),

        allocator: std.mem.Allocator,

        /// Initialize the system using the number of dimensions, an array of masses, and an allocator. The positions and velocities are undefined after initialization.
        pub fn init(ndim: usize, masses: []const T, allocator: std.mem.Allocator) !@This() {
            const particle = @This(){
                .atoms = null,
                .acceleration = try RealVector(T).init(ndim, allocator),
                .charge = 0,
                .masses = try RealVector(T).init(ndim, allocator),
                .ndim = ndim,
                .position = try RealVector(T).init(ndim, allocator),
                .velocity = try RealVector(T).init(ndim, allocator),
                .allocator = allocator
            };

            @memcpy(particle.masses.data, masses);

            return particle;
        }

        /// Initialize the system using the number of dimensions, an array of masses, and an allocator. The positions and velocities are initialized to zero.
        pub fn initZero(ndim: usize, masses: []const T, allocator: std.mem.Allocator) !@This() {
            var particle = try @This().init(ndim, masses, allocator);

            particle.position.zero();
            particle.velocity.zero();
            particle.acceleration.zero();

            return particle;
        }

        /// Free the memory allocated for the system.
        pub fn deinit(self: @This()) void {
            if (self.atoms) |atoms| self.allocator.free(atoms);

            self.position.deinit();
            self.velocity.deinit();
            self.acceleration.deinit();
            self.masses.deinit();
        }

        /// Calculate the kinetic energy of the classical particle.
        pub fn kineticEnergy(self: @This()) T {
            var kinetic_energy: T = 0;

            for (0..self.ndim) |i| {
                kinetic_energy += 0.5 * self.masses.at(i) * self.velocity.at(i) * self.velocity.at(i);
            }

            return kinetic_energy;
        }

        /// Propagate the classical particle using velocity verlet algorithm.
        pub fn propagateVelocityVerlet(self: *@This(), potential: ElectronicPotential(T), potential_matrix: *RealMatrix(T), time: T, current_state: usize, time_step: T) !void {
            for (0..self.ndim) |i| {
                self.position.ptr(i).* += (self.velocity.at(i) + 0.5 * self.acceleration.at(i) * time_step) * time_step;
            }

            for (0..self.ndim) |i| {

                const force = try potential.forceAdiabatic(potential_matrix, self.position, time, current_state, i);

                const previous_acceleration = self.acceleration.at(i);

                self.acceleration.ptr(i).* = force / self.masses.at(i);

                self.velocity.ptr(i).* += 0.5 * (previous_acceleration + self.acceleration.at(i)) * time_step;
            }
        }

        /// Sets the position of the particle from normal distribution with given mean and standard deviation using the provided random number generator state.
        pub fn setPositionRandn(self: *@This(), mean: []const T, stdev: []const T, random: *std.Random) !void {
            if (mean.len != self.ndim or stdev.len != self.ndim) return error.DimensionMismatch;

            for (0..self.ndim) |i| {
                self.position.ptr(i).* = mean[i] + stdev[i] * random.floatNorm(T);
            }
        }

        /// Sets the momentum of the particle from normal distribution with given mean and standard deviation using the provided random number generator state.
        pub fn setMomentumRandn(self: *@This(), mean: []const T, stdev: []const T, random: *std.Random) !void {
            if (mean.len != self.ndim or stdev.len != self.ndim) return error.DimensionMismatch;

            for (0..self.ndim) |i| {
                self.velocity.ptr(i).* = (mean[i] + stdev[i] * random.floatNorm(T)) / self.masses.at(i);
            }
        }
    };
}

/// Read the system from a .xyz file.
pub fn read(comptime T: type, path: []const u8, charge: i32, allocator: std.mem.Allocator) !ClassicalParticle(T) {
    const file = try std.fs.cwd().openFile(path, .{}); defer file.close();

    var buffer: [1024]u8 = undefined; var reader = file.reader(&buffer); var reader_interface = &reader.interface;

    const natom = try std.fmt.parseInt(u32, try reader_interface.takeDelimiterExclusive('\n'), 10);

    var atoms = try allocator.alloc(usize, natom);

    _ = try reader_interface.discardDelimiterInclusive('\n');

    var position = try RealVector(T).init(3 * natom, allocator); defer position.deinit();

    for (0..natom) |i| {

        var it = std.mem.tokenizeAny(u8, try reader_interface.takeDelimiterExclusive('\n'), " "); 

        atoms[i] = SM2AN.get(it.next().?).?;

        for (0..3) |j| {
            position.ptr(3 * i + j).* = try std.fmt.parseFloat(T, it.next().?) * A2AU;
        }
    }

    var masses = try allocator.alloc(T, 3 * natom); defer allocator.free(masses);

    for (0..3 * natom) |i| masses[i] = AN2M[atoms[i / 3]];

    var system = try ClassicalParticle(T).initZero(3 * natom, masses, allocator);

    system.atoms = atoms; system.charge = charge; @memcpy(system.position.data, position.data);

    return system;
}
