//! General system struct that holds all information about the classical like coordinates, velocities, etc.

const std = @import("std");

const electronic_potential = @import("electronic_potential.zig");
const global_variables = @import("global_variables.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const ElectronicPotential = electronic_potential.ElectronicPotential;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const FINITE_DIFFERENCES_STEP = global_variables.FINITE_DIFFERENCES_STEP;

/// Classical particle struct that holds all information about the classical-like coordinates, velocities and masses.
pub fn ClassicalParticle(comptime T: type) type {
    return struct {
        ndim: usize,
        position: RealVector(T),
        velocity: RealVector(T),
        acceleration: RealVector(T),
        masses: RealVector(T),
        allocator: std.mem.Allocator,

        /// Initialize the system using the number of dimensions, an array of masses, and an allocator. The positions and velocities are undefined after initialization.
        pub fn init(ndim: usize, masses: []const T, allocator: std.mem.Allocator) !@This() {
            const particle = @This(){
                .ndim = ndim,
                .position = try RealVector(T).init(ndim, allocator),
                .velocity = try RealVector(T).init(ndim, allocator),
                .acceleration = try RealVector(T).init(ndim, allocator),
                .masses = try RealVector(T).init(ndim, allocator),
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
