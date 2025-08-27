//! General system struct that holds all information about the classical like coordinates, velocities, etc.

const std = @import("std");

const classical_dynamics = @import("classical_dynamics.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

/// Classical particle struct that holds all information about the classical-like coordinates, velocities and masses.
pub fn ClassicalParticle(comptime T: type) type {
    return struct {
        ndim: usize,
        position: RealVector(T),
        velocity: RealVector(T),
        masses: RealVector(T),
        allocator: std.mem.Allocator,

        /// Initialize the system using the number of dimensions, an array of masses, and an allocator. The positions and velocities are undefined after initialization.
        pub fn init(ndim: usize, masses: []const T, allocator: std.mem.Allocator) !ClassicalParticle(T) {
            const particle = ClassicalParticle(T){
                .ndim = ndim,
                .position = try RealVector(T).init(ndim, allocator),
                .velocity = try RealVector(T).init(ndim, allocator),
                .masses = try RealVector(T).init(ndim, allocator),
                .allocator = allocator
            };

            @memcpy(particle.masses.data, masses);

            return particle;
        }

        /// Initialize the system using the number of dimensions, an array of masses, and an allocator. The positions and velocities are initialized to zero.
        pub fn initZero(ndim: usize, masses: []const T, allocator: std.mem.Allocator) !ClassicalParticle(T) {
            var particle = try ClassicalParticle(T).init(ndim, masses, allocator);

            particle.position.zero();
            particle.velocity.zero();

            return particle;
        }

        /// Free the memory allocated for the system.
        pub fn deinit(self: ClassicalParticle(T)) void {
            self.position.deinit();
            self.velocity.deinit();
            self.masses.deinit();
        }

        /// Sets the position of the particle from normal distribution with given mean and standard deviation using the provided random number generator state.
        pub fn setPositionRandn(self: *ClassicalParticle(T), mean: []const T, stdev: []const T, seed: usize) !void {
            if (mean.len != self.ndim or stdev.len != self.ndim) return error.DimensionMismatch;

            var rng = std.Random.DefaultPrng.init(seed); var random = rng.random();

            for (0..self.ndim) |i| {
                self.position.ptr(i).* = mean[i] + stdev[i] * random.floatNorm(T);
            }
        }

        /// Sets the momentum of the particle from normal distribution with given mean and standard deviation using the provided random number generator state.
        pub fn setMomentumRandn(self: *ClassicalParticle(T), mean: []const T, stdev: []const T, seed: usize) !void {
            if (mean.len != self.ndim or stdev.len != self.ndim) return error.DimensionMismatch;

            var rng = std.Random.DefaultPrng.init(seed); var random = rng.random();

            for (0..self.ndim) |i| {
                self.velocity.ptr(i).* = (mean[i] + stdev[i] * random.floatNorm(T)) / self.masses.at(i);
            }
        }
    };
}
