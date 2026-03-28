//! DFTGrid module for generating integration grid points and weights for DFT calculations.

const std = @import("std");

const basis_set = @import("basis_set.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const BasisSet = basis_set.BasisSet;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

/// Enumeration of available DFT functionals.
pub fn DFTGrid(comptime T: type) type {
    return union(enum) {
        uniform: Uniform(T)
    };
}

/// Return type of the DFT functional.
pub fn DFTGridOutput(comptime T: type) type {
    return struct {points: RealMatrix(T), weights: RealVector(T)};
}

/// Generate the integration grid points and weights for DFT calculations based on the provided basis set.
pub fn getGrid(comptime T: type, grid: DFTGrid(T), basis: BasisSet(T), allocator: std.mem.Allocator) !DFTGridOutput(T) {
    switch (grid) {
        inline else => |g| return g.get(basis, allocator)
    }
}

/// Uniform grid for DFT integration, which generates a regular grid of points and corresponding weights within specified limits.
pub fn Uniform(comptime T: type) type {
    return struct {
        limits: [2]T = .{-6, 6},
        points: u32 = 64,

        /// Generate the integration grid points and weights for DFT calculations based on the provided basis set.
        pub fn get(self: @This(), _: BasisSet(T), allocator: std.mem.Allocator) !DFTGridOutput(T) {
            var points = try RealMatrix(T).init(self.points * self.points * self.points, 3, allocator); errdefer points.deinit(allocator);
            var weights = try RealVector(T).init(self.points * self.points * self.points, allocator); errdefer weights.deinit(allocator);

            const dx = (self.limits[1] - self.limits[0]) / @as(T, @floatFromInt(self.points - 1));

            for (0..self.points) |i| {
                for (0..self.points) |j| {
                    for (0..self.points) |k| {

                        const idx = i * self.points * self.points + j * self.points + k;

                        points.ptr(idx, 0).* = self.limits[0] + (self.limits[1] - self.limits[0]) * @as(T, @floatFromInt(i)) / @as(T, @floatFromInt(self.points - 1));
                        points.ptr(idx, 1).* = self.limits[0] + (self.limits[1] - self.limits[0]) * @as(T, @floatFromInt(j)) / @as(T, @floatFromInt(self.points - 1));
                        points.ptr(idx, 2).* = self.limits[0] + (self.limits[1] - self.limits[0]) * @as(T, @floatFromInt(k)) / @as(T, @floatFromInt(self.points - 1));

                        weights.ptr(idx).* = dx * dx * dx;
                    }
                }
            }
            
            return .{.points = points, .weights = weights};
        }
    };
}

