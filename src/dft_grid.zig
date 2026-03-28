//! DFTGrid module for generating integration grid points and weights for DFT calculations.

const std = @import("std");

const basis_set = @import("basis_set.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const BasisSet = basis_set.BasisSet;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

/// Enumeration of available DFT functionals.
pub fn FunctionalGrid(comptime T: type) type {
    return union(enum) {
        uniform: Uniform(T),
        becke: Becke(T),
    };
}

/// Generate the integration grid points and weights for DFT calculations based on the provided basis set.
pub fn getGrid(comptime T: type, grid: FunctionalGrid(T), basis: BasisSet(T), allocator: std.mem.Allocator) !struct {RealMatrix(T), RealVector(T)} {
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
        pub fn get(self: @This(), _: BasisSet(T), allocator: std.mem.Allocator) !struct {RealMatrix(T), RealVector(T)} {
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
            
            return .{points, weights};
        }
    };
}

/// This is the standard grid topology used in almost all production DFT codes.
pub fn Becke(comptime T: type) type {
    return struct {
        n_radial: usize = 50,
        n_theta: usize = 30,
        n_phi: usize = 60,
        atomic_radius: T = 1.0,

        // Becke's 3rd order smoothing polynomial
        pub fn beckeStep(mu: T) T {
            var p = 1.5 * mu - 0.5 * mu * mu * mu;

            p = 1.5 * p - 0.5 * p * p * p;
            p = 1.5 * p - 0.5 * p * p * p;

            return 0.5 * (1.0 - p);
        }

        /// Generate the integration grid points and weights for DFT calculations based on the provided basis set.
        pub fn get(self: @This(), basis: BasisSet(T), allocator: std.mem.Allocator) !struct {RealMatrix(T), RealVector(T)} {
            var centers = std.ArrayList([3]T){}; defer centers.deinit(allocator);

            for (basis.contracted_gaussians) |cg| {

                var is_unique = true;

                for (centers.items) |uc| {

                    const dx = cg.center[0] - uc[0];
                    const dy = cg.center[1] - uc[1];
                    const dz = cg.center[2] - uc[2];

                    if (dx * dx + dy * dy + dz * dz < 1e-8) {
                        is_unique = false; break;
                    }
                }

                if (is_unique) try centers.append(allocator, cg.center);
            }

            const total_points = centers.items.len * self.n_radial * self.n_theta * self.n_phi;

            var points = try RealMatrix(T).init(total_points, 3, allocator); errdefer points.deinit(allocator);
            var weights = try RealVector(T).init(total_points, allocator); errdefer weights.deinit(allocator);

            var W_vals = try allocator.alloc(T, centers.items.len); defer allocator.free(W_vals);
            var R_AB = try allocator.alloc(T,  centers.items.len * centers.items.len); defer allocator.free(R_AB);

            for (0..centers.items.len) |A| for (0..centers.items.len) |B| {

                const dx = centers.items[A][0] - centers.items[B][0];
                const dy = centers.items[A][1] - centers.items[B][1];
                const dz = centers.items[A][2] - centers.items[B][2];

                R_AB[A * centers.items.len + B] = std.math.sqrt(dx * dx + dy * dy + dz * dz);
            };

            var idx: usize = 0;

            for (0..centers.items.len) |A| for (0..self.n_radial) |ir| {

                const x = std.math.cos(std.math.pi * (@as(T, @floatFromInt(ir)) + 0.5) / @as(T, @floatFromInt(self.n_radial)));
                const w_cheb = std.math.pi / @as(T, @floatFromInt(self.n_radial)) * std.math.sqrt(1 - x*x);

                const r = self.atomic_radius * (1 + x) / (1 - x);
                const w_r = 2 * w_cheb * self.atomic_radius / ((1 - x) * (1 - x)) * r * r;

                for (0..self.n_theta) |it| {

                    const theta = std.math.pi * (@as(T, @floatFromInt(it)) + 0.5) / @as(T, @floatFromInt(self.n_theta));
                    const w_theta = std.math.pi / @as(T, @floatFromInt(self.n_theta)) * @sin(theta);

                    for (0..self.n_phi) |ip| {

                        const phi = 2 * std.math.pi * @as(T, @floatFromInt(ip)) / @as(T, @floatFromInt(self.n_phi));
                        const w_phi = 2 * std.math.pi / @as(T, @floatFromInt(self.n_phi));

                        points.ptr(idx, 0).* = centers.items[A][0] + r * std.math.sin(theta) * std.math.cos(phi);
                        points.ptr(idx, 1).* = centers.items[A][1] + r * std.math.sin(theta) * std.math.sin(phi);
                        points.ptr(idx, 2).* = centers.items[A][2] + r * std.math.cos(theta);

                        var P_sum: T = 0;

                        for (0..centers.items.len) |I| {

                            var W_I: T = 1;

                            const r_I = std.math.sqrt(
                                std.math.pow(T, points.at(idx, 0) - centers.items[I][0], 2) + 
                                std.math.pow(T, points.at(idx, 1) - centers.items[I][1], 2) + 
                                std.math.pow(T, points.at(idx, 2) - centers.items[I][2], 2)
                            );

                            for (0..centers.items.len) |J| {

                                if (I == J) continue;

                                const r_J = std.math.sqrt(
                                    std.math.pow(T, points.at(idx, 0) - centers.items[J][0], 2) + 
                                    std.math.pow(T, points.at(idx, 1) - centers.items[J][1], 2) + 
                                    std.math.pow(T, points.at(idx, 2) - centers.items[J][2], 2)
                                );

                                const mu = (r_I - r_J) / R_AB[I * centers.items.len + J];

                                W_I *= beckeStep(mu);
                            }

                            W_vals[I] = W_I; P_sum += W_I;
                        }

                        weights.ptr(idx).* = w_r * w_theta * w_phi * W_vals[A] / P_sum; idx += 1;
                    }
                }
            };

            return .{points, weights};
        }
    };
}
