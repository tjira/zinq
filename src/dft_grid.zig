//! DFTGrid module for generating integration grid points and weights for DFT calculations.

const std = @import("std");

const basis_set = @import("basis_set.zig");
const lebedev_quadrature_nodes = @import("lebedev_quadrature_nodes.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const BasisSet = basis_set.BasisSet;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const getLebedevGrid = lebedev_quadrature_nodes.getLebedevGrid;

/// Enumeration of available DFT functionals.
pub fn FunctionalGrid(comptime T: type) type {
    return union(enum) {
        uniform: Uniform(T),
        becke: Becke(T),
    };
}

/// Returns the Bragg-Slater radius (in Angstroms) for a given atomic number (Z).
pub fn getBraggSlaterRadius(comptime T: type, Z: usize) T {
    const r_f64: f64 = switch (Z) {
        1...2 => 0.35,
        3 => 1.45, 4 => 1.05, 5 => 0.85, 6 => 0.65, 7...8 => 0.60, 9...10 => 0.50,
        11 => 1.80, 12 => 1.50, 13 => 1.25, 14 => 1.15, 15...18 => 1.00,
        else => 1.50
    };
    return @as(T, @floatCast(r_f64));
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
        n_lebedev: usize = 302,

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

            var atomic_numbers = std.ArrayList(usize){}; defer atomic_numbers.deinit(allocator);

            const lebedev_points = try getLebedevGrid(self.n_lebedev);

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

                if (is_unique) {
                    try centers.append(allocator, cg.center); try atomic_numbers.append(allocator, cg.atomic_number);
                }
            }

            const total_points = centers.items.len * self.n_radial * self.n_lebedev;

            var points = try RealMatrix(T).init(total_points, 3, allocator); errdefer points.deinit(allocator);
            var weights = try RealVector(T).init(total_points, allocator); errdefer weights.deinit(allocator);

            var W_vals = try allocator.alloc(T, centers.items.len); defer allocator.free(W_vals);

            var R_AB = try allocator.alloc(T, centers.items.len * centers.items.len); defer allocator.free(R_AB);
            var a_AB = try allocator.alloc(T, centers.items.len * centers.items.len); defer allocator.free(a_AB);

            for (0..centers.items.len) |A| for (0..centers.items.len) |B| {

                const dx = centers.items[A][0] - centers.items[B][0];
                const dy = centers.items[A][1] - centers.items[B][1];
                const dz = centers.items[A][2] - centers.items[B][2];

                R_AB[A * centers.items.len + B] = std.math.sqrt(dx * dx + dy * dy + dz * dz);

                if (A == B) a_AB[A * centers.items.len + B] = 0;

                if (A != B) {

                    const radius_A = getBraggSlaterRadius(T, atomic_numbers.items[A]);
                    const radius_B = getBraggSlaterRadius(T, atomic_numbers.items[B]);

                    const u = (radius_A - radius_B) / (radius_A + radius_B); var a = u / (u * u - 1);
                    
                    if (a > 0.5) a = 0.5; if (a < -0.5) a = -0.5;
                    
                    a_AB[A * centers.items.len + B] = a;
                }
            };

            var idx: usize = 0;

            for (0..centers.items.len) |A| for (0..self.n_radial) |ir| {

                const x = std.math.cos(std.math.pi * (@as(T, @floatFromInt(ir)) + 0.5) / @as(T, @floatFromInt(self.n_radial)));
                const w_cheb = std.math.pi / @as(T, @floatFromInt(self.n_radial)) * std.math.sqrt(1 - x*x);

                const r = getBraggSlaterRadius(T, atomic_numbers.items[A]) * (1 + x) / (1 - x);
                const w_r = 2 * w_cheb * getBraggSlaterRadius(T, atomic_numbers.items[A]) / ((1 - x) * (1 - x)) * r * r;

                for (lebedev_points.nodes, lebedev_points.weights) |lebedev_node, lebedev_weight| {

                    points.ptr(idx, 0).* = centers.items[A][0] + r * lebedev_node[0];
                    points.ptr(idx, 1).* = centers.items[A][1] + r * lebedev_node[1];
                    points.ptr(idx, 2).* = centers.items[A][2] + r * lebedev_node[2];

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

                            W_I *= beckeStep(mu + a_AB[I * centers.items.len + J] * (1 - mu * mu));
                        }

                        W_vals[I] = W_I; P_sum += W_I;
                    }

                    weights.ptr(idx).* = 4 * std.math.pi * w_r * lebedev_weight * W_vals[A] / P_sum; idx += 1;
                }
            };

            return .{points, weights};
        }
    };
}
