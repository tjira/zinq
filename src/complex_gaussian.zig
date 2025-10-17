//! Complex gaussian function implementation in Zig.

const std = @import("std");

const array_functions = @import("array_functions.zig");
const complex_matrix = @import("complex_matrix.zig");
const complex_vector = @import("complex_vector.zig");
const electronic_potential = @import("electronic_potential.zig");
const error_handling = @import("error_handling.zig");
const grid_generator = @import("grid_generator.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const Complex = std.math.Complex;
const ComplexMatrix = complex_matrix.ComplexMatrix;
const ComplexVector = complex_vector.ComplexVector;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const positionAtRow = grid_generator.positionAtRow;
const prod = array_functions.prod;
const throw = error_handling.throw;

/// Complex Gaussian function implementation in Zig.
pub fn ComplexGaussian(comptime T: type) type {
    return struct {
        center: []T,
        gamma: []T,
        momentum: []T,
        norm: T,

        /// Initialize a new complex Gaussian with given parameters.
        pub fn init(center: []const T, gamma: []const T, momentum: []const T, allocator: std.mem.Allocator) !@This() {
            if (center.len != gamma.len or center.len != momentum.len or gamma.len != momentum.len) {
                return throw(ComplexGaussian(T), "CENTER, GAMMA, AND MOMENTUM MUST HAVE THE SAME DIMENSIONALITY", .{});
            }

            var cg = @This(){
                .center = try allocator.alloc(T, center.len),
                .gamma = try allocator.alloc(T, center.len),
                .momentum = try allocator.alloc(T, center.len),
                .norm = std.math.pow(T, std.math.pi, @as(T, @floatFromInt(center.len)) / 4) / std.math.pow(T, prod(T, gamma), 0.25),
            };

            for (center, 0..) |c, i| cg.center[i] = c;
            for (gamma, 0..) |g, i| cg.gamma[i] = g;
            for (momentum, 0..) |p, i| cg.momentum[i] = p;

            return cg;
        }

        /// Free allocated memory for the complex Gaussian.
        pub fn deinit(self: @This(), allocator: std.mem.Allocator) void {
            allocator.free(self.center);
            allocator.free(self.gamma);
            allocator.free(self.momentum);
        }

        /// Evaluate the complex Gaussian function at a given point x.
        pub fn evaluate(self: @This(), q: []const T) Complex(T) {
            var result = Complex(T).init(1, 0);

            for (self.center, 0..) |c, i| {

                const dq = q[i] - c;

                const exponent = Complex(T).init(-0.5 * self.gamma[i] * dq * dq, self.momentum[i] * dq);

                result = result.mul(std.math.complex.exp(exponent));
            }

            return result.div(Complex(T).init(self.norm, 0));
        }

        /// Calculate the kinetic matrix element between this complex Gaussian and another.
        pub fn kinetic(self: @This(), other: @This(), mass: []const T) !Complex(T) {
            if (self.center.len != other.center.len or self.center.len != mass.len or other.center.len != mass.len) {
                return throw(Complex(T), "BOTH COMPLEX GAUSSIANS AND MASS ARRAY MUST HAVE THE SAME DIMENSIONALITY", .{});
            }

            var result = Complex(T).init(1, 0);

            for (self.center, other.center, self.gamma, other.gamma, self.momentum, other.momentum, mass) |q1, q2, g1, g2, p1, p2, m| {

                const dq = q1 - q2;
                const dp = p1 - p2;

                const factor = std.math.sqrt(0.5 * std.math.pi) / std.math.pow(T, g1 + g2, 2.5) / m;

                const factorRe = (std.math.pow(T, g2 * p1 + g1 * p2, 2) + g1 * g2 * (g1 + g2 - g1 * g2 * dq * dq));
                const factorIm = 2 * g1 * g2 * dq * (g2 * p1 + g1 * p2);

                const expRe = -0.5 * (dp * dp + g1 * g2 * dq * dq) / (g1 + g2);
                const expIm = dq * (g2 * p1 + g1 * p2) / (g1 + g2);

                const factor_full = Complex(T).init(factor * factorRe, factor * factorIm);
                const exponent_full = Complex(T).init(expRe, expIm);

                result = result.mul(std.math.complex.exp(exponent_full).mul(factor_full));
            }

            return result.div(Complex(T).init(self.norm * other.norm, 0));
        }

        /// Compute the overlap integral between this complex Gaussian and another.
        pub fn overlap(self: @This(), other: @This()) !Complex(T) {
            if (self.center.len != other.center.len) return throw(Complex(T), "BOTH COMPLEX GAUSSIANS MUST HAVE THE SAME DIMENSIONALITY", .{});

            var result = Complex(T).init(1, 0);

            for (self.center, other.center, self.gamma, other.gamma, self.momentum, other.momentum) |q1, q2, g1, g2, p1, p2| {

                const dq = q1 - q2;
                const dp = p1 - p2;

                const factor = std.math.sqrt(2 * std.math.pi / (g1 + g2));

                const numRe = dp * dp + g1 * g2 * dq * dq;
                const numIm = -2 * dq * (g2 * p1 + g1 * p2);
                const denRe = -2 * (g1 + g2);

                const exponent = std.math.Complex(T).init(numRe / denRe, numIm / denRe);

                result = result.mul(std.math.complex.exp(exponent).mul(Complex(T).init(factor, 0)));
            }

            return result.div(Complex(T).init(self.norm * other.norm, 0));
        }

        /// Compute the potential energy matrix element between this complex Gaussian and another for a given potential function.
        pub fn potential(self: @This(), other: @This(), pot: ElectronicPotential(T), limits: []const []const T, npoint: usize, time: T, allocator: std.mem.Allocator) !ComplexMatrix(T) {
            var V = try ComplexMatrix(T).initZero(pot.nstate(), pot.nstate(), allocator);
            var v = try    RealMatrix(T).initZero(pot.nstate(), pot.nstate(), allocator);

            var q = try RealVector(T).init(limits.len, allocator); defer q.deinit(); var dq: T = 1;

            for (0..limits.len) |i| {
                dq *= (limits[i][1] - limits[i][0]) / @as(T, @floatFromInt(npoint - 1));
            }

            for (0..std.math.pow(usize, npoint, limits.len)) |k| {

                positionAtRow(T, &q, k, limits.len, npoint, limits);

                pot.evaluateDiabatic(&v, q, time);

                const g1_value = self .evaluate(q.data);
                const g2_value = other.evaluate(q.data);

                for (0..V.rows) |i| for (0..V.cols) |j| {
                    V.ptr(i, j).* = V.at(i, j).add(g1_value.conjugate().mul(Complex(T).init(v.at(i, j), 0)).mul(g2_value));
                };
            }

            V.muls(Complex(T).init(dq, 0));

            v.deinit(); return V;
        }
    };
}
