//! Complex gaussian function implementation in Zig.

const std = @import("std");

const array_functions = @import("array_functions.zig");
const complex_matrix = @import("complex_matrix.zig");
const complex_vector = @import("complex_vector.zig");
const device_write = @import("device_write.zig");
const electronic_potential = @import("electronic_potential.zig");
const error_handling = @import("error_handling.zig");
const grid_generator = @import("grid_generator.zig");
const hermite_quadrature_nodes = @import("hermite_quadrature_nodes.zig");
const math_functions = @import("math_functions.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const object_array = @import("object_array.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const Complex = std.math.Complex;
const ComplexMatrix = complex_matrix.ComplexMatrix;
const ComplexMatrixArray = object_array.ComplexMatrixArray;
const ComplexVector = complex_vector.ComplexVector;
const ComplexVectorArray = object_array.ComplexVectorArray;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const mm = matrix_multiplication.mm;
const positionAtRow = grid_generator.positionAtRow;
const powi = math_functions.powi;
const printComplexMatrix = device_write.printComplexMatrix;
const prod = array_functions.prod;
const throw = error_handling.throw;

const HERMITE_NODES = hermite_quadrature_nodes.HERMITE_NODES;
const HERMITE_WEIGHTS = hermite_quadrature_nodes.HERMITE_WEIGHTS;

/// Complex Gaussian function implementation in Zig.
pub fn ComplexGaussian(comptime T: type) type {
    return struct {
        position: []T,
        gamma: []T,
        momentum: []T,
        norm: T,
        allocator: std.mem.Allocator,

        /// Initialize a new complex Gaussian with given parameters.
        pub fn init(position: []const T, gamma: []const T, momentum: []const T, allocator: std.mem.Allocator) !@This() {
            if (position.len != gamma.len or position.len != momentum.len or gamma.len != momentum.len) {
                return throw(ComplexGaussian(T), "POSITION, GAMMA, AND MOMENTUM MUST HAVE THE SAME DIMENSIONALITY", .{});
            }

            var cg = @This(){
                .position = try allocator.alloc(T, position.len),
                .gamma = try allocator.alloc(T, position.len),
                .momentum = try allocator.alloc(T, position.len),
                .norm = std.math.pow(T, std.math.pi, @as(T, @floatFromInt(position.len)) / 4) / std.math.pow(T, prod(T, gamma), 0.25),
                .allocator = allocator
            };

            for (position, 0..) |c, i| cg.position[i] = c;
            for (gamma, 0..) |g, i| cg.gamma[i] = g;
            for (momentum, 0..) |p, i| cg.momentum[i] = p;

            return cg;
        }

        /// Free allocated memory for the complex Gaussian.
        pub fn deinit(self: @This()) void {
            self.allocator.free(self.position);
            self.allocator.free(self.gamma);
            self.allocator.free(self.momentum);
        }

        /// Clone the complex Gaussian.
        pub fn clone(self: @This()) !@This() {
            return try @This().init(self.position, self.gamma, self.momentum, self.allocator);
        }

        /// Returns the derivative of the electronic coefficients where this gaussian is shared between them.
        pub fn coefficientDerivative(self: @This(), coefs: ComplexVector(T), pot: ElectronicPotential(T), n_nodes: usize, time: T) !ComplexVector(T) {
            var dc = try ComplexVector(T).initZero(coefs.len, self.allocator);

            const V = try self.potential(self, pot, n_nodes, time); defer V.deinit();

            var rho = try ComplexMatrix(T).initZero(V.rows, V.cols, self.allocator); defer rho.deinit();
            var VC = try ComplexVector(T).init(V.cols, self.allocator); defer VC.deinit();

            for (0..coefs.len) |i| for (0..coefs.len) |j| {
                rho.ptr(i, j).* = coefs.at(i).mul(coefs.at(j).conjugate());
            };

            for (0..rho.cols) |i| for (0..rho.rows) |j| {
                rho.ptr(i, j).* = if (i == j) Complex(T).init(1, 0).sub(rho.at(i, i)) else rho.at(i, j).neg();
            };

            var VC_matrix = VC.asMatrix(); try mm(T, &VC_matrix, V, false, coefs.asMatrix(), false);

            for (0..dc.len) |i| for (0..coefs.len) |j| {
                dc.ptr(i).* = dc.at(i).add(rho.at(i, j).mul(VC.at(j)).mulbyi().neg());
            };

            return dc;
        }

        /// Evaluate the complex Gaussian function at a given point x.
        pub fn evaluate(self: @This(), q: []const T) Complex(T) {
            var result = Complex(T).init(1, 0);

            for (self.position, 0..) |c, i| {

                const dq = q[i] - c;

                const exponent = Complex(T).init(-0.5 * self.gamma[i] * dq * dq, self.momentum[i] * dq);

                result = result.mul(std.math.complex.exp(exponent));
            }

            return result.div(Complex(T).init(self.norm, 0));
        }

        /// Calculate the kinetic matrix element between this complex Gaussian and another.
        pub fn kinetic(self: @This(), other: @This(), mass: []const T) !Complex(T) {
            if (self.position.len != other.position.len or self.position.len != mass.len or other.position.len != mass.len) {
                return throw(Complex(T), "BOTH POSITION AND MASS ARRAY MUST HAVE THE SAME DIMENSIONALITY", .{});
            }

            var result = Complex(T).init(0, 0);

            for (self.position, other.position, self.gamma, other.gamma, self.momentum, other.momentum, mass) |q1, q2, g1, g2, p1, p2, m| {

                const dq = q1 - q2;
                const dp = p1 - p2;

                const factor = std.math.sqrt(0.25 * (g1 + g2)) / std.math.pow(T, g1 + g2, 2.5) / m;

                const factorRe = (std.math.pow(T, g2 * p1 + g1 * p2, 2) + g1 * g2 * (g1 + g2 - g1 * g2 * dq * dq));
                const factorIm = 2 * g1 * g2 * dq * (g2 * p1 + g1 * p2);

                const expRe = -0.5 * (dp * dp + g1 * g2 * dq * dq) / (g1 + g2);
                const expIm = dq * (g2 * p1 + g1 * p2) / (g1 + g2);

                const factor_full = Complex(T).init(factor * factorRe, factor * factorIm);
                const exponent_full = Complex(T).init(expRe, expIm);

                result = result.add(std.math.complex.exp(exponent_full).mul(factor_full));
            }

            return result.mul(try self.overlap(other));
        }

        /// Returns the kinetic energy expectation value of this complex Gaussian.
        pub fn kineticEnergy(self: @This(), mass: []const T) !T {
            if (self.position.len != mass.len) return throw(T, "POSITION AND MASS ARRAY MUST HAVE THE SAME DIMENSIONALITY IN COMPLEX GAUSSIAN", .{});

            return (try kinetic(self, self, mass)).re;
        }

        /// Calculates the derivative of the momentum.
        pub fn momentumDerivative(self: @This(), pot: ElectronicPotential(T), coefs: ComplexVector(T), n_nodes: usize, time: T) !RealVector(T) {
            const dV = try ComplexMatrixArray(T).init(self.position.len, .{.rows = pot.nstate(), .cols = pot.nstate()}, self.allocator); defer dV.deinit();

            for (0..dV.len) |i| {
                dV.ptr(i).deinit(); dV.ptr(i).* = try self.potentialDerivative(self, pot, i, n_nodes, time);
            }

            var F = try RealVector(T).initZero(dV.len, self.allocator);

            for (0..F.len) |i| for (0..coefs.len) |j| for (0..coefs.len) |k| {
                F.ptr(i).* -= coefs.at(j).conjugate().mul(dV.at(i).at(j, k)).mul(coefs.at(k)).re;
            };

            return F;
        }

        /// Compute the overlap integral between this complex Gaussian and another.
        pub fn overlap(self: @This(), other: @This()) !Complex(T) {
            if (self.position.len != other.position.len) return throw(Complex(T), "BOTH COMPLEX GAUSSIANS MUST HAVE THE SAME DIMENSIONALITY", .{});

            var result = Complex(T).init(1, 0);

            for (self.position, other.position, self.gamma, other.gamma, self.momentum, other.momentum) |q1, q2, g1, g2, p1, p2| {

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

        /// Calculates the derivative of the position.
        pub fn positionDerivative(self: @This(), mass: []const T) !RealVector(T) {
            var dq = try RealVector(T).init(self.position.len, self.allocator);

            for (0..dq.len) |i| {
                dq.ptr(i).* = self.momentum[i] / mass[i];
            }

            return dq;
        }

        /// Compute the potential energy matrix element between this complex Gaussian and another for a given potential function.
        pub fn potential(self: @This(), other: @This(), pot: ElectronicPotential(T), n_nodes: usize, time: T) !ComplexMatrix(T) {
            var V = try ComplexMatrix(T).initZero(pot.nstate(), pot.nstate(), self.allocator);

            var variables = try self.allocator.alloc(T, pot.ndim()); defer self.allocator.free(variables);
            var sigmas = try self.allocator.alloc(T, pot.ndim()); defer self.allocator.free(sigmas);
            var mus = try self.allocator.alloc(T, pot.ndim()); defer self.allocator.free(mus);

            for (0..pot.ndim()) |i| {

                const s1 = 1 / std.math.sqrt(self.gamma[i]); const s2 = 1 / std.math.sqrt(other.gamma[i]);

                mus[i] = (self.position[i] * s2 * s2 + other.position[i] * s1 * s1) / (s1 * s1 + s2 * s2);

                sigmas[i] = (s1 * s2) / std.math.sqrt(s1 * s1 + s2 * s2);
            }

            for (0..std.math.pow(usize, n_nodes, pot.ndim())) |i| {

                var temp = i; var weight: T = 1; var complex_exponent = Complex(T).init(0, 0);

                for (0..pot.ndim()) |j| {

                    variables[j] = std.math.sqrt2 * sigmas[j] * HERMITE_NODES[n_nodes][temp % n_nodes] + mus[j]; 

                    complex_exponent = complex_exponent.add(Complex(T).init(0, self.momentum[j] * (self.position[j] - variables[j]) + other.momentum[j] * (variables[j] - other.position[j])));

                    weight *= HERMITE_WEIGHTS[n_nodes][temp % n_nodes]; temp /= n_nodes;
                }

                for (0..V.rows) |j| for (j..V.cols) |k| {

                    const value = try pot.evaluateDiabaticElement(j, k, RealVector(T){.data = variables, .len = variables.len, .allocator = null}, time);

                    V.ptr(j, k).* = V.at(j, k).add(Complex(T).init(weight * value, 0).mul(std.math.complex.exp(complex_exponent))); V.ptr(k, j).* = V.at(j, k);
                };
            }

            for (0..pot.ndim()) |i| {
                V.muls(Complex(T).init(std.math.sqrt2 * sigmas[i], 0));
            }

            V.divs(Complex(T).init(self.norm * other.norm, 0));

            return V;
        }

        /// Compute the potential energy derivative matrix element between this complex Gaussian and another for a given potential function and a specified coordinate.
        pub fn potentialDerivative(self: @This(), other: @This(), pot: ElectronicPotential(T), index: usize, n_nodes: usize, time: T) !ComplexMatrix(T) {
            var V = try ComplexMatrix(T).initZero(pot.nstate(), pot.nstate(), self.allocator);

            var variables = try self.allocator.alloc(T, pot.ndim()); defer self.allocator.free(variables);
            var sigmas = try self.allocator.alloc(T, pot.ndim()); defer self.allocator.free(sigmas);
            var mus = try self.allocator.alloc(T, pot.ndim()); defer self.allocator.free(mus);

            for (0..pot.ndim()) |i| {

                const s1 = 1 / std.math.sqrt(self.gamma[i]); const s2 = 1 / std.math.sqrt(other.gamma[i]);

                mus[i] = (self.position[i] * s2 * s2 + other.position[i] * s1 * s1) / (s1 * s1 + s2 * s2);

                sigmas[i] = (s1 * s2) / std.math.sqrt(s1 * s1 + s2 * s2);
            }

            for (0..std.math.pow(usize, n_nodes, pot.ndim())) |i| {

                var temp = i; var weight: T = 1;  var complex_exponent = Complex(T).init(0, 0);

                for (0..pot.ndim()) |j| {

                    variables[j] = std.math.sqrt2 * sigmas[j] * HERMITE_NODES[n_nodes][temp % n_nodes] + mus[j]; 

                    complex_exponent = complex_exponent.add(Complex(T).init(0, self.momentum[j] * (self.position[j] - variables[j]) + other.momentum[j] * (variables[j] - other.position[j])));

                    weight *= HERMITE_WEIGHTS[n_nodes][temp % n_nodes]; temp /= n_nodes;
                }

                for (0..V.rows) |j| for (j..V.cols) |k| {

                    const value = try pot.evaluateDiabaticElementDerivative1(j, k, RealVector(T){.data = variables, .len = variables.len, .allocator = null}, time, index, 1e-8);

                    V.ptr(j, k).* = V.at(j, k).add(Complex(T).init(weight * value, 0).mul(std.math.complex.exp(complex_exponent))); V.ptr(k, j).* = V.at(j, k);
                };
            }

            for (0..pot.ndim()) |i| {
                V.muls(Complex(T).init(std.math.sqrt2 * sigmas[i], 0));
            }

            V.divs(Complex(T).init(self.norm * other.norm, 0));

            return V;
        }

        /// Returns the potential energy expectation value of this complex Gaussian for a given potential operator.
        pub fn potentialEnergy(self: @This(), pot: ElectronicPotential(T), coefs: ComplexVector(T), n_nodes: usize, time: T) !T {
            const V = try self.potential(self, pot, n_nodes, time); defer V.deinit();

            var potential_energy: Complex(T) = Complex(T).init(0, 0);

            for (0..coefs.len) |i| for (0..coefs.len) |j| {
                potential_energy = potential_energy.add(coefs.at(i).conjugate().mul(V.at(i, j)).mul(coefs.at(j)));
            };

            return potential_energy.re;
        }
    };
}
