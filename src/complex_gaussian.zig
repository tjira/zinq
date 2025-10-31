//! Complex gaussian function implementation in Zig.

const std = @import("std");

const array_functions = @import("array_functions.zig");
const complex_matrix = @import("complex_matrix.zig");
const complex_vector = @import("complex_vector.zig");
const electronic_potential = @import("electronic_potential.zig");
const error_handling = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");
const grid_generator = @import("grid_generator.zig");
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
const prod = array_functions.prod;
const throw = error_handling.throw;

const FINITE_DIFFERENCES_STEP = global_variables.FINITE_DIFFERENCES_STEP;

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
        pub fn coefficientDerivative(self: @This(), coefs: ComplexVector(T), pot: ElectronicPotential(T), limits: []const []const T, npoint: usize, time: T) !ComplexVector(T) {
            var dc = try ComplexVector(T).initZero(coefs.len, self.allocator);

            const V = try self.potential(self, pot, limits, npoint, time); defer V.deinit();

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
                return throw(Complex(T), "BOTH COMPLEX GAUSSIANS AND MASS ARRAY MUST HAVE THE SAME DIMENSIONALITY", .{});
            }

            var result = Complex(T).init(1, 0);

            for (self.position, other.position, self.gamma, other.gamma, self.momentum, other.momentum, mass) |q1, q2, g1, g2, p1, p2, m| {

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

        /// Returns the kinetic energy expectation value of this complex Gaussian.
        pub fn kineticEnergy(self: @This(), mass: []const T) !T {
            if (self.position.len != mass.len) return throw(T, "COMPLEX GAUSSIAN AND MASS ARRAY MUST HAVE THE SAME DIMENSIONALITY", .{});

            return (try kinetic(self, self, mass)).re;
        }

        /// Calculates the derivative of the momentum.
        pub fn momentumDerivative(self: @This(), pot: ElectronicPotential(T), coefs: ComplexVector(T), limits: []const []const T, npoint: usize, time: T) !RealVector(T) {
            const dV = try ComplexMatrixArray(T).init(self.position.len, .{.rows = pot.nstate(), .cols = pot.nstate()}, self.allocator); defer dV.deinit();

            for (0..dV.len) |i| {
                dV.ptr(i).deinit(); dV.ptr(i).* = try self.potentialDerivative(self, pot, i, limits, npoint, time);
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
        pub fn potential(self: @This(), other: @This(), pot: ElectronicPotential(T), limits: []const []const T, npoint: usize, time: T) !ComplexMatrix(T) {
            var V = try ComplexMatrix(T).initZero(pot.nstate(), pot.nstate(), self.allocator);
            var v = try    RealMatrix(T).initZero(pot.nstate(), pot.nstate(), self.allocator);

            var q = try RealVector(T).init(limits.len, self.allocator); defer q.deinit(); var dq: T = 1;

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

        /// Compute the potential energy derivative matrix element between this complex Gaussian and another for a given potential function and a specified coordinate.
        pub fn potentialDerivative(self: @This(), other: @This(), pot: ElectronicPotential(T), index: usize, limits: []const []const T, npoint: usize, time: T) !ComplexMatrix(T) {
            var V = try ComplexMatrix(T).initZero(pot.nstate(), pot.nstate(), self.allocator);

            var vDerivative = try RealMatrix(T).initZero(pot.nstate(), pot.nstate(), self.allocator); defer vDerivative.deinit();
            var vPlus  = try RealMatrix(T).initZero(pot.nstate(), pot.nstate(), self.allocator); defer vPlus.deinit();
            var vMinus = try RealMatrix(T).initZero(pot.nstate(), pot.nstate(), self.allocator); defer vMinus.deinit();

            var q = try RealVector(T).init(limits.len, self.allocator); defer q.deinit();
            var qPlus  = try RealVector(T).init(limits.len, self.allocator); defer qPlus.deinit();
            var qMinus = try RealVector(T).init(limits.len, self.allocator); defer qMinus.deinit();

            var dq: T = 1;

            for (0..limits.len) |i| {
                dq *= (limits[i][1] - limits[i][0]) / @as(T, @floatFromInt(npoint - 1));
            }

            for (0..std.math.pow(usize, npoint, limits.len)) |k| {

                positionAtRow(T, &q, k, limits.len, npoint, limits);

                try q.copyTo(&qPlus); qPlus.ptr(index).* += FINITE_DIFFERENCES_STEP;
                try q.copyTo(&qMinus); qMinus.ptr(index).* -= FINITE_DIFFERENCES_STEP;

                pot.evaluateDiabatic(&vPlus, qPlus, time);
                pot.evaluateDiabatic(&vMinus, qMinus, time);

                for (0..vDerivative.rows) |i| for (0..vDerivative.cols) |j| {
                    vDerivative.ptr(i, j).* = (vPlus.at(i, j) - vMinus.at(i, j)) / (2 * FINITE_DIFFERENCES_STEP);
                };

                const g1_value = self.evaluate(q.data);
                const g2_value = other.evaluate(q.data);

                for (0..V.rows) |i| for (0..V.cols) |j| {
                    V.ptr(i, j).* = V.at(i, j).add(g1_value.conjugate().mul(Complex(T).init(vDerivative.at(i, j), 0)).mul(g2_value));
                };
            }

            V.muls(Complex(T).init(dq, 0));

            return V;
        }

        /// Returns the potential energy expectation value of this complex Gaussian for a given potential operator.
        pub fn potentialEnergy(self: @This(), pot: ElectronicPotential(T), coefs: ComplexVector(T), limits: []const []const T, npoint: usize, time: T) !T {
            const V = try self.potential(self, pot, limits, npoint, time); defer V.deinit();

            var potential_energy: Complex(T) = Complex(T).init(0, 0);

            for (0..coefs.len) |i| for (0..coefs.len) |j| {
                potential_energy = potential_energy.add(coefs.at(i).conjugate().mul(V.at(i, j)).mul(coefs.at(j)));
            };

            return potential_energy.re;
        }
    };
}
