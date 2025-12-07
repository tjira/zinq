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
const printComplexVector = device_write.printComplexVector;
const prod = array_functions.prod;
const throw = error_handling.throw;

const HERMITE_NODES = hermite_quadrature_nodes.HERMITE_NODES;
const HERMITE_WEIGHTS = hermite_quadrature_nodes.HERMITE_WEIGHTS;

/// Complex Gaussian function implementation in Zig.
pub fn ComplexGaussian(comptime T: type) type {
    return struct {
        position: []T,
        gamma: []Complex(T),
        momentum: []T,
        allocator: std.mem.Allocator,

        /// Initialize a new complex Gaussian with given parameters.
        pub fn init(position: []const T, gamma: []const T, momentum: []const T, allocator: std.mem.Allocator) !@This() {
            if (position.len != gamma.len or position.len != momentum.len or gamma.len != momentum.len) {
                return throw(ComplexGaussian(T), "POSITION, GAMMA, AND MOMENTUM MUST HAVE THE SAME DIMENSIONALITY", .{});
            }

            var cg = @This(){
                .position = try allocator.alloc(T, position.len),
                .gamma = try allocator.alloc(Complex(T), position.len),
                .momentum = try allocator.alloc(T, position.len),
                .allocator = allocator
            };

            for (position, 0..) |c, i| cg.position[i] = c;
            for (gamma, 0..) |g, i| cg.gamma[i] = Complex(T).init(g, 0);
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

        /// Returns the derivative of the electronic coefficients where this gaussian is shared between them.
        pub fn coefficientDerivativeImaginary(self: @This(), coefs: ComplexVector(T), pot: ElectronicPotential(T), n_nodes: usize, time: T) !ComplexVector(T) {
            var dc = try self.coefficientDerivative(coefs, pot, n_nodes, time);

            for (0..dc.len) |i| {
                dc.ptr(i).* = dc.at(i).mulbyi().neg();
            }

            return dc;
        }

        /// Evaluate the complex Gaussian function at a given point x.
        pub fn evaluate(self: @This(), q: []const T) Complex(T) {
            var result = Complex(T).init(1, 0);

            for (self.position, 0..) |c, i| {

                const dq_c = Complex(T).init(q[i] - c, 0); const p_c = Complex(T).init(0, self.momentum[i]);

                const exponent = Complex(T).init(0.5, 0).neg().mul(self.gamma[i].mul(dq_c.mul(dq_c))).add(p_c.mul(dq_c));

                result = result.mul(std.math.complex.exp(exponent));
            }

            return result.div(Complex(T).init(self.norm(), 0));
        }

        /// Calculates the derivative of the gamma parameter.
        pub fn gammaDerivative(self: @This(), pot: ElectronicPotential(T), coefs: ComplexVector(T), mass: []const T, n_nodes: usize, time: T, fdiff_step: T) !ComplexVector(T) {
            const ddV = try ComplexMatrixArray(T).init(self.position.len, .{.rows = pot.nstate(), .cols = pot.nstate()}, self.allocator); defer ddV.deinit();

            for (0..ddV.len) |i| {
                ddV.ptr(i).deinit(); ddV.ptr(i).* = try self.potentialDerivative2(self, pot, i, n_nodes, time, fdiff_step);
            }

            const dg = try ComplexVector(T).initZero(ddV.len, self.allocator);

            for (0..dg.len) |i| {

                dg.ptr(i).* = self.gamma[i].mul(self.gamma[i]).div(Complex(T).init(mass[i], 0));

                for (0..coefs.len) |j| for (0..coefs.len) |k| {
                    dg.ptr(i).* = dg.at(i).sub(coefs.at(j).conjugate().mul(ddV.at(i).at(j, k)).mul(coefs.at(k)));
                };

                dg.ptr(i).* = dg.at(i).mulbyi().neg();
            }

            return dg;
        }

        /// Calculates the derivative of the gamma parameter.
        pub fn gammaDerivativeImaginary(self: @This(), pot: ElectronicPotential(T), coefs: ComplexVector(T), mass: []const T, n_nodes: usize, time: T, fdiff_step: T) !ComplexVector(T) {
            var dg = try self.gammaDerivative(pot, coefs, mass, n_nodes, time, fdiff_step);

            for (0..dg.len) |i| {
                dg.ptr(i).* = dg.at(i).mulbyi().neg();
            }

            return dg;
        }

        /// Calculate the kinetic matrix element between this complex Gaussian and another.
        pub fn kinetic(self: @This(), other: @This(), mass: []const T) !Complex(T) {
            if (self.position.len != other.position.len or self.position.len != mass.len or other.position.len != mass.len) {
                return throw(Complex(T), "BOTH POSITION AND MASS ARRAY MUST HAVE THE SAME DIMENSIONALITY", .{});
            }

            var result = Complex(T).init(0, 0);

            for (self.position, other.position, self.gamma, other.gamma, self.momentum, other.momentum, mass) |q1, q2, g1, g2, p1, p2, m| {

                const g1_c = g1.conjugate();

                const p1_c = Complex(T).init(p1, 0); const p2_c = Complex(T).init(p2, 0);

                const dq_c = Complex(T).init(q1 - q2, 0);

                const m_c = Complex(T).init(m, 0);

                const factor1 = std.math.complex.sqrt(Complex(T).init(0.25, 0).mul(g1_c.add(g2))).div(std.math.complex.pow(g1_c.add(g2), Complex(T).init(2.5, 0))).div(m_c);
                const factor2 = std.math.complex.pow(g2.mul(p1_c).add(g1_c.mul(p2_c)), Complex(T).init(2, 0)).add(g1_c.mul(g2).mul(g1_c.add(g2).sub(g1_c.mul(g2).mul(dq_c).mul(dq_c))));
                const factor3 = Complex(T).init(2, 0).mul(g1_c).mul(g2).mul(dq_c).mul(g2.mul(p1_c).add(g1_c.mul(p2_c))).mulbyi();

                result = result.add(factor2.add(factor3).mul(factor1));
            }

            return result.mul(try self.overlap(other));
        }

        /// Returns the kinetic energy expectation value of this complex Gaussian.
        pub fn kineticEnergy(self: @This(), mass: []const T) !T {
            if (self.position.len != mass.len) return throw(T, "POSITION AND MASS ARRAY MUST HAVE THE SAME DIMENSIONALITY IN COMPLEX GAUSSIAN", .{});

            return (try kinetic(self, self, mass)).re;
        }

        /// Calculates the derivative of the momentum.
        pub fn momentumDerivative(self: @This(), pot: ElectronicPotential(T), coefs: ComplexVector(T), n_nodes: usize, time: T, fdiff_step: T) !RealVector(T) {
            const dV = try ComplexMatrixArray(T).init(self.position.len, .{.rows = pot.nstate(), .cols = pot.nstate()}, self.allocator); defer dV.deinit();

            for (0..dV.len) |i| {
                dV.ptr(i).deinit(); dV.ptr(i).* = try self.potentialDerivative1(self, pot, i, n_nodes, time, fdiff_step);
            }

            var F = try RealVector(T).initZero(dV.len, self.allocator);

            for (0..F.len) |i| for (0..coefs.len) |j| for (0..coefs.len) |k| {
                F.ptr(i).* -= coefs.at(j).conjugate().mul(dV.at(i).at(j, k)).mul(coefs.at(k)).re;
            };

            return F;
        }

        /// Calculates the derivative of the position when performing imaginaty time propagation.
        pub fn momentumDerivativeImaginary(self: @This(), mass: []const T) !RealVector(T) {
            var dp = try self.positionDerivative(mass); dp.muls(-1);

            return dp;
        }

        /// Calculates the matrix element of the momentum operator between this complex Gaussian and another.
        pub fn momentumMatrixElement(self: @This(), other: @This()) !ComplexVector(T) {
            if (self.position.len != other.position.len) return throw(ComplexVector(T), "BOTH COMPLEX GAUSSIANS MUST HAVE THE SAME DIMENSIONALITY", .{});

            var result = try ComplexVector(T).initZero(self.position.len, self.allocator);

            for (self.position, other.position, self.gamma, other.gamma, self.momentum, other.momentum, 0..) |q1, q2, g1, g2, p1, p2, i| {

                const g1_c = g1.conjugate();

                const dq_c = Complex(T).init(q1 - q2, 0);

                result.ptr(i).* = g1_c.mul(Complex(T).init(p2, 0)).add(g2.mul(Complex(T).init(p1, 0))).add(g1_c.mul(g2).mul(dq_c).mulbyi()).div(g1_c.add(g2));
            }

            result.muls(try self.overlap(other));

            return result;
        }

        /// Calculates the matrix element of the position operator between this complex Gaussian and another.
        pub fn positionMatrixElement(self: @This(), other: @This()) !ComplexVector(T) {
            if (self.position.len != other.position.len) return throw(ComplexVector(T), "BOTH COMPLEX GAUSSIANS MUST HAVE THE SAME DIMENSIONALITY", .{});

            var result = try ComplexVector(T).initZero(self.position.len, self.allocator);

            for (self.position, other.position, self.gamma, other.gamma, self.momentum, other.momentum, 0..) |q1, q2, g1, g2, p1, p2, i| {

                const g1_c = g1.conjugate();

                const dq_c = Complex(T).init(q1 - q2, 0);
                const dp_c = Complex(T).init(p2 - p1, 0);
                
                result.ptr(i).* = Complex(T).init(q2, 0).add(g1_c.div(g1_c.add(g2)).mul(dq_c)).add(dp_c.div(g1_c.add(g2)).mulbyi());
            }

            result.muls(try self.overlap(other));

            return result;
        }

        /// Returns the norm of the gaussian.
        pub fn norm(self: @This()) T {
            var gamma_prod: T = 1; for (self.gamma) |g| gamma_prod *= g.re;

            return std.math.pow(T, std.math.pi, @as(T, @floatFromInt(self.position.len)) / 4) / std.math.pow(T, gamma_prod, 0.25);
        }

        /// Compute the overlap integral between this complex Gaussian and another.
        pub fn overlap(self: @This(), other: @This()) !Complex(T) {
            if (self.position.len != other.position.len) return throw(Complex(T), "BOTH COMPLEX GAUSSIANS MUST HAVE THE SAME DIMENSIONALITY", .{});

            var result = Complex(T).init(1, 0);

            for (self.position, other.position, self.gamma, other.gamma, self.momentum, other.momentum) |q1, q2, g1, g2, p1, p2| {

                const p1_c = Complex(T).init(p1, 0); const p2_c = Complex(T).init(p2, 0);

                const dq_c = Complex(T).init(q1 - q2, 0); const dp_c = Complex(T).init(p1 - p2, 0);

                const factor = std.math.complex.sqrt(Complex(T).init(2 * std.math.pi, 0).div(g1.conjugate().add(g2)));

                const num1 = dp_c.mul(dp_c).add(g1.conjugate().mul(g2).mul(dq_c).mul(dq_c));
                const num2 = Complex(T).init(-2, 0).mul(dq_c).mul(g2.mul(p1_c).add(g1.conjugate().mul(p2_c))).mulbyi();

                const exponent = num1.add(num2).div(Complex(T).init(-2, 0).mul(g1.conjugate().add(g2)));

                result = result.mul(std.math.complex.exp(exponent).mul(factor));
            }

            return result.div(Complex(T).init(self.norm() * other.norm(), 0));
        }

        /// Compute the overlap integral between this complex Gaussian and another which is differentiated with respect to momentum.
        pub fn overlapDiffMomentum(self: @This(), other: @This()) !Complex(T) {
            if (self.position.len != other.position.len) return throw(Complex(T), "BOTH COMPLEX GAUSSIANS MUST HAVE THE SAME DIMENSIONALITY", .{});

            var result = Complex(T).init(0, 0);

            for (self.position, other.position, self.gamma, other.gamma, self.momentum, other.momentum) |q1, q2, g1, g2, p1, p2| {

                const g1_c = g1.conjugate();

                const p1_c = Complex(T).init(p1, 0); const p2_c = Complex(T).init(p2, 0);

                const dq_c = Complex(T).init(q1 - q2, 0);

                result = result.add(p1_c.sub(p2_c).add(g1_c.mul(dq_c).mulbyi()).div(g1_c.add(g2)));
            }

            return result.mul(try self.overlap(other));
        }

        /// Compute the overlap integral between this complex Gaussian and another which is differentiated with respect to position.
        pub fn overlapDiffPosition(self: @This(), other: @This()) !Complex(T) {
            if (self.position.len != other.position.len) return throw(Complex(T), "BOTH COMPLEX GAUSSIANS MUST HAVE THE SAME DIMENSIONALITY", .{});

            var result = Complex(T).init(0, 0);

            for (self.position, other.position, self.gamma, other.gamma, self.momentum, other.momentum) |q1, q2, g1, g2, p1, p2| {

                const g1_c = g1.conjugate();

                const p1_c = Complex(T).init(p1, 0); const p2_c = Complex(T).init(p2, 0);

                const dq_c = Complex(T).init(q1 - q2, 0);

                result = result.add(g1_c.mul(p2_c).neg().mulbyi().add(g2.mul(p1_c.mulbyi().neg().add(g1_c.mul(dq_c)))).div(g1_c.add(g2)));
            }

            return result.mul(try self.overlap(other));
        }

        /// Compute the overlap integral between this complex Gaussian and another which is differentiated with respect to time.
        pub fn overlapDiffTime(self: @This(), other: @This(), dq: RealVector(T), dp: RealVector(T)) !Complex(T) {
            if (self.position.len != other.position.len) return throw(Complex(T), "BOTH COMPLEX GAUSSIANS MUST HAVE THE SAME DIMENSIONALITY", .{});

            var result = Complex(T).init(0, 0);

            for (self.position, other.position, self.gamma, other.gamma, self.momentum, other.momentum, 0..) |q1, q2, g1, g2, p1, p2, i| {

                const g1_c = g1.conjugate();

                const p1_c = Complex(T).init(p1, 0); const p2_c = Complex(T).init(p2, 0);

                const dq_c = Complex(T).init(q1 - q2, 0);

                result = result.add(g1_c.mul(p2_c).neg().mulbyi().add(g2.mul(p1_c.mulbyi().neg().add(g1_c.mul(dq_c)))).div(g1_c.add(g2)).mul(Complex(T).init(dq.at(i), 0)));
            }

            for (self.position, other.position, self.gamma, other.gamma, self.momentum, other.momentum, 0..) |q1, q2, g1, g2, p1, p2, i| {

                const g1_c = g1.conjugate();

                const p1_c = Complex(T).init(p1, 0); const p2_c = Complex(T).init(p2, 0);

                const dq_c = Complex(T).init(q1 - q2, 0);

                result = result.add(p1_c.sub(p2_c).add(g1_c.mul(dq_c).mulbyi()).div(g1_c.add(g2)).mul(Complex(T).init(dp.at(i), 0)));
            }

            return result.mul(try self.overlap(other));
        }

        /// Calculates the derivative of the position.
        pub fn positionDerivative(self: @This(), mass: []const T) !RealVector(T) {
            var dq = try RealVector(T).init(self.position.len, self.allocator);

            for (0..dq.len) |i| {
                dq.ptr(i).* = self.momentum[i] / mass[i];
            }

            return dq;
        }

        /// Calculates the derivative of the position. in imaginary time propagation.
        pub fn positionDerivativeImaginary(self: @This(), pot: ElectronicPotential(T), coefs: ComplexVector(T), n_nodes: usize, time: T, fdiff_step: T) !RealVector(T) {
            return try self.momentumDerivative(pot, coefs, n_nodes, time, fdiff_step);
        }

        /// Compute the potential energy matrix element between this complex Gaussian and another for a given potential function.
        pub fn potential(self: @This(), other: @This(), pot: ElectronicPotential(T), n_nodes: usize, time: T) !ComplexMatrix(T) {
            var V = try ComplexMatrix(T).initZero(pot.nstate(), pot.nstate(), self.allocator);

            var variables = try self.allocator.alloc(T, pot.ndim()); defer self.allocator.free(variables);
            var sigmas = try self.allocator.alloc(T, pot.ndim()); defer self.allocator.free(sigmas);
            var mus = try self.allocator.alloc(T, pot.ndim()); defer self.allocator.free(mus);

            for (0..pot.ndim()) |i| {

                const s1 = 1 / std.math.sqrt(self.gamma[i].re); const s2 = 1 / std.math.sqrt(other.gamma[i].re);

                mus[i] = (self.position[i] * s2 * s2 + other.position[i] * s1 * s1) / (s1 * s1 + s2 * s2);

                sigmas[i] = (s1 * s2) / std.math.sqrt(s1 * s1 + s2 * s2);
            }

            for (0..std.math.pow(usize, n_nodes, pot.ndim())) |i| {

                var temp = i; var weight: T = 1; var complex_exponent = Complex(T).init(0, 0);

                for (0..pot.ndim()) |j| {

                    variables[j] = std.math.sqrt2 * sigmas[j] * HERMITE_NODES[n_nodes][temp % n_nodes] + mus[j]; 

                    const real_exponent = -0.5 * self.gamma[j].re * other.gamma[j].re / (self.gamma[j].re + other.gamma[j].re) * std.math.pow(T, self.position[j] - other.position[j], 2);

                    const imag_exponent_self = self.momentum[j] * (self.position[j] - variables[j]);
                    const imag_exponent_other = other.momentum[j] * (variables[j] - other.position[j]);

                    const imag_gamma_self = 0.5 * self.gamma[j].im * std.math.pow(T, variables[j] - self.position[j], 2);
                    const imag_gamma_other = -0.5 * other.gamma[j].im * std.math.pow(T, variables[j] - other.position[j], 2);

                    complex_exponent = complex_exponent.add(Complex(T).init(real_exponent, imag_exponent_self + imag_exponent_other + imag_gamma_self + imag_gamma_other));

                    weight *= std.math.sqrt2 * sigmas[j] * HERMITE_WEIGHTS[n_nodes][temp % n_nodes]; temp /= n_nodes;
                }

                for (0..V.rows) |j| for (j..V.cols) |k| {

                    const variable_vector = RealVector(T){.data = variables, .len = variables.len, .allocator = null};

                    const value = try pot.evaluateDiabaticElement(j, k, variable_vector, time);

                    V.ptr(j, k).* = V.at(j, k).add(Complex(T).init(weight * value, 0).mul(std.math.complex.exp(complex_exponent))); V.ptr(k, j).* = V.at(j, k);
                };
            }

            V.divs(Complex(T).init(self.norm() * other.norm(), 0));

            return V;
        }

        /// Compute the potential energy derivative matrix element between this complex Gaussian and another for a given potential function and a specified coordinate.
        pub fn potentialDerivative1(self: @This(), other: @This(), pot: ElectronicPotential(T), index: usize, n_nodes: usize, time: T, fdiff_step: T) !ComplexMatrix(T) {
            var V = try ComplexMatrix(T).initZero(pot.nstate(), pot.nstate(), self.allocator);

            var variables = try self.allocator.alloc(T, pot.ndim()); defer self.allocator.free(variables);
            var sigmas = try self.allocator.alloc(T, pot.ndim()); defer self.allocator.free(sigmas);
            var mus = try self.allocator.alloc(T, pot.ndim()); defer self.allocator.free(mus);

            for (0..pot.ndim()) |i| {

                const s1 = 1 / std.math.sqrt(self.gamma[i].re); const s2 = 1 / std.math.sqrt(other.gamma[i].re);

                mus[i] = (self.position[i] * s2 * s2 + other.position[i] * s1 * s1) / (s1 * s1 + s2 * s2);

                sigmas[i] = (s1 * s2) / std.math.sqrt(s1 * s1 + s2 * s2);
            }

            for (0..std.math.pow(usize, n_nodes, pot.ndim())) |i| {

                var temp = i; var weight: T = 1; var complex_exponent = Complex(T).init(0, 0);

                for (0..pot.ndim()) |j| {

                    variables[j] = std.math.sqrt2 * sigmas[j] * HERMITE_NODES[n_nodes][temp % n_nodes] + mus[j]; 

                    const real_exponent = -0.5 * self.gamma[j].re * other.gamma[j].re / (self.gamma[j].re + other.gamma[j].re) * std.math.pow(T, self.position[j] - other.position[j], 2);

                    const imag_exponent_self = self.momentum[j] * (self.position[j] - variables[j]);
                    const imag_exponent_other = other.momentum[j] * (variables[j] - other.position[j]);

                    const imag_gamma_self = 0.5 * self.gamma[j].im * std.math.pow(T, variables[j] - self.position[j], 2);
                    const imag_gamma_other = -0.5 * other.gamma[j].im * std.math.pow(T, variables[j] - other.position[j], 2);

                    complex_exponent = complex_exponent.add(Complex(T).init(real_exponent, imag_exponent_self + imag_exponent_other + imag_gamma_self + imag_gamma_other));

                    weight *= std.math.sqrt2 * sigmas[j] * HERMITE_WEIGHTS[n_nodes][temp % n_nodes]; temp /= n_nodes;
                }

                for (0..V.rows) |j| for (j..V.cols) |k| {

                    const variable_vector = RealVector(T){.data = variables, .len = variables.len, .allocator = null};

                    const value = try pot.evaluateDiabaticElementDerivative1(j, k, variable_vector, time, index, fdiff_step);

                    V.ptr(j, k).* = V.at(j, k).add(Complex(T).init(weight * value, 0).mul(std.math.complex.exp(complex_exponent))); V.ptr(k, j).* = V.at(j, k);
                };
            }

            V.divs(Complex(T).init(self.norm() * other.norm(), 0));

            return V;
        }

        /// Compute the potential energy second derivative matrix element between this complex Gaussian and another for a given potential function and a specified coordinate.
        pub fn potentialDerivative2(self: @This(), other: @This(), pot: ElectronicPotential(T), index: usize, n_nodes: usize, time: T, fdiff_step: T) !ComplexMatrix(T) {
            var V = try ComplexMatrix(T).initZero(pot.nstate(), pot.nstate(), self.allocator);

            var variables = try self.allocator.alloc(T, pot.ndim()); defer self.allocator.free(variables);
            var sigmas = try self.allocator.alloc(T, pot.ndim()); defer self.allocator.free(sigmas);
            var mus = try self.allocator.alloc(T, pot.ndim()); defer self.allocator.free(mus);

            for (0..pot.ndim()) |i| {

                const s1 = 1 / std.math.sqrt(self.gamma[i].re); const s2 = 1 / std.math.sqrt(other.gamma[i].re);

                mus[i] = (self.position[i] * s2 * s2 + other.position[i] * s1 * s1) / (s1 * s1 + s2 * s2);

                sigmas[i] = (s1 * s2) / std.math.sqrt(s1 * s1 + s2 * s2);
            }

            for (0..std.math.pow(usize, n_nodes, pot.ndim())) |i| {

                var temp = i; var weight: T = 1; var complex_exponent = Complex(T).init(0, 0);

                for (0..pot.ndim()) |j| {

                    variables[j] = std.math.sqrt2 * sigmas[j] * HERMITE_NODES[n_nodes][temp % n_nodes] + mus[j]; 

                    const real_exponent = -0.5 * self.gamma[j].re * other.gamma[j].re / (self.gamma[j].re + other.gamma[j].re) * std.math.pow(T, self.position[j] - other.position[j], 2);

                    const imag_exponent_self = self.momentum[j] * (self.position[j] - variables[j]);
                    const imag_exponent_other = other.momentum[j] * (variables[j] - other.position[j]);

                    const imag_gamma_self = 0.5 * self.gamma[j].im * std.math.pow(T, variables[j] - self.position[j], 2);
                    const imag_gamma_other = -0.5 * other.gamma[j].im * std.math.pow(T, variables[j] - other.position[j], 2);

                    complex_exponent = complex_exponent.add(Complex(T).init(real_exponent, imag_exponent_self + imag_exponent_other + imag_gamma_self + imag_gamma_other));

                    weight *= std.math.sqrt2 * sigmas[j] * HERMITE_WEIGHTS[n_nodes][temp % n_nodes]; temp /= n_nodes;
                }

                for (0..V.rows) |j| for (j..V.cols) |k| {

                    const variable_vector = RealVector(T){.data = variables, .len = variables.len, .allocator = null};

                    const value = try pot.evaluateDiabaticElementDerivative2(j, k, variable_vector, time, index, fdiff_step);

                    V.ptr(j, k).* = V.at(j, k).add(Complex(T).init(weight * value, 0).mul(std.math.complex.exp(complex_exponent))); V.ptr(k, j).* = V.at(j, k);
                };
            }

            V.divs(Complex(T).init(self.norm() * other.norm(), 0));

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
