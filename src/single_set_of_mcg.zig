//! File that wraps a set of complex Gaussian functions into a single structure for sharing between multiple wavefunctions.

const std = @import("std");

const complex_gaussian = @import("complex_gaussian.zig");
const complex_matrix = @import("complex_matrix.zig");
const complex_vector = @import("complex_vector.zig");
const electronic_potential = @import("electronic_potential.zig");
const error_handlingg = @import("error_handling.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const Complex = std.math.Complex;
const ComplexGaussian = complex_gaussian.ComplexGaussian;
const ComplexMatrix = complex_matrix.ComplexMatrix;
const ComplexVector = complex_vector.ComplexVector;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const mm = matrix_multiplication.mm;
const throw = error_handlingg.throw;

/// A linear combination of complex Gaussian functions shared between multiple wavefunctions.
pub fn SingleSetOfMCG(comptime T: type) type {
    return struct {
        gaussians: []ComplexGaussian(T),
        coefs: ComplexVector(T),
        allocator: std.mem.Allocator,

        /// Initialize a new complex Gaussian with given parameters.
        pub fn init(position: []const []const T, gamma: []const []const T, momentum: []const []const T, state: usize, nstate: usize, allocator: std.mem.Allocator) !@This() {
            if (position.len != gamma.len or position.len != momentum.len or gamma.len != momentum.len) {
                return throw(@This(), "POSITION, GAMMA, AND MOMENTUM ARRAYS MUST HAVE THE SAME LENGTH", .{});
            }

            var gaussians = try allocator.alloc(ComplexGaussian(T), position.len); var coefs = try ComplexVector(T).init(nstate * gaussians.len, allocator);

            for (0..gaussians.len) |i| {
                gaussians[i] = try ComplexGaussian(T).init(position[i], gamma[i], momentum[i], allocator);
            }

            coefs.ptr(state * gaussians.len).* = Complex(T).init(1, 0);

            return @This(){
                .gaussians = gaussians,
                .coefs = coefs,
                .allocator = allocator,
            };
        }

        /// Deinitialize the complex Gaussian array.
        pub fn deinit(self: @This()) void {
            for (self.gaussians) |gaussian| gaussian.deinit();
            self.allocator.free(self.gaussians);
            self.coefs.deinit();
        }

        /// Calculates the expectation value of position or momentum.
        pub fn coordinateExpectation(self: @This(), coord: enum {momentum, position}) !RealVector(T) {
            var coordinate = try RealVector(T).initZero(self.gaussians[0].momentum.len, self.allocator); var total_norm: T = 0;

            const S = try self.overlap(); defer S.deinit();

            for (0..self.coefs.len / self.gaussians.len) |i| {

                for (0..self.gaussians.len) |j| {

                    const index_1 = i * self.gaussians.len + j;

                    for (0..self.gaussians.len) |k| {

                        const index_2 = i * self.gaussians.len + k;

                        const c_1 = self.coefs.at(index_1);
                        const c_2 = self.coefs.at(index_2);

                        total_norm += c_1.conjugate().mul(c_2).mul(S.at(index_1, index_2)).re;

                        const me = switch (coord) {
                            .momentum => try self.gaussians[j].momentumMatrixElement(self.gaussians[k]),
                            .position => try self.gaussians[j].positionMatrixElement(self.gaussians[k])
                        }; defer me.deinit();

                        for (0..coordinate.len) |l| {
                            coordinate.ptr(l).* += c_1.conjugate().mul(c_2).mul(me.at(l)).re;
                        }
                    }
                }
            }

            coordinate.muls(total_norm);

            return coordinate;
        }

        /// Calculates the kinetic matrix.
        pub fn kinetic(self: @This(), mass: []const T) !ComplexMatrix(T) {
            var K = try ComplexMatrix(T).initZero(self.coefs.len, self.coefs.len, self.allocator);

            for (0..K.rows) |i| for (0..K.cols) |j| {

                const ig = i % self.gaussians.len; const is = i / self.gaussians.len;
                const jg = j % self.gaussians.len; const js = j / self.gaussians.len;

                if (is == js) {
                    K.ptr(i, j).* = try self.gaussians[ig].kinetic(self.gaussians[jg], mass);
                }
            };

            return K;
        }

        /// Calculates the kinetic energy.
        pub fn kineticEnergy(self: @This(), mass: []const T) !T {
            const K = try self.kinetic(mass); defer K.deinit();

            return try self.operatorExpectation(K);
        }

        /// Calculate the expectation value of an operator given its matrix representation.
        pub fn operatorExpectation(self: @This(), M: ComplexMatrix(T)) !T {
            var Mexp: T = 0; var norm: T = 0;

            const S = try self.overlap(); defer S.deinit();

            var temp1 = try ComplexVector(T).initZero(self.coefs.len, self.allocator); defer temp1.deinit(); var temp1_matrix = temp1.asMatrix();
            var temp2 = try ComplexVector(T).initZero(self.coefs.len, self.allocator); defer temp2.deinit(); var temp2_matrix = temp2.asMatrix();

            try mm(T, &temp1_matrix, M, false, self.coefs.asMatrix(), false);
            try mm(T, &temp2_matrix, S, false, self.coefs.asMatrix(), false);

            for (0..self.coefs.len) |i| {
                Mexp += self.coefs.at(i).conjugate().mul(temp1.at(i)).re;
                norm += self.coefs.at(i).conjugate().mul(temp2.at(i)).re;
            }

            return Mexp / norm;
        }

        /// Calculates the overlap matrix.
        pub fn overlap(self: @This()) !ComplexMatrix(T) {
            var S = try ComplexMatrix(T).initZero(self.coefs.len, self.coefs.len, self.allocator);

            for (0..S.rows) |i| for (0..S.cols) |j| {

                const ig = i % self.gaussians.len; const is = i / self.gaussians.len;
                const jg = j % self.gaussians.len; const js = j / self.gaussians.len;

                if (is == js) {
                    S.ptr(i, j).* = try self.gaussians[ig].overlap(self.gaussians[jg]);
                }
            };

            return S;
        }

        /// Calculates the population of a given state. If the coefficients are passed, those are used instead of the internal ones.
        pub fn population(self : @This(), state: usize, coefs_to_use: ?ComplexVector(T)) !T {
            var pop: T = 0;

            const coefs = coefs_to_use orelse self.coefs;

            const S = try self.overlap(); defer S.deinit();

            for (0..self.gaussians.len) |i| {

                const index1 = state * self.gaussians.len + i;

                for (0..self.gaussians.len) |j| {

                    const index2 = state * self.gaussians.len + j;

                    const c1 = coefs.at(index1);
                    const c2 = coefs.at(index2);

                    pop += c1.conjugate().mul(c2).mul(S.at(index1, index2)).re;
                }
            }

            return pop;
        }

        /// Calculates the potential matrix for a given potential function.
        pub fn potential(self: @This(), pot: ElectronicPotential(T), integration_nodes: usize, time: T) !ComplexMatrix(T) {
            var V = try ComplexMatrix(T).initZero(self.coefs.len, self.coefs.len, self.allocator);

            for (0..self.gaussians.len) |i| for (0..self.gaussians.len) |j| {

                var Vij = try self.gaussians[i].potential(self.gaussians[j], pot, integration_nodes, time); defer Vij.deinit();

                for (0..Vij.rows) |k| for (0..Vij.cols) |l| {

                    const row = k * self.gaussians.len + i;
                    const col = l * self.gaussians.len + j;

                    V.ptr(row, col).* = V.at(row, col).add(Vij.at(k, l));
                };
            };

            return V;
        }

        /// Calculates the potential energy.
        pub fn potentialEnergy(self: @This(), pot: ElectronicPotential(T), integration_nodes: usize, time: T) !T {
            const V = try self.potential(pot, integration_nodes, time); defer V.deinit();

            return try self.operatorExpectation(V);
        }

        /// Returns the transformed coefficients. If to_adiabatic is true, transforms to adiabatic representation; otherwise, to diabatic.
        pub fn transformedCoefs(self: @This(), pot: ElectronicPotential(T), time: T, to_adiabatic: bool) !ComplexVector(T) {
            const coefs_after = try ComplexVector(T).initZero(self.coefs.len, self.allocator);

            var diabatic_potential = try RealMatrix(T).init(self.coefs.len / self.gaussians.len, self.coefs.len / self.gaussians.len, self.allocator); defer diabatic_potential.deinit();
            var adiabatic_potential = try RealMatrix(T).init(self.coefs.len / self.gaussians.len, self.coefs.len / self.gaussians.len, self.allocator); defer adiabatic_potential.deinit();
            var adiabatic_eigenvectors = try RealMatrix(T).init(self.coefs.len / self.gaussians.len, self.coefs.len / self.gaussians.len, self.allocator); defer adiabatic_eigenvectors.deinit();

            for (self.gaussians, 0..) |gaussian, i| {

                const position = RealVector(T){.data = gaussian.position, .len = gaussian.position.len, .allocator = null};

                const coefs_before_i = try ComplexVector(T).initZero(self.coefs.len / self.gaussians.len, self.allocator); defer coefs_before_i.deinit(); const coefs_before_i_matrix = coefs_before_i.asMatrix();
                const coefs_after_i  = try ComplexVector(T).initZero(self.coefs.len / self.gaussians.len, self.allocator); defer  coefs_after_i.deinit(); var    coefs_after_i_matrix =  coefs_after_i.asMatrix();

                for (0..self.coefs.len / self.gaussians.len) |j| {
                    coefs_before_i.ptr(j).* = self.coefs.at(j * self.gaussians.len + i);
                }

                try pot.evaluateEigensystem(&diabatic_potential, &adiabatic_potential, &adiabatic_eigenvectors, position, time);

                if (to_adiabatic) {
                    try mm(T, &coefs_after_i_matrix, adiabatic_eigenvectors, true, coefs_before_i_matrix, false);
                } else {
                    try mm(T, &coefs_after_i_matrix, adiabatic_eigenvectors, false, coefs_before_i_matrix, false);
                }

                for (0..coefs_after_i.len) |j| {
                    coefs_after.ptr(j * self.gaussians.len + i).* = coefs_after_i.at(j);
                }
            }

            return coefs_after;
        }
    };
}
