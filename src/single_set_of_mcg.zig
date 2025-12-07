//! File that wraps a set of complex Gaussian functions into a single structure for sharing between multiple wavefunctions.

const std = @import("std");

const complex_gaussian = @import("complex_gaussian.zig");
const complex_matrix = @import("complex_matrix.zig");
const complex_vector = @import("complex_vector.zig");
const electronic_potential = @import("electronic_potential.zig");
const error_handlingg = @import("error_handling.zig");
const matrix_inverse = @import("matrix_inverse.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const object_array = @import("object_array.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const Complex = std.math.Complex;
const ComplexGaussian = complex_gaussian.ComplexGaussian;
const ComplexMatrix = complex_matrix.ComplexMatrix;
const ComplexMatrixArray = object_array.ComplexMatrixArray;
const ComplexVector = complex_vector.ComplexVector;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const inverseHermitianAlloc = matrix_inverse.inverseHermitianAlloc;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const mm = matrix_multiplication.mm;
const mmAlloc = matrix_multiplication.mmAlloc;
const throw = error_handlingg.throw;

/// A linear combination of complex Gaussian functions shared between multiple wavefunctions.
pub fn SingleSetOfMCG(comptime T: type) type {
    return struct {
        gaussians: []ComplexGaussian(T),
        coefs: ComplexVector(T),
        allocator: std.mem.Allocator,

        /// Initialize a new complex Gaussian with given parameters.
        pub fn init(position: []const []const T, gamma: []const []const T, momentum: []const []const T, state: usize, bf_spread: []const u32, nstate: usize, allocator: std.mem.Allocator) !@This() {
            if (position.len != gamma.len or position.len != momentum.len or gamma.len != momentum.len) {
                return throw(@This(), "POSITION, GAMMA, AND MOMENTUM ARRAYS MUST HAVE THE SAME LENGTH", .{});
            }

            if (position.len != bf_spread.len) {
                return throw(@This(), "BASIS FUNCTION SPREAD LENGTH MUST MATCH NUMBER OF GAUSSIANS", .{});
            }

            var gaussians = try allocator.alloc(ComplexGaussian(T), position.len); var coefs = try ComplexVector(T).init(nstate * gaussians.len, allocator);

            for (0..gaussians.len) |i| {
                gaussians[i] = try ComplexGaussian(T).init(position[i], gamma[i], momentum[i], allocator);
            }

            for (bf_spread, 0..) |contr, i| {
                coefs.ptr(state * gaussians.len + i).* = Complex(T).init(@as(T, @floatFromInt(contr)), 0);
            }

            var mcg = @This(){
                .gaussians = gaussians,
                .coefs = coefs,
                .allocator = allocator,
            };

            try mcg.normalize();

            return mcg;
        }

        /// Deinitialize the complex Gaussian array.
        pub fn deinit(self: @This()) void {
            for (self.gaussians) |gaussian| gaussian.deinit();
            self.allocator.free(self.gaussians);
            self.coefs.deinit();
        }

        /// Calculates the derivative of the coefficients with respect to time.
        pub fn coefficientDerivative(self: @This(), pot: ElectronicPotential(T), dq: RealVector(T), dp: RealVector(T), mass: []const T, integration_nodes: usize, time: T) !ComplexVector(T) {
            const S = try self.overlap(); defer S.deinit();
            const H = try self.hamiltonian(mass, pot, integration_nodes, time); defer H.deinit();
            const tau = try self.overlapDiffTime(dq, dp); defer tau.deinit();
            const Sinv = try inverseHermitianAlloc(T, S, self.allocator); defer Sinv.deinit();

            var Heff = try ComplexMatrix(T).initZero(self.coefs.len, self.coefs.len, self.allocator); defer Heff.deinit();

            for (0..Heff.rows) |i| for (0..Heff.cols) |j| {
                Heff.ptr(i, j).* = H.at(i, j).sub(tau.at(i, j).mulbyi());
            };

            const Hc = try mmAlloc(T, Heff, false, self.coefs.asMatrix(), false, self.allocator); defer Hc.deinit();
            var dc = try mmAlloc(T, Sinv, false, Hc, false, self.allocator);

            dc.muls(Complex(T).init(0, -1));

            return ComplexVector(T){.data = dc.data, .len = dc.rows, .allocator = self.allocator};
        }

        /// Calculates the derivative of the coefficients with respect to imaginary time.
        pub fn coefficientDerivativeImaginary(self: @This(), pot: ElectronicPotential(T), dq: RealVector(T), dp: RealVector(T), mass: []const T, integration_nodes: usize, time: T) !ComplexVector(T) {
            const S = try self.overlap(); defer S.deinit();
            const H = try self.hamiltonian(mass, pot, integration_nodes, time); defer H.deinit();
            const tau = try self.overlapDiffTime(dq, dp); defer tau.deinit();
            const Sinv = try inverseHermitianAlloc(T, S, self.allocator); defer Sinv.deinit();

            var Heff = try ComplexMatrix(T).initZero(self.coefs.len, self.coefs.len, self.allocator); defer Heff.deinit();

            for (0..Heff.rows) |i| for (0..Heff.cols) |j| {
                Heff.ptr(i, j).* = H.at(i, j).add(tau.at(i, j));
            };

            const Hc = try mmAlloc(T, Heff, false, self.coefs.asMatrix(), false, self.allocator); defer Hc.deinit();
            var dc = try mmAlloc(T, Sinv, false, Hc, false, self.allocator);

            dc.muls(Complex(T).init(-1, 0));

            return ComplexVector(T){.data = dc.data, .len = dc.rows, .allocator = self.allocator};
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

        /// Exports the parameters of the Gaussians and coefficients into a single complex vector. The parameters are in order: positions, momentas, gammas for each Gaussian followed by the coefficients for each state and Gaussian.
        pub fn exportParameterVector(self: @This()) !ComplexVector(T) {
            const params = try ComplexVector(T).initZero(3 * self.gaussians[0].position.len * self.gaussians.len + self.coefs.len, self.allocator);

            for (self.gaussians, 0..) |gaussian, i| for (0..gaussian.position.len) |j| {

                const index = 3 * i * gaussian.position.len + j;

                params.ptr(index).* = Complex(T).init(gaussian.position[j], 0);
                params.ptr(gaussian.position.len + index).* = Complex(T).init(gaussian.momentum[j], 0);
                params.ptr(2 * gaussian.position.len + index).* = gaussian.gamma[j];
            };

            for (0..self.coefs.len) |i| {
                params.ptr(3 * self.gaussians[0].position.len * self.gaussians.len + i).* = self.coefs.at(i);
            }

            return params;
        }

        /// Calculates the Ehrenfest-like gamma derivative.
        pub fn gammaDerivativeEhrenfest(self: @This(), pot: ElectronicPotential(T), mass: []const T, integration_nodes: usize, time: T, fdiff_step: T) !ComplexVector(T) {
            const ddV = try ComplexMatrixArray(T).init(self.gaussians[0].position.len, .{.rows = pot.nstate(), .cols = pot.nstate()}, self.allocator); defer ddV.deinit();

            var dg = try ComplexVector(T).initZero(self.gaussians.len * self.gaussians[0].gamma.len, self.allocator); 

            for (self.gaussians, 0..) |gaussian, i| {

                var dg_i = try ComplexVector(T).initZero(gaussian.gamma.len, self.allocator); defer dg_i.deinit();

                for (0..gaussian.gamma.len) |j| {
                    dg_i.ptr(j).* = gaussian.gamma[j].mul(gaussian.gamma[j]).div(Complex(T).init(mass[j], 0));
                }

                for (self.gaussians, 0..) |other, j| {

                    for (0..ddV.len) |k| {
                        ddV.ptr(k).deinit(); ddV.ptr(k).* = try self.gaussians[i].potentialDerivative2(other, pot, k, integration_nodes, time, fdiff_step);
                    }

                    for (0..dg_i.len) |k| for (0..self.coefs.len / self.gaussians.len) |l| for (0..self.coefs.len / self.gaussians.len) |m| {

                        const c1 = self.coefs.at(l * self.gaussians.len + i); const c2 = self.coefs.at(m * self.gaussians.len + j);

                        dg_i.ptr(k).* = dg_i.at(k).sub(c1.conjugate().mul(ddV.at(k).at(l, m)).mul(c2));
                    };
                }

                dg_i.muls(Complex(T).init(0, -1));

                for (0..gaussian.momentum.len) |j| {
                    dg.ptr(gaussian.momentum.len * i + j).* = dg_i.at(j);
                }
            }

            return dg;
        }

        /// Calculates the derivative of the gamma parameter with respect to imaginary time in an Ehrenfest-like manner.
        pub fn gammaDerivativeImaginaryEhrenfest(self: @This(), pot: ElectronicPotential(T), mass: []const T, integration_nodes: usize, time: T, fdiff_step: T) !ComplexVector(T) {
            var dg = try self.gammaDerivativeEhrenfest(pot, mass, integration_nodes, time, fdiff_step);

            for (0..dg.len) |i| {
                dg.ptr(i).* = dg.at(i).mulbyi().neg();
            }

            return dg;
        }

        /// Calculates the Hamiltonian matrix.
        pub fn hamiltonian(self: @This(), mass: []const T, pot: ElectronicPotential(T), integration_nodes: usize, time: T) !ComplexMatrix(T) {
            var H = try self.potential(pot, integration_nodes, time);

            const K = try self.kinetic(mass); defer K.deinit();

            try H.add(K);

            return H;
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

        /// Loads the structure from a complex vector with parameters. The parameters are in order: positions, momentas, gammas for each Gaussian followed by the coefficients for each state and Gaussian.
        pub fn loadParameterVector(self: *@This(), params: ComplexVector(T)) !void {
            if (params.len != 3 * self.gaussians[0].position.len * self.gaussians.len + self.coefs.len) return throw(void, "PARAMETER VECTOR LENGTH MISMATCH", .{});

            for (self.gaussians, 0..) |gaussian, i| for (0..gaussian.position.len) |j| {

                const index = 3 * i * gaussian.position.len + j;

                gaussian.position[j] = params.at(index).re;
                gaussian.momentum[j] = params.at(gaussian.position.len + index).re;
                gaussian.gamma[j] = params.at(2 * gaussian.position.len + index);
            };

            for (0..self.coefs.len) |i| {
                self.coefs.ptr(i).* = params.at(3 * self.gaussians[0].position.len * self.gaussians.len + i);
            }
        }

        /// Calculates the Ehrenfest-like momentum derivative.
        pub fn momentumDerivativeEhrenfest(self: @This(), pot: ElectronicPotential(T), integration_nodes: usize, time: T, fdiff_step: T) !RealVector(T) {
            const dV = try ComplexMatrixArray(T).init(self.gaussians[0].position.len, .{.rows = pot.nstate(), .cols = pot.nstate()}, self.allocator); defer dV.deinit();

            var F = try RealVector(T).initZero(self.gaussians.len * self.gaussians[0].momentum.len, self.allocator); 

            for (self.gaussians, 0..) |gaussian, i| {

                var F_i = try RealVector(T).initZero(self.gaussians[0].momentum.len, self.allocator); defer F_i.deinit();

                for (self.gaussians, 0..) |other, j| {

                    for (0..dV.len) |k| {
                        dV.ptr(k).deinit(); dV.ptr(k).* = try self.gaussians[i].potentialDerivative1(other, pot, k, integration_nodes, time, fdiff_step);
                    }

                    for (0..F_i.len) |k| for (0..self.coefs.len / self.gaussians.len) |l| for (0..self.coefs.len / self.gaussians.len) |m| {

                        const c1 = self.coefs.at(l * self.gaussians.len + i); const c2 = self.coefs.at(m * self.gaussians.len + j);

                        F_i.ptr(k).* -= c1.conjugate().mul(dV.at(k).at(l, m)).mul(c2).re;
                    };
                }

                for (0..gaussian.momentum.len) |j| {
                    F.ptr(gaussian.momentum.len * i + j).* = F_i.at(j);
                }
            }

            return F;
        }

        /// Calculates the Ehrenfest-like momentum derivative in imaginary time.
        pub fn momentumDerivativeImaginaryEhrenfest(self: @This(), mass: []const T) !RealVector(T) {
            var dp = try self.positionDerivativeEhrenfest(mass); dp.muls(-1);

            return dp;
        }

        /// Normalizes the coefficient vector in place such that <Psi|Psi> = 1.
        pub fn normalize(self: *@This()) !void {
            const S = try self.overlap(); defer S.deinit();

            var norm_sq: T = 0;

            const SC = try mmAlloc(T, S, false, self.coefs.asMatrix(), false, self.allocator); defer SC.deinit();

            for (0..self.coefs.len) |i| {
                norm_sq += self.coefs.at(i).conjugate().mul(SC.at(i, 0)).re;
            }

            self.coefs.divs(Complex(T).init(std.math.sqrt(norm_sq), 0));
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

        /// Compute the overlap integral matrix where the ket vector is differentiated in time.
        pub fn overlapDiffTime(self: @This(), dq: RealVector(T), dp: RealVector(T)) !ComplexMatrix(T) {
            var tau = try ComplexMatrix(T).initZero(self.coefs.len, self.coefs.len, self.allocator);

            for (0..tau.rows) |i| for (0..tau.cols) |j| {

                const ig = i % self.gaussians.len; const is = i / self.gaussians.len;
                const jg = j % self.gaussians.len; const js = j / self.gaussians.len;

                if (is == js) {

                    const dq_j = dq.slice(jg * self.gaussians[0].position.len, (jg + 1) * self.gaussians[0].position.len);
                    const dp_j = dp.slice(jg * self.gaussians[0].position.len, (jg + 1) * self.gaussians[0].position.len);

                    tau.ptr(i, j).* = try self.gaussians[ig].overlapDiffTime(self.gaussians[jg], dq_j, dp_j);
                }
            };

            return tau;
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

        /// Calculates the Ehrenfest-like position derivative.
        pub fn positionDerivativeEhrenfest(self: @This(), mass: []const T) !RealVector(T) {
            var dq = try RealVector(T).initZero(self.gaussians.len, self.allocator);

            for (self.gaussians, 0..) |gaussian, i| {

                const dq_i = try gaussian.positionDerivative(mass); defer dq_i.deinit();

                for (0..gaussian.position.len) |j| {
                    dq.ptr(gaussian.position.len * i + j).* = dq_i.at(j);
                }
            }

            return dq;
        }

        /// Calculates the Ehrenfest-like position derivative in imaginary time.
        pub fn positionDerivativeImaginaryEhrenfest(self: @This(), pot: ElectronicPotential(T), integration_nodes: usize, time: T, fdiff_step: T) !RealVector(T) {
            return try self.momentumDerivativeEhrenfest(pot, integration_nodes, time, fdiff_step);
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
