//! File that wraps a set of complex Gaussian functions into a single structure for sharing between multiple wavefunctions.

const std = @import("std");

const complex_gaussian = @import("complex_gaussian.zig");
const complex_matrix = @import("complex_matrix.zig");
const complex_vector = @import("complex_vector.zig");
const device_write = @import("device_write.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const electronic_potential = @import("electronic_potential.zig");
const error_handlingg = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");
const matrix_inverse = @import("matrix_inverse.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const object_array = @import("object_array.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const Complex = std.math.Complex;
const printRealMatrix = device_write.printRealMatrix;
const ComplexGaussian = complex_gaussian.ComplexGaussian;
const ComplexMatrix = complex_matrix.ComplexMatrix;
const ComplexMatrixArray = object_array.ComplexMatrixArray;
const ComplexVector = complex_vector.ComplexVector;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const printComplexMatrix = device_write.printComplexMatrix;
const pseudoInverseHermitianAlloc = matrix_inverse.pseudoInverseHermitianAlloc;
const inverseHermitianAlloc = matrix_inverse.inverseHermitianAlloc;
const fixGauge = eigenproblem_solver.fixGauge;
const mm = matrix_multiplication.mm;
const mmAlloc = matrix_multiplication.mmAlloc;
const throw = error_handlingg.throw;

const SINGULARITY_TOLERANCE = global_variables.SINGULARITY_TOLERANCE;

/// A linear combination of complex Gaussian functions shared between multiple wavefunctions.
pub fn SingleSetOfMCG(comptime T: type) type {
    return struct {
        gaussians: []ComplexGaussian(T),
        coefs: ComplexVector(T),

        /// Initialize a new complex Gaussian with given parameters.
        pub fn init(position: []const []const T, gamma: []const []const T, momentum: []const []const T, state: usize, bf_spread: []const u32, nstate: usize, allocator: std.mem.Allocator) !@This() {
            if (position.len != gamma.len or position.len != momentum.len or gamma.len != momentum.len) {
                return throw(@This(), "POSITION, GAMMA, AND MOMENTUM ARRAYS MUST HAVE THE SAME LENGTH", .{});
            }

            if (position.len != bf_spread.len) {
                return throw(@This(), "BASIS FUNCTION SPREAD LENGTH MUST MATCH NUMBER OF GAUSSIANS", .{});
            }

            var gaussians = try allocator.alloc(ComplexGaussian(T), position.len); var coefs = try ComplexVector(T).initZero(nstate * gaussians.len, allocator);

            for (0..gaussians.len) |i| {
                gaussians[i] = try ComplexGaussian(T).init(position[i], gamma[i], momentum[i], allocator);
            }

            for (bf_spread, 0..) |contr, i| {
                coefs.ptr(state * gaussians.len + i).* = Complex(T).init(@as(T, @floatFromInt(contr)), 0);
            }

            var mcg = @This(){
                .gaussians = gaussians,
                .coefs = coefs
            };

            var S = try ComplexMatrix(T).init(mcg.coefs.len, mcg.coefs.len, allocator); defer S.deinit(allocator);

            try mcg.overlap(&S);

            try mcg.normalize(S);

            return mcg;
        }

        /// Deinitialize the complex Gaussian array.
        pub fn deinit(self: @This(), allocator: std.mem.Allocator) void {
            for (self.gaussians) |gaussian| gaussian.deinit(allocator);
            allocator.free(self.gaussians);
            self.coefs.deinit(allocator);
        }

        /// Clone the SingleSetOfMCG structure.
        pub fn clone(self: @This(), allocator: std.mem.Allocator) !@This() {
            var gaussians = try allocator.alloc(ComplexGaussian(T), self.gaussians.len);

            for (0..self.gaussians.len) |i| {
                gaussians[i] = try self.gaussians[i].clone(allocator);
            }

            var coefs = try ComplexVector(T).init(self.coefs.len, allocator);

            for (0..self.coefs.len) |i| {
                coefs.ptr(i).* = self.coefs.at(i);
            }

            return @This(){
                .gaussians = gaussians,
                .coefs = coefs
            };
        }

        /// Calculates the derivative of the coefficients with respect to time.
        pub fn coefficientDerivative(self: @This(), dc: *ComplexVector(T), matrix_eom: anytype, allocator: std.mem.Allocator) !void {
            const Sinv = try pseudoInverseHermitianAlloc(T, matrix_eom.S, SINGULARITY_TOLERANCE, allocator); defer Sinv.deinit(allocator);

            var Heff = try ComplexMatrix(T).initZero(self.coefs.len, self.coefs.len, allocator); defer Heff.deinit(allocator);

            for (0..Heff.rows) |i| for (0..Heff.cols) |j| {
                Heff.ptr(i, j).* = matrix_eom.T.at(i, j).add(matrix_eom.V.at(i, j)).sub(matrix_eom.tau.at(i, j).mulbyi());
            };

            const Hc = try mmAlloc(T, Heff, false, self.coefs.asMatrix(), false, allocator); defer Hc.deinit(allocator);

            var dc_matrix = dc.asMatrix();

            try mm(T, &dc_matrix, Sinv, false, Hc, false);

            dc.muls(Complex(T).init(0, -1));
        }

        /// Calculates the derivative of the coefficients with respect to imaginary time.
        pub fn coefficientDerivativeImaginary(self: @This(), dc: *ComplexVector(T), matrix_eom: anytype, allocator: std.mem.Allocator) !void {
            const Sinv = try pseudoInverseHermitianAlloc(T, matrix_eom.S, SINGULARITY_TOLERANCE, allocator); defer Sinv.deinit(allocator);

            var Heff = try ComplexMatrix(T).initZero(self.coefs.len, self.coefs.len, allocator); defer Heff.deinit(allocator);

            for (0..Heff.rows) |i| for (0..Heff.cols) |j| {
                Heff.ptr(i, j).* = matrix_eom.T.at(i, j).add(matrix_eom.V.at(i, j)).add(matrix_eom.tau.at(i, j));
            };

            const Hc = try mmAlloc(T, Heff, false, self.coefs.asMatrix(), false, allocator); defer Hc.deinit(allocator);

            var dc_matrix = dc.asMatrix();

            try mm(T, &dc_matrix, Sinv, false, Hc, false);

            dc.muls(Complex(T).init(-1, 0));
        }

        /// Calculates the expectation value of position or momentum.
        pub fn coordinateExpectation(self: @This(), coord: enum {momentum, position}, S: ComplexMatrix(T), allocator: std.mem.Allocator) !RealVector(T) {
            var coordinate = try RealVector(T).initZero(self.gaussians[0].momentum.len, allocator); var total_norm: T = 0;

            for (0..self.coefs.len / self.gaussians.len) |i| {

                for (0..self.gaussians.len) |j| {

                    const index_1 = i * self.gaussians.len + j;

                    for (0..self.gaussians.len) |k| {

                        const index_2 = i * self.gaussians.len + k;

                        const c_1 = self.coefs.at(index_1);
                        const c_2 = self.coefs.at(index_2);

                        total_norm += c_1.conjugate().mul(c_2).mul(S.at(index_1, index_2)).re;

                        const me = switch (coord) {
                            .momentum => try self.gaussians[j].momentumMatrixElement(self.gaussians[k], allocator),
                            .position => try self.gaussians[j].positionMatrixElement(self.gaussians[k], allocator),
                        }; defer me.deinit(allocator);

                        for (0..coordinate.len) |l| {
                            coordinate.ptr(l).* += c_1.conjugate().mul(c_2).mul(me.at(l)).re;
                        }
                    }
                }
            }

            coordinate.divs(total_norm);

            return coordinate;
        }

        /// Exports the parameters of the Gaussians and coefficients into a single complex vector. The parameters are in order: positions, momentas, gammas for each Gaussian followed by the coefficients for each state and Gaussian.
        pub fn exportParameterVector(self: @This(), params: *ComplexVector(T)) void {
            for (self.gaussians, 0..) |gaussian, i| for (0..gaussian.position.len) |j| {

                const index = 3 * i * gaussian.position.len + j;

                params.ptr(index).* = Complex(T).init(gaussian.position[j], 0);
                params.ptr(gaussian.position.len + index).* = Complex(T).init(gaussian.momentum[j], 0);
                params.ptr(2 * gaussian.position.len + index).* = gaussian.gamma[j];
            };

            for (0..self.coefs.len) |i| {
                params.ptr(3 * self.gaussians[0].position.len * self.gaussians.len + i).* = self.coefs.at(i);
            }
        }

        /// Calculates the Ehrenfest-like gamma derivative.
        pub fn gammaDerivativeEhrenfest(self: @This(), dg: *ComplexVector(T), pot: ElectronicPotential(T), mass: []const T, n_nodes: usize, time: T, fdiff_step: T, allocator: std.mem.Allocator) !void {
            for (self.gaussians, 0..) |gaussian, i| {

                var c_i = try ComplexVector(T).initZero(self.coefs.len / self.gaussians.len, allocator); defer c_i.deinit(allocator);

                for (0..c_i.len) |j| {
                    c_i.ptr(j).* = self.coefs.at(j * self.gaussians.len + i);
                }

                var dg_i = try gaussian.gammaDerivative(pot, c_i, mass, n_nodes, time, fdiff_step, allocator); defer dg_i.deinit(allocator);

                for (0..gaussian.gamma.len) |j| {
                    dg.ptr(gaussian.gamma.len * i + j).* = dg_i.at(j);
                }
            }
        }

        /// Calculates the derivative of the gamma parameter with respect to imaginary time in an Ehrenfest-like manner.
        pub fn gammaDerivativeImaginaryEhrenfest(self: @This(), dg: *ComplexVector(T), pot: ElectronicPotential(T), mass: []const T, n_nodes: usize, time: T, fdiff_step: T, allocator: std.mem.Allocator) !void {
            try self.gammaDerivativeEhrenfest(dg, pot, mass, n_nodes, time, fdiff_step, allocator);

            for (0..dg.len) |i| {
                dg.ptr(i).* = dg.at(i).mulbyi().neg();
            }
        }

        /// Calculates the kinetic matrix.
        pub fn kinetic(self: @This(), K: *ComplexMatrix(T), mass: []const T) !void {
            for (0..K.rows) |i| for (0..K.cols) |j| {

                const ig = i % self.gaussians.len; const is = i / self.gaussians.len;
                const jg = j % self.gaussians.len; const js = j / self.gaussians.len;

                if (is == js) {
                    K.ptr(i, j).* = try self.gaussians[ig].kinetic(self.gaussians[jg], mass);
                }
            };
        }

        /// Calculates the kinetic energy.
        pub fn kineticEnergy(self: @This(), K: ComplexMatrix(T), S: ComplexMatrix(T)) !T {
            return try self.operatorExpectation(K, S);
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
        pub fn momentumDerivativeEhrenfest(self: @This(), dp: *RealVector(T), pot: ElectronicPotential(T), n_nodes: usize, time: T, fdiff_step: T, allocator: std.mem.Allocator) !void {
            for (self.gaussians, 0..) |gaussian, i| {

                var c_i = try ComplexVector(T).initZero(self.coefs.len / self.gaussians.len, allocator); defer c_i.deinit(allocator);

                for (0..c_i.len) |j| {
                    c_i.ptr(j).* = self.coefs.at(j * self.gaussians.len + i);
                }

                var dp_i = try gaussian.momentumDerivative(pot, c_i, n_nodes, time, fdiff_step, allocator); defer dp_i.deinit(allocator);

                for (0..gaussian.momentum.len) |j| {
                    dp.ptr(gaussian.momentum.len * i + j).* = dp_i.at(j);
                }
            }
        }

        /// Calculates the Ehrenfest-like momentum derivative in imaginary time.
        pub fn momentumDerivativeImaginaryEhrenfest(self: @This(), dp: *RealVector(T), mass: []const T) !void {
            try self.positionDerivativeEhrenfest(dp, mass); dp.muls(-1);
        }

        /// Normalizes the coefficient vector in place such that <Psi|Psi> = 1.
        pub fn normalize(self: *@This(), S: ComplexMatrix(T)) !void {
            var norm_sq: T = 0;

            for (0..self.coefs.len) |i| {

                var S_psi_i = Complex(T).init(0, 0);
                
                for (0..self.coefs.len) |j| {
                    S_psi_i = S_psi_i.add(S.at(i, j).mul(self.coefs.at(j)));
                }

                norm_sq += self.coefs.at(i).conjugate().mul(S_psi_i).re;
            }

            self.coefs.divs(Complex(T).init(std.math.sqrt(norm_sq), 0));
        }

        /// Calculate the expectation value of an operator given its matrix representation.
        pub fn operatorExpectation(self: @This(), M: ComplexMatrix(T), S: ComplexMatrix(T)) !T {
            var Mexp: T = 0; var norm: T = 0;

            for (0..self.coefs.len) |i| {

                var row_acc_M = Complex(T).init(0, 0); 
                var row_acc_S = Complex(T).init(0, 0);

                for (0..self.coefs.len) |j| {
                    row_acc_M = row_acc_M.add(M.at(i, j).mul(self.coefs.at(j)));
                    row_acc_S = row_acc_S.add(S.at(i, j).mul(self.coefs.at(j)));
                }

                Mexp += self.coefs.at(i).conjugate().mul(row_acc_M).re;
                norm += self.coefs.at(i).conjugate().mul(row_acc_S).re;
            }

            return Mexp / norm;
        }

        /// Calculates the overlap matrix.
        pub fn overlap(self: @This(), S: *ComplexMatrix(T)) !void {
            for (0..S.rows) |i| for (0..S.cols) |j| {

                const ig = i % self.gaussians.len; const is = i / self.gaussians.len;
                const jg = j % self.gaussians.len; const js = j / self.gaussians.len;

                if (is == js) {
                    S.ptr(i, j).* = try self.gaussians[ig].overlap(self.gaussians[jg]);
                }
            };
        }

        /// Compute the overlap integral matrix where the ket vector is differentiated in time.
        pub fn overlapDiffTime(self: @This(), tau: *ComplexMatrix(T), dq: RealVector(T), dp: RealVector(T), dg: ComplexVector(T)) !void {
            for (0..tau.rows) |i| for (0..tau.cols) |j| {

                const ig = i % self.gaussians.len; const is = i / self.gaussians.len;
                const jg = j % self.gaussians.len; const js = j / self.gaussians.len;

                if (is == js) {

                    const dq_j = dq.slice(jg * self.gaussians[0].position.len, (jg + 1) * self.gaussians[0].position.len);
                    const dp_j = dp.slice(jg * self.gaussians[0].position.len, (jg + 1) * self.gaussians[0].position.len);
                    const dg_j = dg.slice(jg * self.gaussians[0].position.len, (jg + 1) * self.gaussians[0].position.len);

                    tau.ptr(i, j).* = try self.gaussians[ig].overlapDiffTime(self.gaussians[jg], dq_j, dp_j, dg_j);
                }
            };
        }

        /// Calculates the population of a given state. If the coefficients are passed, those are used instead of the internal ones.
        pub fn population(self : @This(), S: ComplexMatrix(T), state: usize, coefs_to_use: ?ComplexVector(T)) !T {
            var pop: T = 0;

            const coefs = coefs_to_use orelse self.coefs;

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
        pub fn positionDerivativeEhrenfest(self: @This(), dq: *RealVector(T), mass: []const T) !void {
            for (self.gaussians, 0..) |gaussian, i| for (0..gaussian.position.len) |j| {
                dq.ptr(gaussian.position.len * i + j).* = gaussian.positionDerivativeIndex(mass, j);
            };
        }

        /// Calculates the Ehrenfest-like position derivative in imaginary time.
        pub fn positionDerivativeImaginaryEhrenfest(self: @This(), dq: *RealVector(T), pot: ElectronicPotential(T), n_nodes: usize, time: T, fdiff_step: T, allocator: std.mem.Allocator) !void {
            try self.momentumDerivativeEhrenfest(dq, pot, n_nodes, time, fdiff_step, allocator);
        }

        /// Calculates the potential matrix for a given potential function.
        pub fn potential(self: @This(), V: *ComplexMatrix(T), pot: ElectronicPotential(T), n_nodes: usize, time: T, allocator: std.mem.Allocator) !void {
            V.zero();

            for (0..self.gaussians.len) |i| for (0..self.gaussians.len) |j| {

                var Vij = try self.gaussians[i].potential(self.gaussians[j], pot, n_nodes, time, allocator); defer Vij.deinit(allocator);

                for (0..Vij.rows) |k| for (0..Vij.cols) |l| {

                    const row = k * self.gaussians.len + i;
                    const col = l * self.gaussians.len + j;

                    V.ptr(row, col).* = V.at(row, col).add(Vij.at(k, l));
                };
            };
        }

        /// Calculates the potential energy.
        pub fn potentialEnergy(self: @This(), V: ComplexMatrix(T), S: ComplexMatrix(T)) !T {
            return try self.operatorExpectation(V, S);
        }

        /// Calculates the total scalar overlap between this wavefunction and another wavefunction.
        pub fn selfOverlap(self: @This(), other: @This(), allocator: std.mem.Allocator) !Complex(T) {
            const n_states_self = self.coefs.len / self.gaussians.len;
            const n_states_other = other.coefs.len / other.gaussians.len;

            if (self.coefs.len / self.gaussians.len != other.coefs.len / other.gaussians.len) {
                return throw(Complex(T), "CANNOT CALCULATE OVERLAP: STATE COUNTS DO NOT MATCH (SELF: {}, OTHER: {})", .{n_states_self, n_states_other});
            }

            var S = try ComplexMatrix(T).initZero(self.gaussians.len, other.gaussians.len, allocator); defer S.deinit(allocator);

            for (0..self.gaussians.len) |i| {
                for (0..other.gaussians.len) |j| {
                    S.ptr(i, j).* = try self.gaussians[i].overlap(other.gaussians[j]);
                }
            }

            var total_overlap = Complex(T).init(0, 0);

            for (0..n_states_self) |state| {
                for (0..self.gaussians.len) |i| {

                    const c1 = self.coefs.at(state * self.gaussians.len + i);

                    if (c1.re == 0 and c1.im == 0) continue;

                    for (0..other.gaussians.len) |j| {

                        const c2 = other.coefs.at(state * other.gaussians.len + j);

                        total_overlap = total_overlap.add(c1.conjugate().mul(c2).mul(S.at(i, j)));
                    }
                }
            }

            return total_overlap;
        }

        /// Returns the transformed coefficients. If to_adiabatic is true, transforms to adiabatic representation; otherwise, to diabatic.
        pub fn transformedCoefs(self: @This(), pot: ElectronicPotential(T), time: T, to_adiabatic: bool, allocator: std.mem.Allocator) !ComplexVector(T) {
            const coefs_after = try ComplexVector(T).initZero(self.coefs.len, allocator);

            var diabatic_potential = try RealMatrix(T).init(self.coefs.len / self.gaussians.len, self.coefs.len / self.gaussians.len, allocator); defer diabatic_potential.deinit(allocator);
            var adiabatic_potential = try RealMatrix(T).init(self.coefs.len / self.gaussians.len, self.coefs.len / self.gaussians.len, allocator); defer adiabatic_potential.deinit(allocator);
            var adiabatic_eigenvectors = try RealMatrix(T).init(self.coefs.len / self.gaussians.len, self.coefs.len / self.gaussians.len, allocator); defer adiabatic_eigenvectors.deinit(allocator);
            var previous_eigenvectors = try RealMatrix(T).init(self.coefs.len / self.gaussians.len, self.coefs.len / self.gaussians.len, allocator); defer previous_eigenvectors.deinit(allocator);

            for (self.gaussians, 0..) |gaussian, i| {

                const position = RealVector(T){.data = gaussian.position, .len = gaussian.position.len};

                const coefs_before_i = try ComplexVector(T).initZero(self.coefs.len / self.gaussians.len, allocator); defer coefs_before_i.deinit(allocator);
                const coefs_after_i  = try ComplexVector(T).initZero(self.coefs.len / self.gaussians.len, allocator); defer  coefs_after_i.deinit(allocator);

                const coefs_before_i_matrix = coefs_before_i.asMatrix();
                var    coefs_after_i_matrix =  coefs_after_i.asMatrix();

                for (0..self.coefs.len / self.gaussians.len) |j| {
                    coefs_before_i.ptr(j).* = self.coefs.at(j * self.gaussians.len + i);
                }

                try pot.evaluateEigensystem(&diabatic_potential, &adiabatic_potential, &adiabatic_eigenvectors, position, time);

                if (i > 0) try fixGauge(T, &adiabatic_eigenvectors, previous_eigenvectors);

                try adiabatic_eigenvectors.copyTo(&previous_eigenvectors);

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
