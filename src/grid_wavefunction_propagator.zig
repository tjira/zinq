//! A propagator for grid-based wavefunctions in quantum mechanics.

const std = @import("std");

const complex_vector = @import("complex_vector.zig");
const linear_algebra = @import("linear_algebra.zig");
const complex_matrix = @import("complex_matrix.zig");
const real_matrix = @import("real_matrix.zig");
const grid_wavefunction = @import("grid_wavefunction.zig");
const electronic_potential = @import("electronic_potential.zig");

const RealMatrix = real_matrix.RealMatrix;
const Complex = std.math.complex.Complex;
const ComplexVector = complex_vector.ComplexVector;
const ComplexMatrix = complex_matrix.ComplexMatrix;
const GridWavefunction = grid_wavefunction.GridWavefunction;
const eigensystemSymmetric = linear_algebra.eigensystemSymmetric;
const mmComplex = linear_algebra.mmReal;
const ElectronicPotential = electronic_potential.ElectronicPotential;

/// Propagator for grid-based wavefunctions.
pub fn GridWavefunctionPropagator(comptime T: type) type {
    return struct {
        position_propagators: []ComplexMatrix(T),
        momentum_propagators: []ComplexMatrix(T),
        allocator: std.mem.Allocator,

        pub fn init(npoint: usize, ndim: usize, nstate: usize, allocator: std.mem.Allocator) !@This() {
            var propagator = @This(){
                .position_propagators = try allocator.alloc(ComplexMatrix(T), std.math.pow(usize, npoint, ndim)),
                .momentum_propagators = try allocator.alloc(ComplexMatrix(T), std.math.pow(usize, npoint, ndim)),
                .allocator = allocator
            };

            for (0..std.math.pow(usize, npoint, ndim)) |i| {
                propagator.position_propagators[i] = try ComplexMatrix(T).initZero(nstate, nstate, allocator);
                propagator.momentum_propagators[i] = try ComplexMatrix(T).initZero(nstate, nstate, allocator);
            }

            return propagator;
        }

        /// Free the memory allocated for the propagator.
        pub fn deinit(self: *@This()) void {
            for (self.position_propagators) |element| element.deinit();
            for (self.momentum_propagators) |element| element.deinit();

            self.allocator.free(self.position_propagators);
            self.allocator.free(self.momentum_propagators);
        }

        pub fn generate(self: *@This(), wavefunction: GridWavefunction(T), potential: ElectronicPotential(T), time: T, time_step: T, imaginary: bool) !void {
            const unit = std.math.Complex(T).init(if (imaginary) 1 else 0, if (imaginary) 0 else 1);

            var diabatic_potential = try RealMatrix(T).init(wavefunction.nstate, wavefunction.nstate, self.allocator); defer diabatic_potential.deinit();
            var adiabatic_potential = try RealMatrix(T).init(wavefunction.nstate, wavefunction.nstate, self.allocator); defer adiabatic_potential.deinit();
            var adiabatic_eigenvectors = try RealMatrix(T).init(wavefunction.nstate, wavefunction.nstate, self.allocator); defer adiabatic_eigenvectors.deinit();

            var temporary_multiplication_matrix = try ComplexMatrix(T).initZero(wavefunction.nstate, wavefunction.nstate, self.allocator); defer temporary_multiplication_matrix.deinit();

            for (0..wavefunction.data.rows) |i| {

                potential.evaluateDiabatic(&diabatic_potential, wavefunction.position_grid_pointer.?.row(i), time);
                try eigensystemSymmetric(T, &adiabatic_potential, &adiabatic_eigenvectors, diabatic_potential);

                for (0..wavefunction.nstate) |j| {

                    for (0..wavefunction.ndim) |k| {

                        const momentum = wavefunction.momentum_grid_pointer.?.at(i, k);

                        self.momentum_propagators[i].ptr(j, j).* = self.momentum_propagators[i].at(j, j).add(Complex(T).init(momentum * momentum, 0));
                    }

                    self.momentum_propagators[i].ptr(j, j).* = std.math.complex.exp(self.momentum_propagators[i].at(j, j).mul(Complex(T).init(-0.5 * time_step / wavefunction.mass, 0)).mul(unit));

                    self.position_propagators[i].ptr(j, j).* = std.math.complex.exp(Complex(T).init(adiabatic_potential.at(j, j), 0).mul(Complex(T).init(-0.5 * time_step, 0)).mul(unit));
                }

                for (0..adiabatic_eigenvectors.rows) |k| for (0..self.position_propagators[i].cols) |l| {

                    var sum = Complex(T).init(0, 0);

                    for (0..adiabatic_eigenvectors.cols) |m| {
                        sum = sum.add(Complex(T).init(adiabatic_eigenvectors.at(k, m), 0).mul(self.position_propagators[i].at(m, l)));
                    }

                    temporary_multiplication_matrix.ptr(k, l).* = sum;
                };

                for (0..temporary_multiplication_matrix.rows) |k| for (0..adiabatic_eigenvectors.rows) |l| {

                    var sum = Complex(T).init(0, 0);

                    for (0..temporary_multiplication_matrix.cols) |m| {
                        sum = sum.add(temporary_multiplication_matrix.at(k, m).mul(Complex(T).init(adiabatic_eigenvectors.at(l, m), 0)));
                    }

                    self.position_propagators[i].ptr(k, l).* = sum;
                };
            }
        }
    };
}
