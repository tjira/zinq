//! Struct containing a wavefunction defined on a grid, along with methods for manipulating it.

const std = @import("std");

const complex_matrix = @import("complex_matrix.zig");
const complex_vector = @import("complex_vector.zig");
const electronic_potential = @import("electronic_potential.zig");
const fourier_transform = @import("fourier_transform.zig");
const linear_algebra = @import("linear_algebra.zig");
const real_matrix = @import("real_matrix.zig");
const grid_wavefunction_propagator = @import("grid_wavefunction_propagator.zig");
const real_vector = @import("real_vector.zig");

const Complex = std.math.complex.Complex;
const GridWavefunctionPropagator = grid_wavefunction_propagator.GridWavefunctionPropagator;
const ComplexMatrix = complex_matrix.ComplexMatrix;
const ComplexVector = complex_vector.ComplexVector;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const FourierTransform = fourier_transform.FourierTransform;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const mmRealAlloc = linear_algebra.mmRealAlloc;
const mmReal = linear_algebra.mmReal;
const mmComplexAlloc = linear_algebra.mmComplexAlloc;
const eigensystemSymmetric = linear_algebra.eigensystemSymmetric;
const cfftn = fourier_transform.cfftn;

/// A wavefunction defined on a grid.
pub fn GridWavefunction(comptime T: type) type {
    return struct {
        data: ComplexMatrix(T),
        npoint: usize,
        nstate: usize,
        ndim: usize,
        mass: T,

        position_grid_pointer: ?*const RealMatrix(T) = null,
        momentum_grid_pointer: ?*const RealMatrix(T) = null,

        allocator: std.mem.Allocator,

        /// Allocate a wavefunction on a grid.
        pub fn init(npoint: usize, nstate: usize, ndim: usize, mass: T, allocator: std.mem.Allocator) !@This() {
            return @This(){
                .data = try ComplexMatrix(T).init(std.math.pow(usize, npoint, ndim), nstate, allocator),
                .npoint = npoint,
                .nstate = nstate,
                .ndim = ndim,
                .mass = mass,
                .allocator = allocator
            };
        }

        /// Free the memory allocated for the wavefunction.
        pub fn deinit(self: @This()) void {
            self.data.deinit();
        }

        /// Attach position and momentum grids to the wavefunction.
        pub fn attachGrids(self: *@This(), position_grid: *const RealMatrix(T), momentum_grid: *const RealMatrix(T)) void {
            self.position_grid_pointer = position_grid;
            self.momentum_grid_pointer = momentum_grid;
        }

        /// Calculate the density matrix of the wavefunction. The matrix is allocated and returned.
        pub fn density(self: @This()) !ComplexMatrix(T) {
            var density_matrix = try ComplexMatrix(T).init(self.nstate, self.nstate, self.allocator);

            for (0..self.nstate) |i| for (0..self.nstate) |j| for (0..self.npoint) |k| {
                density_matrix.ptr(i, j).* = density_matrix.at(i, j).add(self.data.at(k, i).mul(self.data.at(k, j).conjugate()));
            };

            density_matrix.muls(Complex(T).init(self.getIntegrationElement(), 0));

            return density_matrix;
        }

        /// Get the size of the elemet of the grid.
        pub fn getIntegrationElement(self: @This()) T {
            const dr = self.position_grid_pointer.?.at(1, self.ndim - 1) - self.position_grid_pointer.?.at(0, self.ndim - 1);

            return std.math.pow(T, dr, @as(T, @floatFromInt(self.ndim))) ;
        }

        /// Initialize the position of the wavefunction as a Gaussian wavepacket.
        pub fn initialGaussian(self: *@This(), position: []const T, momentum: []const T, state: usize, gamma: T) void {
            for (0..self.npoint) |i| for (0..self.nstate) |j| if (j == state) {

                var exp = Complex(T).init(0, 0);

                for (0..self.ndim) |k| {

                    exp.re -= 0.5 * gamma * (self.position_grid_pointer.?.at(i, k) - position[k]) * (self.position_grid_pointer.?.at(i, k) - position[k]);

                    exp.im += momentum[k] * (self.position_grid_pointer.?.at(i, k) - position[k]);
                }

                self.data.ptr(i, j).* = std.math.complex.exp(exp);
            };

            self.normalize();
        }

        /// Kinetic energy operator in momentum space.
        pub fn kineticEnergy(self: @This(), temporary_column: *ComplexVector(T)) !T {
            const shape = try self.allocator.alloc(usize, self.ndim); defer self.allocator.free(shape);

            @memset(shape, self.npoint);

            var kinetic_energy: T = 0;

            for (0..self.nstate) |i| {

                for (0..self.data.rows) |j| temporary_column.ptr(j).* = self.data.at(j, i);

                try cfftn(T, temporary_column, shape, -1);

                for (0..self.data.rows) |j| {
                    
                    var k_sum_squared: T = 0;

                    for (0..self.ndim) |k| {
                        k_sum_squared += self.momentum_grid_pointer.?.at(j, k) * self.momentum_grid_pointer.?.at(j, k);
                    }

                    temporary_column.ptr(j).* = temporary_column.at(j).mul(Complex(T).init(k_sum_squared, 0));
                }

                try cfftn(T, temporary_column, shape, 1);

                for (0..self.data.rows) |j| {
                    kinetic_energy += temporary_column.at(j).mul(self.data.at(j, i).conjugate()).re;
                }
            }

            return 0.5 * kinetic_energy * self.getIntegrationElement() / self.mass;
        }

        /// Calculate the momentum of the wavefunction. The resulting vector is allocated inside the function and returned.
        pub fn momentumMean(self: @This(), temporary_column: *ComplexVector(T)) !RealVector(T) {
            const shape = try self.allocator.alloc(usize, self.ndim); defer self.allocator.free(shape);

            @memset(shape, self.npoint);

            var momentum = try RealVector(T).initZero(self.ndim, self.allocator);

            for (0..self.ndim) |i| for (0..self.nstate) |j| {

                for (0..self.data.rows) |k| temporary_column.ptr(k).* = self.data.at(k, j);

                try cfftn(f64, temporary_column, shape, -1);

                for (0..self.data.rows) |k| temporary_column.ptr(k).* = temporary_column.at(k).mul(Complex(T).init(self.momentum_grid_pointer.?.at(k, i), 0));

                try cfftn(f64, temporary_column, shape, 1);

                for (0..self.data.rows) |k| momentum.ptr(i).* += temporary_column.at(k).mul(self.data.at(k, j).conjugate()).re;
            };

            momentum.muls(self.getIntegrationElement());

            return momentum;
        }

        /// Normalize the wavefunction.
        pub fn normalize(self: *@This()) void {
            self.data.divs(Complex(T).init(std.math.sqrt(self.overlap(self.*).magnitude()), 0));
        }

        /// Calculate the overlap with another wavefunction.
        pub fn overlap(self: @This(), other: @This()) Complex(T) {
            var value = std.math.Complex(T).init(0, 0);

            for (0..self.npoint) |i| for (0..self.nstate) |j| for (0..self.nstate) |k| {
                value = value.add(self.data.at(i, j).conjugate().mul(other.data.at(i, k)));
            };

            return value.mul(Complex(T).init(self.getIntegrationElement(), 0));
        }

        /// Calculate the position of the wavefunction. The resulting vector is allocated inside the function and returned.
        pub fn positionMean(self: @This()) !RealVector(T) {
            var position = try RealVector(T).initZero(self.ndim, self.allocator);

            for (0..self.ndim) |i| for (0..self.nstate) |j| for (0..self.data.rows) |k| {
                position.ptr(i).* += self.data.at(k, j).conjugate().mul(Complex(T).init(self.position_grid_pointer.?.at(k, i), 0)).mul(self.data.at(k, j)).re;
            };

            position.muls(self.getIntegrationElement());

            return position;
        }

        /// Calculate potential energy of the wavefunction. The temporary matrix for storing the potential matrix is allocated inside the function.
        pub fn potentialEnergy(self: *@This(), potential: ElectronicPotential(T), time: T) !T {
            var diabatic_potential_matrix = try RealMatrix(T).init(self.nstate, self.nstate, self.allocator); defer diabatic_potential_matrix.deinit();

            var potential_energy: T = 0;

            for (0..self.data.rows) |i| {

                potential.evaluateDiabatic(&diabatic_potential_matrix, self.position_grid_pointer.?.row(i), time);

                for (0..self.nstate) |j| for (0..self.nstate) |k| {
                    potential_energy += self.data.at(i, j).conjugate().mul(Complex(T).init(diabatic_potential_matrix.at(j, k), 0).mul(self.data.at(i, k))).re;
                };
            }

            return potential_energy * self.getIntegrationElement();
        }

        /// Propagate the wavefunction in time using the split-operator method.
        pub fn propagate(self: *@This(), propagator: GridWavefunctionPropagator(T), imaginary: bool, temporary_column: *ComplexVector(T)) !void {
            const shape = try self.allocator.alloc(usize, self.ndim); defer self.allocator.free(shape);

            @memset(shape, self.npoint);

            for (0..self.data.rows) |i| {
                for (0..self.nstate) |j| temporary_column.ptr(j).* = std.math.Complex(T).init(0, 0);
                for (0..self.nstate) |j| for (0..self.nstate) |k| {
                    temporary_column.ptr(j).* = temporary_column.at(j).add(propagator.position_propagators[i].at(j, k).mul(self.data.at(i, k)));
                };
                for (0..self.nstate) |j| self.data.ptr(i, j).* = temporary_column.at(j);
            }

            for (0..self.nstate) |j| {
                for (0..self.data.rows) |i| temporary_column.data[i] = self.data.at(i, j);
                try cfftn(f64, temporary_column, shape, -1);
                for (0..self.data.rows) |i| self.data.ptr(i, j).* = temporary_column.at(i);
            }

            for (0..self.data.rows) |i| {
                for (0..self.nstate) |j| temporary_column.ptr(j).* = std.math.Complex(T).init(0, 0);
                for (0..self.nstate) |j| for (0..self.nstate) |k| {
                    temporary_column.ptr(j).* = temporary_column.at(j).add(propagator.momentum_propagators[i].at(j, k).mul(self.data.at(i, k)));
                };
                for (0..self.nstate) |j| self.data.ptr(i, j).* = temporary_column.at(j);
            }

            for (0..self.nstate) |j| {
                for (0..self.data.rows) |i| temporary_column.data[i] = self.data.at(i, j);
                try cfftn(f64, temporary_column, shape, 1);
                for (0..self.data.rows) |i| self.data.ptr(i, j).* = temporary_column.at(i);
            }

            for (0..self.data.rows) |i| {
                for (0..self.nstate) |j| temporary_column.ptr(j).* = std.math.Complex(T).init(0, 0);
                for (0..self.nstate) |j| for (0..self.nstate) |k| {
                    temporary_column.ptr(j).* = temporary_column.at(j).add(propagator.position_propagators[i].at(j, k).mul(self.data.at(i, k)));
                };
                for (0..self.nstate) |j| self.data.ptr(i, j).* = temporary_column.at(j);
            }

            if (imaginary) self.normalize();
        }
    };
}
