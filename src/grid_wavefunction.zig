//! Struct containing a wavefunction defined on a grid, along with methods for manipulating it.

const std = @import("std");

const complex_matrix = @import("complex_matrix.zig");
const complex_vector = @import("complex_vector.zig");
const eigensystem_solver = @import("eigenproblem_solver.zig");
const electronic_potential = @import("electronic_potential.zig");
const error_handling = @import("error_handling.zig");
const fourier_transform = @import("fourier_transform.zig");
const grid_generator = @import("grid_generator.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const Complex = std.math.complex.Complex;
const ComplexMatrix = complex_matrix.ComplexMatrix;
const ComplexVector = complex_vector.ComplexVector;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const FourierTransform = fourier_transform.FourierTransform;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const cfftn = fourier_transform.cfftn;
const eigensystemSymmetric = eigensystem_solver.eigensystemSymmetric;
const mm = matrix_multiplication.mm;
const momentumAtRow = grid_generator.momentumAtRow;
const positionAtRow = grid_generator.positionAtRow;
const throw = error_handling.throw;

/// A wavefunction defined on a grid.
pub fn GridWavefunction(comptime T: type) type {
    return struct {
        data: ComplexMatrix(T),
        limits: []const []const T,
        npoint: usize,
        nstate: usize,
        ndim: usize,
        mass: T,

        allocator: std.mem.Allocator,

        /// Allocate a wavefunction on a grid.
        pub fn init(npoint: usize, nstate: usize, ndim: usize, limits: []const []const T, mass: T, allocator: std.mem.Allocator) !@This() {
            if (limits.len != ndim) return throw(@This(), "LIMITS LENGTH MUST BE EQUAL TO NUMBER OF DIMENSIONS", .{});
            for (0..ndim) |i| if (limits[i].len != 2) return throw(@This(), "EACH LIMIT MUST HAVE A LENGTH OF 2", .{});
            if (npoint < 2) return throw(@This(), "NUMBER OF POINTS MUST BE AT LEAST 2", .{});

            return @This(){
                .data = try ComplexMatrix(T).init(std.math.pow(usize, npoint, ndim), nstate, allocator),
                .limits = limits,
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

        /// Calculate the density matrix of the wavefunction. The matrix is allocated and returned.
        pub fn density(self: @This(), potential: ElectronicPotential(T), time: T, adiabatic: bool) !ComplexMatrix(T) {
            var density_matrix = try ComplexMatrix(T).init(self.nstate, self.nstate, self.allocator);

            var diabatic_potential = try RealMatrix(T).init(self.nstate, self.nstate, self.allocator); defer diabatic_potential.deinit();
            var adiabatic_potential = try RealMatrix(T).init(self.nstate, self.nstate, self.allocator); defer adiabatic_potential.deinit();
            var adiabatic_eigenvectors = try RealMatrix(T).init(self.nstate, self.nstate, self.allocator); defer adiabatic_eigenvectors.deinit();

            var position_at_row = try RealVector(T).init(self.ndim, self.allocator); defer position_at_row.deinit();

            var wavefunction_row = try ComplexMatrix(T).init(self.nstate, 1, self.allocator); defer wavefunction_row.deinit();
            
            for (0..self.data.rows) |i| {

                for (0..self.nstate) |j| wavefunction_row.ptr(j, 0).* = self.data.at(i, j);

                if (adiabatic) {

                    positionAtRow(T, &position_at_row, i, self.ndim, self.npoint, self.limits);

                    try potential.evaluateEigensystem(&diabatic_potential, &adiabatic_potential, &adiabatic_eigenvectors, position_at_row, time);

                    try mm(T, &wavefunction_row, adiabatic_eigenvectors, true, self.data.row(i).asMatrix(), false);
                }

                for (0..self.nstate) |j| for (0..self.nstate) |k| {
                    density_matrix.ptr(j, k).* = density_matrix.at(j, k).add(wavefunction_row.at(j, 0).mul(wavefunction_row.at(k, 0).conjugate()));
                };
            }

            density_matrix.muls(Complex(T).init(self.getIntegrationElement(), 0));

            return density_matrix;
        }

        /// Get the size of the elemet of the grid.
        pub fn getIntegrationElement(self: @This()) T {
            var dr: T = 1;

            for (0..self.ndim) |i| {
                dr *= (self.limits[i][1] - self.limits[i][0]) / @as(T, @floatFromInt(self.npoint - 1));
            }

            return dr;
        }

        /// Get the shape of the wavefunction as an array of usize.
        pub fn getShape(self: @This()) ![]usize {
            const shape = try self.allocator.alloc(usize, self.ndim);

            for (0..self.ndim) |i| shape[i] = self.npoint;

            return shape;
        }

        /// Initialize the position of the wavefunction as a Gaussian wavepacket.
        pub fn initialGaussian(self: *@This(), position: []const T, momentum: []const T, state: usize, gamma: T) !void {
            if (position.len != self.ndim) return throw(void, "POSITION LENGTH MUST BE EQUAL TO NUMBER OF DIMENSIONS", .{});
            if (momentum.len != self.ndim) return throw(void, "MOMENTUM LENGTH MUST BE EQUAL TO NUMBER OF DIMENSIONS", .{});
            if (state >= self.nstate) return throw(void, "STATE INDEX OUT OF BOUNDS", .{});

            var position_at_row = try RealVector(T).init(self.ndim, self.allocator); defer position_at_row.deinit();

            for (0..self.data.rows) |i| {

                positionAtRow(T, &position_at_row, i, self.ndim, self.npoint, self.limits);

                for (0..self.nstate) |j| if (j == state) {

                    var exp = Complex(T).init(0, 0);

                    for (0..self.ndim) |k| {

                        exp.re -= 0.5 * gamma * (position_at_row.at(k) - position[k]) * (position_at_row.at(k) - position[k]);

                        exp.im += momentum[k] * (position_at_row.at(k) - position[k]);
                    }

                    self.data.ptr(i, j).* = std.math.complex.exp(exp);
                };
            }

            self.normalize();
        }

        /// Kinetic energy operator in momentum space.
        pub fn kineticEnergy(self: @This(), temporary_column: *ComplexVector(T)) !T {
            const shape = try self.getShape(); defer self.allocator.free(shape);

            var momentum_at_row = try RealVector(T).init(self.ndim, self.allocator); defer momentum_at_row.deinit();

            var kinetic_energy: T = 0;

            for (0..self.nstate) |i| {

                for (0..self.data.rows) |j| temporary_column.ptr(j).* = self.data.at(j, i);

                try cfftn(T, temporary_column, shape, -1);

                for (0..self.data.rows) |j| {

                    momentumAtRow(T, &momentum_at_row, j, self.ndim, self.npoint, self.limits);
                    
                    var k_sum_squared: T = 0;

                    for (0..self.ndim) |k| {
                        k_sum_squared += momentum_at_row.at(k) * momentum_at_row.at(k);
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
            const shape = try self.getShape(); defer self.allocator.free(shape);

            var momentum_at_row = try RealVector(T).init(self.ndim, self.allocator); defer momentum_at_row.deinit();

            var momentum = try RealVector(T).initZero(self.ndim, self.allocator);

            for (0..self.ndim) |i| for (0..self.nstate) |j| {

                for (0..self.data.rows) |k| temporary_column.ptr(k).* = self.data.at(k, j);

                try cfftn(f64, temporary_column, shape, -1);

                for (0..self.data.rows) |k| {

                    momentumAtRow(T, &momentum_at_row, k, self.ndim, self.npoint, self.limits);

                    temporary_column.ptr(k).* = temporary_column.at(k).mul(Complex(T).init(momentum_at_row.at(i), 0));
                }

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
            var value = Complex(T).init(0, 0);

            for (0..self.data.rows) |i| for (0..self.nstate) |j| for (0..self.nstate) |k| {
                value = value.add(self.data.at(i, j).conjugate().mul(other.data.at(i, k)));
            };

            return value.mul(Complex(T).init(self.getIntegrationElement(), 0));
        }

        /// Calculate the position of the wavefunction. The resulting vector is allocated inside the function and returned.
        pub fn positionMean(self: @This()) !RealVector(T) {
            var position_at_row = try RealVector(T).init(self.ndim, self.allocator); defer position_at_row.deinit();

            var position = try RealVector(T).initZero(self.ndim, self.allocator);

            for (0..self.data.rows) |i| {

                positionAtRow(T, &position_at_row, i, self.ndim, self.npoint, self.limits);

                for (0..self.ndim) |j| for (0..self.nstate) |k| {
                    position.ptr(j).* += self.data.at(i, k).conjugate().mul(Complex(T).init(position_at_row.at(j), 0)).mul(self.data.at(i, k)).re;
                };
            }

            position.muls(self.getIntegrationElement());

            return position;
        }

        /// Calculate potential energy of the wavefunction. The temporary matrix for storing the potential matrix is allocated inside the function.
        pub fn potentialEnergy(self: *@This(), potential: ElectronicPotential(T), time: T) !T {
            var position_at_row = try RealVector(T).init(self.ndim, self.allocator); defer position_at_row.deinit();

            var diabatic_potential_matrix = try RealMatrix(T).init(self.nstate, self.nstate, self.allocator); defer diabatic_potential_matrix.deinit();

            var potential_energy: T = 0;

            for (0..self.data.rows) |i| {

                positionAtRow(T, &position_at_row, i, self.ndim, self.npoint, self.limits);
                potential.evaluateDiabatic(&diabatic_potential_matrix, position_at_row, time);

                for (0..self.nstate) |j| for (0..self.nstate) |k| {
                    potential_energy += self.data.at(i, j).conjugate().mul(Complex(T).init(diabatic_potential_matrix.at(j, k), 0).mul(self.data.at(i, k))).re;
                };
            }

            return potential_energy * self.getIntegrationElement();
        }

        /// Propagate the wavefunction in time using the split-operator method.
        pub fn propagate(self: *@This(), potential: ElectronicPotential(T), time: T, time_step: T, imaginary: bool, temporary_column: *ComplexVector(T)) !void {
            const unit = Complex(T).init(if (imaginary) 1 else 0, if (imaginary) 0 else 1);

            try propagateHalfPosition(self, potential, time, time_step, unit, temporary_column);

            try propagateFullMomentum(self, time_step, unit, temporary_column);

            try propagateHalfPosition(self, potential, time, time_step, unit, temporary_column);

            if (imaginary) self.normalize();
        }

        /// Propagate the wavefunction full time step in momentum space.
        pub fn propagateFullMomentum(self: *@This(), time_step: T, unit: Complex(T), temporary_column: *ComplexVector(T)) !void {
            const shape = try self.getShape(); defer self.allocator.free(shape);

            var momentum_at_row = try RealVector(T).init(self.ndim, self.allocator); defer momentum_at_row.deinit();

            var K = try ComplexMatrix(T).init(self.nstate, self.nstate, self.allocator); defer K.deinit();

            for (0..self.nstate) |j| {

                for (0..self.data.rows) |i| temporary_column.ptr(i).* = self.data.at(i, j);

                try cfftn(f64, temporary_column, shape, -1);

                for (0..self.data.rows) |i| self.data.ptr(i, j).* = temporary_column.at(i);
            }

            for (0..self.data.rows) |i| {

                momentumAtRow(T, &momentum_at_row, i, self.ndim, self.npoint, self.limits);
                getMomentumPropagator(T, &K, momentum_at_row, self.mass, self.nstate, time_step, unit);

                for (0..self.nstate) |j| temporary_column.ptr(j).* = Complex(T).init(0, 0);

                for (0..self.nstate) |j| for (0..self.nstate) |k| {
                    temporary_column.ptr(j).* = temporary_column.at(j).add(K.at(j, k).mul(self.data.at(i, k)));
                };

                for (0..self.nstate) |j| self.data.ptr(i, j).* = temporary_column.at(j);
            }

            for (0..self.nstate) |j| {

                for (0..self.data.rows) |i| temporary_column.ptr(i).* = self.data.at(i, j);

                try cfftn(f64, temporary_column, shape, 1);

                for (0..self.data.rows) |i| self.data.ptr(i, j).* = temporary_column.at(i);
            }
        }

        /// Propagate the wavefunction half a time step in position space.
        pub fn propagateHalfPosition(self: *@This(), potential: ElectronicPotential(T), time: T, time_step: T, unit: Complex(T), temporary_column: *ComplexVector(T)) !void {
            var diabatic_potential = try RealMatrix(T).init(self.nstate, self.nstate, self.allocator); defer diabatic_potential.deinit();
            var adiabatic_potential = try RealMatrix(T).init(self.nstate, self.nstate, self.allocator); defer adiabatic_potential.deinit();
            var adiabatic_eigenvectors = try RealMatrix(T).init(self.nstate, self.nstate, self.allocator); defer adiabatic_eigenvectors.deinit();

            var position_at_row = try RealVector(T).init(self.ndim, self.allocator); defer position_at_row.deinit();

            var R = try ComplexMatrix(T).init(self.nstate, self.nstate, self.allocator); defer R.deinit();

            var mm_temporary = try ComplexMatrix(T).init(self.nstate, self.nstate, self.allocator); defer mm_temporary.deinit();

            for (0..self.data.rows) |i| {

                positionAtRow(T, &position_at_row, i, self.ndim, self.npoint, self.limits);

                try potential.evaluateEigensystem(&diabatic_potential, &adiabatic_potential, &adiabatic_eigenvectors, position_at_row, time);

                try getPositionPropagator(T, &R, adiabatic_potential, adiabatic_eigenvectors, time_step, unit, &mm_temporary);

                for (0..self.nstate) |j| {

                    temporary_column.ptr(j).* = Complex(T).init(0, 0);

                    for (0..self.nstate) |k| {
                        temporary_column.ptr(j).* = temporary_column.at(j).add(R.at(j, k).mul(self.data.at(i, k)));
                    }
                }

                for (0..self.nstate) |j| self.data.ptr(i, j).* = temporary_column.at(j);
            }
        }

        /// Perform the transformation from diabatic to adiabatic representation or vice versa.
        pub fn transformRepresentation(self: *@This(), potential: ElectronicPotential(T), time: T, to_adiabatic: bool) !void {
            var diabatic_potential = try RealMatrix(T).init(self.nstate, self.nstate, self.allocator); defer diabatic_potential.deinit();
            var adiabatic_potential = try RealMatrix(T).init(self.nstate, self.nstate, self.allocator); defer adiabatic_potential.deinit();
            var adiabatic_eigenvectors = try RealMatrix(T).init(self.nstate, self.nstate, self.allocator); defer adiabatic_eigenvectors.deinit();

            var position_at_row = try RealVector(T).init(self.ndim, self.allocator); defer position_at_row.deinit();

            var mm_temporary = try ComplexMatrix(T).init(self.nstate, 1, self.allocator); defer mm_temporary.deinit();

            for (0..self.data.rows) |i| {

                positionAtRow(T, &position_at_row, i, self.ndim, self.npoint, self.limits);

                try potential.evaluateEigensystem(&diabatic_potential, &adiabatic_potential, &adiabatic_eigenvectors, position_at_row, time);

                if (to_adiabatic) {
                    try mm(T, &mm_temporary, adiabatic_eigenvectors, true, self.data.row(i).asMatrix(), false);
                } else {
                    try mm(T, &mm_temporary, adiabatic_eigenvectors, false, self.data.row(i).asMatrix(), false);
                }

                for (0..self.nstate) |j| self.data.ptr(i, j).* = mm_temporary.at(j, 0);
            }
        }
    };
}

/// Returns a wavefunction propagator for a full step in momentum space.
pub fn getMomentumPropagator(comptime T: type, K: *ComplexMatrix(T), momentum: RealVector(T), mass: T, nstate: usize, time_step: T, unit: Complex(T)) void {
    K.zero();

    for (0..nstate) |j| {

        for (0..momentum.len) |k| {
            K.ptr(j, j).* = K.at(j, j).add(Complex(T).init(momentum.at(k) * momentum.at(k), 0));
        }

        K.ptr(j, j).* = std.math.complex.exp(K.at(j, j).mul(Complex(T).init(-0.5 * time_step / mass, 0)).mul(unit));
    }
}

/// Returns a wavefunction propagator for a half step in position space.
pub fn getPositionPropagator(comptime T: type, R: *ComplexMatrix(T), adiabatic: RealMatrix(T), eigenvectors: RealMatrix(T), time_step: T, unit: Complex(T), mm_temporary: *ComplexMatrix(T)) !void {
    R.zero();

    for (0..R.rows) |j| {
        R.ptr(j, j).* = std.math.complex.exp(Complex(T).init(adiabatic.at(j, j), 0).mul(Complex(T).init(-0.5 * time_step, 0)).mul(unit));
    }

    try mm(T, mm_temporary, eigenvectors, false, R.*, false);
    try mm(T, R, mm_temporary.*, false, eigenvectors, true);
}
