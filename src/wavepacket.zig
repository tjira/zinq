//! Wavefunction and Hamiltonian representations for grid-based quantum wavepacket dynamics.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;

const FftPlan = @import("fourier_transform.zig").FftPlan;
const Matrix = @import("tensor.zig").Matrix;
const Potential = @import("potential.zig").Potential;
const Vector = @import("tensor.zig").Vector;

const eighBatch = @import("linear_algebra.zig").eighBatch;

/// Initial parameters of the wavepacket including position, momentum, and Gaussian width.
pub const InitialConditions = struct {
    position: []const f64,
    momentum: []const f64,
    gamma: []const f64,

    state: u32 = 0,
    adiabatic: bool = false,
};

/// Generates a multidimensional discrete coordinate and momentum space grid.
pub fn Grid(comptime T: type) type {
    return struct {
        r: Matrix(T),
        k: Matrix(T),

        dr: T,
        dk: T,

        /// Allocates and initializes grid coordinates and momentum vectors.
        pub fn init(bounds: []const [2]T, npoint: u32, gpa: Allocator) !@This() {
            const ncol = std.math.pow(usize, npoint, bounds.len);

            var r = try Matrix(T).init(ncol, bounds.len, gpa);
            errdefer r.deinit(gpa);

            var k = try Matrix(T).init(ncol, bounds.len, gpa);
            errdefer k.deinit(gpa);

            var dr: T = 1;

            for (0..bounds.len) |i| {
                const min = bounds[i][0];
                const max = bounds[i][1];

                dr *= (max - min) / @as(T, @floatFromInt(npoint));
            }

            const dk = dr / @as(T, @floatFromInt(ncol));

            for (0..ncol) |i| {
                var temp = i;

                for (0..bounds.len) |l| {
                    const j = bounds.len - l - 1;

                    const n = @as(T, @floatFromInt(temp % npoint));

                    const min = bounds[j][0];
                    const max = bounds[j][1];

                    r.ptr(i, j).* = min + n * (max - min) / @as(T, @floatFromInt(npoint));

                    const dki = 2 * std.math.pi / (max - min);

                    k.ptr(i, j).* = (if (temp % npoint < npoint / 2) n else n - @as(T, @floatFromInt(npoint))) * dki;

                    temp /= npoint;
                }
            }

            return .{ .r = r, .k = k, .dr = dr, .dk = dk };
        }

        /// Deallocates coordinate and momentum space grid matrices.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.r.deinit(gpa);
            self.k.deinit(gpa);
        }
    };
}

/// Representation of kinetic energy and potential energy operators on the grid.
pub fn Hamiltonian(comptime T: type) type {
    return struct {
        V: Matrix(T),
        W: Matrix(T),
        U: Matrix(T),
        K: Vector(T),

        /// Allocates and computes kinetic and potential operator matrix elements.
        pub fn init(grid: Grid(T), pot: Potential(T), m: T, gpa: Allocator) !@This() {
            var V = try Matrix(T).init(grid.r.nrow(), pot.nstate() * pot.nstate(), gpa);
            errdefer V.deinit(gpa);

            var U = try Matrix(T).init(grid.r.nrow(), pot.nstate() * pot.nstate(), gpa);
            errdefer U.deinit(gpa);

            var W = try Matrix(T).init(grid.r.nrow(), pot.nstate(), gpa);
            errdefer W.deinit(gpa);

            var K = try Vector(T).initZero(grid.r.nrow(), gpa);
            errdefer K.deinit(gpa);

            for (0..grid.r.nrow()) |i| {
                var sum: T = 0;

                for (0..grid.r.ncol()) |j| {
                    const kij = grid.k.at(i, j);

                    sum += 0.5 * kij * kij / m;
                }

                K.ptr(i).* = sum;
            }

            var ham = @This(){ .V = V, .W = W, .U = U, .K = K };

            try ham.update(grid, pot, 0, gpa);

            return ham;
        }

        /// Deallocates Hamiltonian operator matrices.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.V.deinit(gpa);
            self.W.deinit(gpa);
            self.U.deinit(gpa);
            self.K.deinit(gpa);
        }

        /// Updates potential energy values and diagonalizes to get adiabatic states.
        pub fn update(self: *@This(), grid: Grid(T), pot: Potential(T), t: T, gpa: Allocator) !void {
            var U_prev = if (pot.isTd() and t > 0) try self.U.clone(gpa) else null;
            defer if (U_prev) |*u| u.deinit(gpa);

            pot.evalBatch(T, &self.V, grid.r, t);

            try eighBatch(T, &self.W, &self.U, self.V);

            if (grid.r.ncol() == 1) for (1..grid.r.nrow()) |i| for (0..pot.nstate()) |j| {
                var overlap: T = 0;

                for (0..pot.nstate()) |k| {
                    overlap += self.U.at(i, k * pot.nstate() + j) * self.U.at(i - 1, k * pot.nstate() + j);
                }

                if (overlap < 0) for (0..pot.nstate()) |k| {
                    self.U.ptr(i, k * pot.nstate() + j).* = -self.U.at(i, k * pot.nstate() + j);
                };
            };

            if (U_prev) |prev| for (0..pot.nstate()) |j| {
                var total_overlap: T = 0;

                for (0..grid.r.nrow()) |i| for (0..pot.nstate()) |k| {
                    total_overlap += self.U.at(i, k * pot.nstate() + j) * prev.at(i, k * pot.nstate() + j);
                };

                if (total_overlap < 0) for (0..grid.r.nrow()) |i| for (0..pot.nstate()) |k| {
                    self.U.ptr(i, k * pot.nstate() + j).* = -self.U.at(i, k * pot.nstate() + j);
                };
            };
        }
    };
}

/// Representation of a multi-state wavepacket and its Fourier transform plans.
pub fn Wavefunction(comptime T: type) type {
    return struct {
        W: Matrix(Complex(T)),

        ffft: FftPlan(Complex(T)),
        ifft: FftPlan(Complex(T)),

        /// Allocates wavefunction components and plans forward and backward FFTs.
        pub fn init(ndim: usize, nstate: usize, npoint: usize, plan_mode: u32, gpa: Allocator) !@This() {
            var W = try Matrix(Complex(T)).init(nstate, std.math.pow(usize, npoint, ndim), gpa);
            errdefer W.deinit(gpa);

            const shape = try gpa.alloc(i32, ndim);
            defer gpa.free(shape);

            for (0..shape.len) |i| {
                shape[i] = @as(i32, @intCast(npoint));
            }

            const ffft = try FftPlan(Complex(T)).init(W.rowSlice(0), shape, -1, plan_mode);
            errdefer ffft.deinit();

            const ifft = try FftPlan(Complex(T)).init(W.rowSlice(0), shape, 1, plan_mode);
            errdefer ifft.deinit();

            return .{ .W = W, .ffft = ffft, .ifft = ifft };
        }

        /// Deallocates wavefunction array and FFT plans.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.W.deinit(gpa);

            self.ffft.deinit();
            self.ifft.deinit();
        }

        /// Clones the wavefunction and its FFT plans into a new structure.
        pub fn clone(self: @This(), gpa: Allocator) !@This() {
            const ffft = try self.ffft.clone();
            errdefer ffft.deinit();

            const ifft = try self.ifft.clone();
            errdefer ifft.deinit();

            const W = try self.W.clone(gpa);
            errdefer W.deinit(gpa);

            return .{ .W = W, .ffft = ffft, .ifft = ifft };
        }

        /// Computes kinetic energy expectation value in momentum space.
        pub fn ekin(self: @This(), ham: Hamiltonian(T), grid: Grid(T)) T {
            var value: T = 0;

            for (0..self.W.nrow()) |i| for (0..self.W.rowSlice(i).len) |j| {
                value += self.W.rowSlice(i)[j].squaredMagnitude() * ham.K.at(j);
            };

            return value * grid.dk;
        }

        /// Computes potential energy expectation value in coordinate space.
        pub fn epot(self: @This(), ham: Hamiltonian(T), grid: Grid(T)) T {
            var value: T = 0;

            for (0..self.W.nrow()) |i| for (0..self.W.nrow()) |k| for (0..self.W.ncol()) |j| {
                const psi_i = self.W.at(i, j);
                const psi_k = self.W.at(k, j);

                value += (psi_i.re * psi_k.re + psi_i.im * psi_k.im) * ham.V.at(j, i * self.W.nrow() + k);
            };

            return value * grid.dr;
        }

        /// Transforms the wavefunction between position and momentum space.
        pub fn fft(self: *@This(), comptime sign: i32) void {
            for (0..self.W.nrow()) |i| {
                const slice = self.W.rowSlice(i);

                if (comptime sign == -1) self.ffft.execute(slice) else self.ifft.execute(slice);
            }
        }

        /// Computes momentum expectation value of the wavepacket.
        pub fn mom(self: @This(), grid: Grid(T), gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(grid.r.ncol(), gpa);

            for (0..self.W.nrow()) |i| for (0..self.W.rowSlice(i).len) |j| for (0..grid.r.ncol()) |k| {
                value.ptr(k).* += self.W.rowSlice(i)[j].squaredMagnitude() * grid.k.at(j, k);
            };

            value.muls(grid.dk);

            return value;
        }

        /// Computes spatial norm of the wavefunction.
        pub fn norm(self: @This(), grid: Grid(T)) T {
            var value: T = 0;

            for (0..self.W.nrow()) |i| for (0..self.W.rowSlice(i).len) |j| {
                value += self.W.rowSlice(i)[j].squaredMagnitude();
            };

            return value * grid.dr;
        }

        /// Normalizes the wavefunction to unit norm.
        pub fn normalize(self: *@This(), grid: Grid(T)) void {
            self.W.divs(Complex(T).init(std.math.sqrt(self.norm(grid)), 0));
        }

        /// Computes quantum mechanical overlap integral between two wavefunctions.
        pub fn overlap(self: @This(), other: @This(), grid: Grid(T)) Complex(T) {
            var value = Complex(T).init(0, 0);

            for (0..self.W.nrow()) |i| for (0..self.W.rowSlice(i).len) |j| {
                value = value.add(self.W.rowSlice(i)[j].conjugate().mul(other.W.rowSlice(i)[j]));
            };

            return value.mul(Complex(T).init(grid.dr, 0));
        }

        /// Computes diabatic populations of electronic states.
        pub fn pop(self: @This(), grid: Grid(T), gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(self.W.nrow(), gpa);

            for (0..self.W.nrow()) |i| {
                var sum: T = 0;

                for (0..self.W.rowSlice(i).len) |j| {
                    sum += self.W.rowSlice(i)[j].squaredMagnitude();
                }

                value.ptr(i).* = sum;
            }

            value.muls(grid.dr);

            return value;
        }

        /// Computes adiabatic populations of electronic states.
        pub fn popAdia(self: @This(), ham: Hamiltonian(T), grid: Grid(T), gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(self.W.nrow(), gpa);

            for (0..self.W.nrow()) |i| {
                for (0..self.W.ncol()) |j| {
                    var adia_re: T = 0;
                    var adia_im: T = 0;

                    for (0..self.W.nrow()) |k| {
                        const U_jk = ham.U.at(j, k * self.W.nrow() + i);

                        adia_re += self.W.at(k, j).re * U_jk;
                        adia_im += self.W.at(k, j).im * U_jk;
                    }

                    value.ptr(i).* += adia_re * adia_re + adia_im * adia_im;
                }
            }

            value.muls(grid.dr);

            return value;
        }

        /// Computes position expectation value of the wavepacket.
        pub fn pos(self: @This(), grid: Grid(T), gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(grid.r.ncol(), gpa);

            for (0..self.W.nrow()) |i| for (0..self.W.rowSlice(i).len) |j| for (0..grid.r.ncol()) |k| {
                value.ptr(k).* += self.W.rowSlice(i)[j].squaredMagnitude() * grid.r.at(j, k);
            };

            value.muls(grid.dr);

            return value;
        }

        /// Sets the wavefunction to a Gaussian wavepacket with specified phase.
        pub fn setGaussian(self: *@This(), ic: InitialConditions, grid: Grid(T)) void {
            self.W.fill(Complex(T).init(0, 0));

            for (0..grid.r.nrow()) |i| {
                var exponent = Complex(T).init(0, 0);

                for (0..grid.r.ncol()) |j| {
                    const dx = grid.r.at(i, j) - ic.position[j];

                    exponent = exponent.add(Complex(T).init(-0.5 * ic.gamma[j] * dx * dx, ic.momentum[j] * dx));
                }

                self.W.ptr(ic.state, i).* = std.math.complex.exp(exponent);
            }

            self.normalize(grid);
        }

        /// Transforms the wavepacket from diabatic to adiabatic representation.
        pub fn toAdia(self: *@This(), ham: Hamiltonian(T), gpa: Allocator) !void {
            var temp = try gpa.alloc(Complex(T), self.W.nrow());
            defer gpa.free(temp);

            for (0..self.W.ncol()) |j| {
                for (0..self.W.nrow()) |i| {
                    temp[i] = self.W.at(i, j);
                }

                for (0..self.W.nrow()) |i| {
                    var sum = Complex(T).init(0, 0);

                    for (0..self.W.nrow()) |k| {
                        const u_kj = Complex(T).init(ham.U.at(j, k * self.W.nrow() + i), 0);

                        sum = sum.add(temp[k].mul(u_kj));
                    }

                    self.W.ptr(i, j).* = sum;
                }
            }
        }

        /// Transforms the wavepacket from adiabatic to diabatic representation.
        pub fn toDia(self: *@This(), ham: Hamiltonian(T), gpa: Allocator) !void {
            var temp = try gpa.alloc(Complex(T), self.W.nrow());
            defer gpa.free(temp);

            for (0..self.W.ncol()) |j| {
                for (0..self.W.nrow()) |i| {
                    temp[i] = self.W.at(i, j);
                }

                for (0..self.W.nrow()) |i| {
                    var sum = Complex(T).init(0, 0);

                    for (0..self.W.nrow()) |k| {
                        const u_jk = Complex(T).init(ham.U.at(j, i * self.W.nrow() + k), 0);

                        sum = sum.add(temp[k].mul(u_jk));
                    }

                    self.W.ptr(i, j).* = sum;
                }
            }
        }
    };
}
