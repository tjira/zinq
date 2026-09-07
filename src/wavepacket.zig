//! Wavefunction and Hamiltonian representations for grid-based quantum wavepacket dynamics.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;

const FftPlan = @import("fourier_transform.zig").FftPlan;
const Matrix = @import("tensor.zig").Matrix;
const Potential = @import("potential.zig").Potential;
const Vector = @import("tensor.zig").Vector;

const eighBatch = @import("linear_algebra.zig").eighBatch;
const eighSlice = @import("linear_algebra.zig").eighSlice;

/// Generates a multidimensional discrete coordinate and momentum space grid.
pub fn Grid(comptime T: type) type {
    return struct {
        r: ?Matrix(T),
        k: ?Matrix(T),
        npoint: usize,

        dr: T,
        dk: T,

        lim: []const [2]T,
        cylindrical: bool,

        /// Allocates and initializes grid coordinates and momentum vectors.
        pub fn init(bounds: []const [2]T, npoint: u32, cylindrical: bool, optimize_memory: bool, gpa: Allocator) !@This() {
            const total_points = std.math.pow(usize, npoint, bounds.len);

            const bounds_copy = try gpa.alloc([2]T, bounds.len);
            errdefer gpa.free(bounds_copy);

            @memcpy(bounds_copy, bounds);

            var dr: T = 1;

            for (0..bounds.len) |i| {
                const min = bounds[i][0];
                const max = bounds[i][1];

                dr *= (max - min) / @as(T, @floatFromInt(npoint));
            }

            const dk = dr / @as(T, @floatFromInt(total_points));

            if (optimize_memory) {
                var grid: @This() = undefined;

                grid.r = null;
                grid.k = null;

                grid.lim, grid.npoint = .{ bounds_copy, npoint };

                grid.dr = dr;
                grid.dk = dk;

                grid.cylindrical = cylindrical;

                return grid;
            }

            var r = try Matrix(T).init(total_points, bounds.len, gpa);
            errdefer r.deinit(gpa);

            var k = try Matrix(T).init(total_points, bounds.len, gpa);
            errdefer k.deinit(gpa);

            for (0..total_points) |i| {
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

            var grid: @This() = undefined;

            grid.r = r;
            grid.k = k;

            grid.lim, grid.npoint = .{ bounds_copy, npoint };

            grid.dr = dr;
            grid.dk = dk;

            grid.cylindrical = cylindrical;

            return grid;
        }

        /// Deallocates coordinate and momentum space grid matrices and boundaries.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            if (self.r) |*r| r.deinit(gpa);
            if (self.k) |*k| k.deinit(gpa);

            gpa.free(self.lim);
        }

        /// Calculates the momentum space coordinate along dimension j at grid index i.
        pub fn getK(self: @This(), i: usize, j: usize) T {
            if (self.k) |k| return k.at(i, j);

            const exponent, var div: usize = .{ self.lim.len - j - 1, 1 };

            for (0..exponent) |_| {
                div *= self.npoint;
            }

            const n_idx = (i / div) % self.npoint;

            const min = self.lim[j][0];
            const max = self.lim[j][1];

            const dki, const n = .{ 2 * std.math.pi / (max - min), @as(T, @floatFromInt(n_idx)) };

            return (if (n_idx < self.npoint / 2) n else n - @as(T, @floatFromInt(self.npoint))) * dki;
        }

        /// Calculates the grid coordinate value along dimension j at grid index i.
        pub fn getR(self: @This(), i: usize, j: usize) T {
            if (self.r) |r| return r.at(i, j);

            const exponent, var div: usize = .{ self.lim.len - j - 1, 1 };

            for (0..exponent) |_| {
                div *= self.npoint;
            }

            const n_idx = (i / div) % self.npoint;

            const min = self.lim[j][0];
            const max = self.lim[j][1];

            return min + @as(T, @floatFromInt(n_idx)) * (max - min) / @as(T, @floatFromInt(self.npoint));
        }

        /// Returns the number of dimensions/axes in the multi-dimensional grid.
        pub fn ncol(self: @This()) usize {
            if (self.r) |r| return r.ncol();

            return self.lim.len;
        }

        /// Returns the total number of grid points in the discrete coordinates space.
        pub fn nrow(self: @This()) usize {
            if (self.r) |r| return r.nrow();

            return std.math.pow(usize, self.npoint, self.lim.len);
        }
    };
}

/// Representation of kinetic energy and potential energy operators on the grid.
pub fn Hamiltonian(comptime T: type) type {
    return struct {
        V: ?Matrix(T),
        W: ?Matrix(T),
        U: ?Matrix(T),
        K: ?Vector(T),

        mass: []const T,
        cylindric: bool,
        j_quantn: usize,

        w_buf: ?[]T,
        u_buf: ?[]T,
        v_buf: ?[]T,
        r_buf: ?[]T,

        /// Allocates and computes kinetic and potential operator matrix elements.
        pub fn init(grid: Grid(T), pot: Potential(T), m: []const T, j_quantn: u32, optimize_memory: bool, gpa: Allocator) !@This() {
            const mass = try gpa.alloc(T, m.len);
            errdefer gpa.free(mass);

            @memcpy(mass, m);

            if (optimize_memory) {
                const nstate = pot.nstate();

                const w_buf = try gpa.alloc(T, nstate);
                errdefer gpa.free(w_buf);

                const u_buf = try gpa.alloc(T, nstate * nstate);
                errdefer gpa.free(u_buf);

                const v_buf = try gpa.alloc(T, nstate * nstate);
                errdefer gpa.free(v_buf);

                const r_buf = try gpa.alloc(T, grid.ncol());
                errdefer gpa.free(r_buf);

                var ham: @This() = undefined;

                ham.V = null;
                ham.W = null;
                ham.U = null;
                ham.K = null;

                ham.mass, ham.cylindric, ham.j_quantn = .{ mass, grid.cylindrical, j_quantn };

                ham.w_buf = w_buf;
                ham.u_buf = u_buf;
                ham.v_buf = v_buf;
                ham.r_buf = r_buf;

                return ham;
            }

            var V = try Matrix(T).init(grid.nrow(), pot.nstate() * pot.nstate(), gpa);
            errdefer V.deinit(gpa);

            var U = try Matrix(T).init(grid.nrow(), pot.nstate() * pot.nstate(), gpa);
            errdefer U.deinit(gpa);

            var W = try Matrix(T).init(grid.nrow(), pot.nstate(), gpa);
            errdefer W.deinit(gpa);

            var K = try Vector(T).initZero(grid.nrow(), gpa);
            errdefer K.deinit(gpa);

            for (0..grid.nrow()) |i| {
                var sum: T = 0;

                for (0..grid.ncol()) |j| {
                    const kij = grid.getK(i, j);

                    sum += 0.5 * kij * kij / mass[j];
                }

                K.ptr(i).* = sum;
            }

            var ham: @This() = undefined;

            ham.V = V;
            ham.W = W;
            ham.U = U;
            ham.K = K;

            ham.mass, ham.cylindric, ham.j_quantn = .{ mass, grid.cylindrical, j_quantn };

            ham.w_buf = null;
            ham.u_buf = null;
            ham.v_buf = null;
            ham.r_buf = null;

            try ham.update(grid, pot, 0, gpa);

            return ham;
        }

        /// Deallocates Hamiltonian operator matrices.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            gpa.free(self.mass);

            if (self.w_buf) |buf| gpa.free(buf);
            if (self.u_buf) |buf| gpa.free(buf);
            if (self.v_buf) |buf| gpa.free(buf);
            if (self.r_buf) |buf| gpa.free(buf);

            if (self.V) |*V| V.deinit(gpa);
            if (self.W) |*W| W.deinit(gpa);
            if (self.U) |*U| U.deinit(gpa);
            if (self.K) |*K| K.deinit(gpa);
        }

        /// Returns the kinetic energy expectation value at grid index i.
        pub fn getK(self: @This(), grid: Grid(T), i: usize) T {
            if (self.K) |K| return K.at(i);

            var sum: T = 0;

            for (0..grid.ncol()) |j| {
                const kij = grid.getK(i, j);

                sum += 0.5 * kij * kij / self.mass[j];
            }

            return sum;
        }

        /// Computes or retrieves eigenvalues and eigenvectors at grid index i.
        pub fn getTriple(self: @This(), grid: Grid(T), pot: Potential(T), t: T, i: usize) !struct { []const T, []const T, []const T } {
            if (self.U) |U| {
                return .{ self.W.?.rowSlice(i), U.rowSlice(i), self.V.?.rowSlice(i) };
            }

            const w = self.w_buf.?;
            const u = self.u_buf.?;
            const v = self.v_buf.?;

            _ = self.getV(grid, pot, t, i);

            try eighSlice(T, w, u, v);

            return .{ w, u, v };
        }

        /// Computes or retrieves potential energy matrix elements at grid coordinate index i.
        pub fn getV(self: @This(), grid: Grid(T), pot: Potential(T), t: T, i: usize) []const T {
            if (self.V) |V| return V.rowSlice(i);

            const buffer, const r_coords = .{ self.v_buf.?, self.r_buf.? };

            for (0..grid.ncol()) |j| {
                r_coords[j] = grid.getR(i, j);
            }

            pot.eval(T, buffer, r_coords, t);

            if (self.cylindric) {
                const radial_idx = grid.ncol() - 1;

                const r, const m = .{ grid.getR(i, radial_idx), self.mass[radial_idx] };

                for (0..pot.nstate()) |s| {
                    buffer[s * pot.nstate() + s] -= if (r != 0) 1 / (8 * m * r * r) else 0;
                }
            }

            if (self.j_quantn > 0) {
                const j = @as(T, @floatFromInt(self.j_quantn));

                const r, const m = .{ grid.getR(i, 0), self.mass[0] };

                for (0..pot.nstate()) |s| {
                    buffer[s * pot.nstate() + s] += if (r != 0) j * (j + 1) / (2 * m * r * r) else 0;
                }
            }

            return buffer;
        }

        /// Updates potential energy values and diagonalizes to get adiabatic states.
        pub fn update(self: *@This(), grid: Grid(T), pot: Potential(T), t: T, gpa: Allocator) !void {
            if (self.V == null) return;

            var U_prev = if (pot.isTd() and t > 0) try self.U.?.clone(gpa) else null;
            defer if (U_prev) |*u| u.deinit(gpa);

            if (grid.r) |r| {
                pot.evalBatch(T, &self.V.?, r, t);
            }

            if (grid.r == null) {
                const r_coords = try gpa.alloc(T, grid.ncol());
                defer gpa.free(r_coords);

                for (0..grid.nrow()) |i| {
                    for (0..grid.ncol()) |j| {
                        r_coords[j] = grid.getR(i, j);
                    }

                    pot.eval(T, self.V.?.rowSlice(i), r_coords, t);
                }
            }

            if (self.cylindric) {
                const radial_idx = grid.ncol() - 1;

                for (0..grid.nrow()) |i| {
                    const r, const m = .{ grid.getR(i, radial_idx), self.mass[radial_idx] };

                    for (0..pot.nstate()) |s| {
                        self.V.?.ptr(i, s * pot.nstate() + s).* -= if (r != 0) 1 / (8 * m * r * r) else 0;
                    }
                }
            }

            if (self.j_quantn > 0) {
                const j, const m = .{ @as(T, @floatFromInt(self.j_quantn)), self.mass[0] };

                for (0..grid.nrow()) |i| {
                    const r = grid.getR(i, 0);

                    for (0..pot.nstate()) |s| {
                        self.V.?.ptr(i, s * pot.nstate() + s).* += if (r != 0) j * (j + 1) / (2 * m * r * r) else 0;
                    }
                }
            }

            try eighBatch(T, &self.W.?, &self.U.?, self.V.?);

            if (grid.ncol() == 1) for (1..grid.nrow()) |i| for (0..pot.nstate()) |j| {
                var overlap: T = 0;

                for (0..pot.nstate()) |k| {
                    overlap += self.U.?.at(i, k * pot.nstate() + j) * self.U.?.at(i - 1, k * pot.nstate() + j);
                }

                if (overlap < 0) for (0..pot.nstate()) |k| {
                    self.U.?.ptr(i, k * pot.nstate() + j).* = -self.U.?.at(i, k * pot.nstate() + j);
                };
            };

            if (U_prev) |prev| for (0..pot.nstate()) |j| {
                var total_overlap: T = 0;

                for (0..grid.nrow()) |i| for (0..pot.nstate()) |k| {
                    total_overlap += self.U.?.at(i, k * pot.nstate() + j) * prev.at(i, k * pot.nstate() + j);
                };

                if (total_overlap < 0) for (0..grid.nrow()) |i| for (0..pot.nstate()) |k| {
                    self.U.?.ptr(i, k * pot.nstate() + j).* = -self.U.?.at(i, k * pot.nstate() + j);
                };
            };
        }
    };
}

/// Initial parameters of the wavepacket including position, momentum, and Gaussian width.
pub const InitialConditions = struct {
    position: []const f64,
    momentum: []const f64,
    gamma: []const f64,

    state: u32 = 0,
    adiabatic: bool = false,
};

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

            var ffft = try FftPlan(Complex(T)).init(W.rowSlice(0), shape, -1, plan_mode);
            errdefer ffft.deinit();

            var ifft = try FftPlan(Complex(T)).init(W.rowSlice(0), shape, 1, plan_mode);
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
        pub fn ekin(self: @This(), ham: Hamiltonian(T), grid: Grid(T), langer: T) T {
            var value: T = 0;

            for (0..self.W.nrow()) |i| for (0..self.W.rowSlice(i).len) |j| {
                value += self.W.rowSlice(i)[j].squaredMagnitude() * ham.getK(grid, j);
            };

            return value * grid.dk - langer;
        }

        /// Computes potential energy expectation value in coordinate space.
        pub fn epot(self: @This(), ham: Hamiltonian(T), grid: Grid(T), pot: Potential(T), t: T, langer: T) T {
            var value: T = 0;

            for (0..self.W.ncol()) |j| {
                const V_j = ham.getV(grid, pot, t, j);

                for (0..self.W.nrow()) |i| {
                    for (0..self.W.nrow()) |k| {
                        const psi_i = self.W.at(i, j);
                        const psi_k = self.W.at(k, j);

                        value += (psi_i.re * psi_k.re + psi_i.im * psi_k.im) * V_j[i * self.W.nrow() + k];
                    }
                }
            }

            return value * grid.dr + langer;
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
            var value, const radial_idx = .{ try Vector(T).initZero(grid.ncol(), gpa), grid.ncol() - 1 };

            for (0..self.W.nrow()) |i| for (0..self.W.rowSlice(i).len) |j| for (0..grid.ncol()) |k| {
                const val = if (grid.cylindrical and k == radial_idx) @abs(grid.getK(j, k)) else grid.getK(j, k);

                value.ptr(k).* += self.W.rowSlice(i)[j].squaredMagnitude() * val;
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
        pub fn popAdia(self: @This(), ham: Hamiltonian(T), grid: Grid(T), pot: Potential(T), t: T, gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(self.W.nrow(), gpa);
            errdefer value.deinit(gpa);

            for (0..self.W.ncol()) |j| {
                const triple = try ham.getTriple(grid, pot, t, j);

                for (0..self.W.nrow()) |i| {
                    var adia_re: T = 0;
                    var adia_im: T = 0;

                    for (0..self.W.nrow()) |k| {
                        const U_jk = triple[1][k * self.W.nrow() + i];

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
            var value, const radial_idx = .{ try Vector(T).initZero(grid.ncol(), gpa), grid.ncol() - 1 };

            for (0..self.W.nrow()) |i| for (0..self.W.rowSlice(i).len) |j| for (0..grid.ncol()) |k| {
                const val = if (grid.cylindrical and k == radial_idx) @abs(grid.getR(j, k)) else grid.getR(j, k);

                value.ptr(k).* += self.W.rowSlice(i)[j].squaredMagnitude() * val;
            };

            value.muls(grid.dr);

            return value;
        }

        /// Sets the wavefunction to a Gaussian wavepacket with specified phase.
        pub fn setGaussian(self: *@This(), ic: InitialConditions, grid: Grid(T)) void {
            self.W.fill(Complex(T).init(0, 0));

            for (0..grid.nrow()) |i| {
                var exponent = Complex(T).init(0, 0);

                for (0..grid.ncol()) |j| {
                    const dx = grid.getR(i, j) - ic.position[j];

                    exponent = exponent.add(Complex(T).init(-0.5 * ic.gamma[j] * dx * dx, ic.momentum[j] * dx));
                }

                var val = std.math.complex.exp(exponent);

                if (grid.cylindrical) {
                    const r = grid.getR(i, grid.ncol() - 1);

                    val = val.mul(Complex(T).init(std.math.sign(r) * std.math.sqrt(@abs(r)), 0));
                }

                self.W.ptr(ic.state, i).* = val;
            }

            self.normalize(grid);
        }

        /// Transforms the wavepacket from diabatic to adiabatic representation.
        pub fn toAdia(self: *@This(), ham: Hamiltonian(T), grid: Grid(T), pot: Potential(T), t: T, gpa: Allocator) !void {
            var temp = try gpa.alloc(Complex(T), self.W.nrow());
            defer gpa.free(temp);

            for (0..self.W.ncol()) |j| {
                for (0..self.W.nrow()) |i| {
                    temp[i] = self.W.at(i, j);
                }

                _, const u, _ = try ham.getTriple(grid, pot, t, j);

                for (0..self.W.nrow()) |i| {
                    var sum = Complex(T).init(0, 0);

                    for (0..self.W.nrow()) |k| {
                        const u_kj = u[k * self.W.nrow() + i];

                        sum.re += temp[k].re * u_kj;
                        sum.im += temp[k].im * u_kj;
                    }

                    self.W.ptr(i, j).* = sum;
                }
            }
        }

        /// Transforms the wavepacket from adiabatic to diabatic representation.
        pub fn toDia(self: *@This(), ham: Hamiltonian(T), grid: Grid(T), pot: Potential(T), t: T, gpa: Allocator) !void {
            var temp = try gpa.alloc(Complex(T), self.W.nrow());
            defer gpa.free(temp);

            for (0..self.W.ncol()) |j| {
                for (0..self.W.nrow()) |i| {
                    temp[i] = self.W.at(i, j);
                }

                _, const u, _ = try ham.getTriple(grid, pot, t, j);

                for (0..self.W.nrow()) |i| {
                    var sum = Complex(T).init(0, 0);

                    for (0..self.W.nrow()) |k| {
                        const u_jk = u[i * self.W.nrow() + k];

                        sum.re += temp[k].re * u_jk;
                        sum.im += temp[k].im * u_jk;
                    }

                    self.W.ptr(i, j).* = sum;
                }
            }
        }
    };
}
