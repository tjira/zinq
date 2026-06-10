const std = @import("std");

const fftw = @cImport(@cInclude("fftw3.h"));

// zig fmt: off
const Allocator = std.mem .Allocator;
const Complex   = std.math.  Complex;
// zig fmt: on

// zig fmt: off
const FftPlan          = @import("fftw.zig"     ).  FftPlan;
const Matrix           = @import("tensor.zig"   ).   Matrix;
const Potential        = @import("potential.zig").Potential;
const PotentialOptions = @import("potential.zig").  Options;
const Vector           = @import("tensor.zig"   ).   Vector;
// zig fmt: on

// zig fmt: off
const eighBatch         = @import("openblas.zig"  )        .eighBatch;
const printf            = @import("read_write.zig")           .printf;
const writeMatrixHjoin  = @import("read_write.zig") .writeMatrixHjoin;
const writeMatrixLspace = @import("read_write.zig").writeMatrixLspace;
// zig fmt: on

// OPTIONS =============================================================================================================

// zig fmt: off
const InitialConditions = struct {
    position: []const f64,
    momentum: []const f64,
    gamma:    []const f64,

    state:     u32 =      0,
    adiabatic: bool = false,
};
// zig fmt: on

const Imaginary = struct {
    nstate: u32 = 1,
};

// zig fmt: off
const Write = struct {
    kinetic_energy:   ?[]const u8 = null,
    momentum:         ?[]const u8 = null,
    norm:             ?[]const u8 = null,
    population:       ?[]const u8 = null,
    position:         ?[]const u8 = null,
    potential_energy: ?[]const u8 = null,
    total_energy:     ?[]const u8 = null,
    wavefunction:     ?[]const u8 = null,
};
// zig fmt: on

// zig fmt: off
pub const Options = struct {
    const FftOpt = struct {
        plan: enum { estimate, measure, patient, exhaustive } = .measure,
    };
    const GridOpt = struct {
        bounds: []const [2]f64,
        npoint:            u32,
    };

    grid:                         GridOpt,
    initial_conditions: InitialConditions,
    potential:           PotentialOptions,
    time_step:                        f64,
    iterations:                       u32,
    mass:                             f64,

    fft:          FftOpt    =   .{},
    imaginary:   ?Imaginary =  null,
    write:        Write     =   .{},
    adiabatic:    bool      = false,
    log_interval: u32       =     1,
};
// zig fmt: on

// GRID ================================================================================================================

fn Grid(comptime T: type) type {
    return struct {
        r: Matrix(T),
        k: Matrix(T),

        dr: T,
        dk: T,

        pub fn init(bounds: []const [2]T, npoint: u32, gpa: Allocator) !@This() {
            const ncol = std.math.pow(usize, npoint, bounds.len);

            var r = try Matrix(T).init(ncol, bounds.len, gpa);
            var k = try Matrix(T).init(ncol, bounds.len, gpa);

            var dr: T = 1;

            for (bounds) |e| {
                const min = e[0];
                const max = e[1];

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

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.r.deinit(gpa);
            self.k.deinit(gpa);
        }
    };
}

// WAVEFUNCTION ========================================================================================================

fn Wavefunction(comptime T: type) type {
    return struct {
        W: Matrix(Complex(T)),

        ffft: []FftPlan(Complex(T)),
        ifft: []FftPlan(Complex(T)),

        pub fn init(ndim: usize, nstate: usize, npoint: usize, plan_mode: u32, gpa: Allocator) !@This() {
            const W = try Matrix(Complex(T)).init(nstate, std.math.pow(usize, npoint, ndim), gpa);

            const ffft = try gpa.alloc(FftPlan(Complex(T)), nstate);
            const ifft = try gpa.alloc(FftPlan(Complex(T)), nstate);

            const shape = try gpa.alloc(i32, ndim);

            for (shape) |*e| {
                e.* = @as(i32, @intCast(npoint));
            }

            for (0..nstate) |i| {
                // zig fmt: off
                ffft[i] = try FftPlan(Complex(T)).init(W.rowSlice(i), shape, -1, plan_mode);
                ifft[i] = try FftPlan(Complex(T)).init(W.rowSlice(i), shape,  1, plan_mode);
                // zig fmt: on
            }

            gpa.free(shape);

            return .{ .W = W, .ffft = ffft, .ifft = ifft };
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.W.deinit(gpa);

            for (0..self.ffft.len) |i| {
                self.ffft[i].deinit();
                self.ifft[i].deinit();
            }

            gpa.free(self.ffft);
            gpa.free(self.ifft);
        }

        pub fn clone(self: @This(), gpa: Allocator) !@This() {
            var ffft = try gpa.alloc(FftPlan(Complex(T)), self.ffft.len);
            var ifft = try gpa.alloc(FftPlan(Complex(T)), self.ifft.len);

            for (0..self.ffft.len) |i| {
                ffft[i] = try self.ffft[i].clone();
                ifft[i] = try self.ifft[i].clone();
            }

            return .{ .W = try self.W.clone(gpa), .ffft = ffft, .ifft = ifft };
        }

        pub fn ekin(self: @This(), ham: Hamiltonian(T), grid: Grid(T)) T {
            var value: T = 0;

            for (0..self.W.nrow()) |i| for (self.W.rowSlice(i), 0..) |e, j| {
                value += e.squaredMagnitude() * ham.K.at(j);
            };

            return value * grid.dk;
        }

        pub fn epot(self: @This(), ham: Hamiltonian(T), grid: Grid(T)) T {
            var value: T = 0;

            for (0..self.W.nrow()) |i| for (0..self.W.nrow()) |k| for (0..self.W.ncol()) |j| {
                const psi_i = self.W.at(i, j);
                const psi_k = self.W.at(k, j);

                value += (psi_i.re * psi_k.re + psi_i.im * psi_k.im) * ham.V.at(j, i * self.W.nrow() + k);
            };

            return value * grid.dr;
        }

        pub fn fft(self: *@This(), comptime sign: i32) !void {
            for (0..self.W.nrow()) |i| {
                const slice = self.W.rowSlice(i);

                if (comptime sign == -1) self.ffft[i].execute(slice) else self.ifft[i].execute(slice);
            }
        }

        pub fn mom(self: @This(), grid: Grid(T), gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(grid.r.ncol(), gpa);

            for (0..self.W.nrow()) |i| for (self.W.rowSlice(i), 0..) |e, j| {
                const mag = e.squaredMagnitude();

                for (0..grid.r.ncol()) |k| {
                    value.ptr(k).* += mag * grid.k.at(j, k);
                }
            };

            value.muls(grid.dk);

            return value;
        }

        pub fn norm(self: @This(), grid: Grid(T)) T {
            var value: T = 0;

            for (0..self.W.nrow()) |i| for (self.W.rowSlice(i)) |e| {
                value += e.squaredMagnitude();
            };

            return value * grid.dr;
        }

        pub fn normalize(self: *@This(), grid: Grid(T)) void {
            self.W.divs(Complex(T).init(std.math.sqrt(self.norm(grid)), 0));
        }

        pub fn overlap(self: @This(), other: @This(), grid: Grid(T)) Complex(T) {
            var value = Complex(T).init(0, 0);

            for (0..self.W.nrow()) |i| for (self.W.rowSlice(i), other.W.rowSlice(i)) |s, o| {
                value = value.add(s.conjugate().mul(o));
            };

            return value.mul(Complex(T).init(grid.dr, 0));
        }

        pub fn pop(self: @This(), grid: Grid(T), gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(self.W.nrow(), gpa);

            for (0..self.W.nrow()) |i| {
                var sum: T = 0;

                for (self.W.rowSlice(i)) |e| {
                    sum += e.squaredMagnitude();
                }

                value.ptr(i).* = sum;
            }

            value.muls(grid.dr);

            return value;
        }

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

        pub fn pos(self: @This(), grid: Grid(T), gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(grid.r.ncol(), gpa);

            for (0..self.W.nrow()) |i| for (self.W.rowSlice(i), 0..) |e, j| {
                const mag = e.squaredMagnitude();

                for (0..grid.r.ncol()) |k| {
                    value.ptr(k).* += mag * grid.r.at(j, k);
                }
            };

            value.muls(grid.dr);

            return value;
        }

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

        pub fn toAdia(self: *@This(), ham: Hamiltonian(T), gpa: Allocator) !void {
            var temp = try gpa.alloc(Complex(T), self.W.nrow());

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

            gpa.free(temp);
        }

        pub fn toDia(self: *@This(), ham: Hamiltonian(T), gpa: Allocator) !void {
            var temp = try gpa.alloc(Complex(T), self.W.nrow());

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

            gpa.free(temp);
        }
    };
}

// HAMILTONIAN =========================================================================================================

fn Hamiltonian(comptime T: type) type {
    return struct {
        V: Matrix(T),
        W: Matrix(T),
        U: Matrix(T),
        K: Vector(T),

        pub fn init(grid: Grid(T), pot: Potential(T), m: T, gpa: Allocator) !@This() {
            const V = try Matrix(T).init(grid.r.nrow(), pot.nstate() * pot.nstate(), gpa);
            const U = try Matrix(T).init(grid.r.nrow(), pot.nstate() * pot.nstate(), gpa);

            const W = try Matrix(T).init(grid.r.nrow(), pot.nstate(), gpa);

            var K = try Vector(T).initZero(grid.r.nrow(), gpa);

            for (0..grid.r.nrow()) |i| {
                var sum: T = 0;

                for (0..grid.r.ncol()) |j| {
                    const kij = grid.k.at(i, j);

                    sum += 0.5 * kij * kij / m;
                }

                K.ptr(i).* = sum;
            }

            var ham = @This(){ .V = V, .W = W, .U = U, .K = K };

            try ham.update(grid, pot, 0);

            return ham;
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.V.deinit(gpa);
            self.W.deinit(gpa);
            self.U.deinit(gpa);
            self.K.deinit(gpa);
        }

        pub fn update(self: *@This(), grid: Grid(T), pot: Potential(T), t: T) !void {
            pot.evalBatch(T, &self.V, grid.r, t);

            try eighBatch(T, &self.W, &self.U, self.V);
        }
    };
}

// PROPAGATOR ==========================================================================================================

fn Propagator(comptime T: type) type {
    return struct {
        R: Matrix(Complex(T)),
        K: Vector(Complex(T)),

        dt: Complex(T),

        pub fn init(ham: Hamiltonian(T), dt: Complex(T), gpa: Allocator) !@This() {
            const R = try Matrix(Complex(T)).init(ham.V.nrow(), ham.V.ncol(), gpa);

            var K = try Vector(Complex(T)).initZero(ham.K.length(), gpa);

            for (ham.K.data, 0..) |e, i| {
                K.data[i] = std.math.complex.exp(Complex(T).init(0, -e).mul(dt));
            }

            var prop = @This(){ .R = R, .K = K, .dt = dt };

            prop.update(ham);

            return prop;
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.R.deinit(gpa);
            self.K.deinit(gpa);
        }

        pub fn step(self: @This(), wfn: *Wavefunction(T), gpa: Allocator) !void {
            // zig fmt: off
            try self.applyR(wfn, gpa);
            try self.applyK(wfn     );
            try self.applyR(wfn, gpa);
            // zig fmt: on
        }

        pub fn update(self: *@This(), ham: Hamiltonian(T)) void {
            const nstate = std.math.sqrt(ham.V.ncol());

            for (0..self.R.nrow()) |i| for (0..nstate) |k| for (0..nstate) |j| {
                var sum = Complex(T).init(0, 0);

                for (0..nstate) |m| {
                    const U_km = ham.U.at(i, k * nstate + m);
                    const U_jm = ham.U.at(i, j * nstate + m);

                    const phase = std.math.complex.exp(Complex(T).init(0, -0.5 * ham.W.at(i, m)).mul(self.dt));

                    sum = sum.add(Complex(T).init(phase.re * U_km * U_jm, phase.im * U_km * U_jm));
                }

                self.R.ptr(i, k * nstate + j).* = sum;
            };
        }

        fn applyK(self: @This(), wfn: *Wavefunction(T)) !void {
            try wfn.fft(-1);

            for (0..wfn.W.nrow()) |i| for (0..wfn.W.ncol()) |j| {
                wfn.W.ptr(i, j).* = self.K.at(j).mul(wfn.W.at(i, j));
            };

            try wfn.fft(1);
        }

        fn applyR(self: @This(), wfn: *Wavefunction(T), gpa: Allocator) !void {
            var temp = try gpa.alloc(Complex(T), wfn.W.nrow());

            for (0..wfn.W.ncol()) |j| {
                for (0..wfn.W.nrow()) |i| {
                    var sum = Complex(T).init(0, 0);

                    for (0..wfn.W.nrow()) |k| {
                        sum = sum.add(self.R.at(j, i * wfn.W.nrow() + k).mul(wfn.W.at(k, j)));
                    }

                    temp[i] = sum;
                }

                for (0..wfn.W.nrow()) |i| {
                    wfn.W.ptr(i, j).* = temp[i];
                }
            }

            gpa.free(temp);
        }
    };
}

// OBSERVABLES =========================================================================================================

fn Observables(comptime T: type) type {
    return struct {
        pos: ?Vector(T) = null,
        mom: ?Vector(T) = null,
        pop: ?Vector(T) = null,

        epot: ?T = null,
        ekin: ?T = null,
        norm: ?T = null,

        pub fn init(qsys: *QuantumSystem(T), write: Write, adia: bool, log: bool, gpa: Allocator) !@This() {
            var obs = @This(){};

            // zig fmt: off
            var calc = .{
                .pos  = log or write.position         != null,
                .norm = log or write.norm             != null,
                .epot = log or write.potential_energy != null,
                .pop  = log or write.population       != null,
                .mom  = log or write.momentum         != null,
                .ekin = log or write.kinetic_energy   != null,
            };
            // zig fmt: on

            calc.ekin = calc.ekin or write.total_energy != null;
            calc.epot = calc.epot or write.total_energy != null;

            // zig fmt: off
            if (calc.pos ) obs.pos  = try qsys.wfn.pos (          qsys.grid, gpa);
            if (calc.norm) obs.norm =     qsys.wfn.norm(          qsys.grid     );
            if (calc.epot) obs.epot =     qsys.wfn.epot(qsys.ham, qsys.grid     );
            // zig fmt: on

            if (calc.pop) {
                // zig fmt: off
                if (adia == true) obs.pop = try qsys.wfn.popAdia(qsys.ham,  qsys.grid, gpa);
                if (adia != true) obs.pop = try qsys.wfn.pop    (qsys.grid,            gpa);
                // zig fmt: on
            }

            const needs_fft = calc.mom or calc.ekin;

            if (needs_fft) {
                try qsys.wfn.fft(-1);

                // zig fmt: off
                if (calc.mom ) obs.mom  = try qsys.wfn.mom (          qsys.grid, gpa);
                if (calc.ekin) obs.ekin =     qsys.wfn.ekin(qsys.ham, qsys.grid     );
                // zig fmt: on

                try qsys.wfn.fft(1);
            }

            return obs;
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            if (self.pos) |*pos| pos.deinit(gpa);
            if (self.mom) |*mom| mom.deinit(gpa);
            if (self.pop) |*pop| pop.deinit(gpa);
        }
    };
}

// HISTORY =============================================================================================================

fn History(comptime T: type) type {
    return struct {
        // zig fmt: off
        pos:  ?Matrix(T),
        mom:  ?Matrix(T),
        pop:  ?Matrix(T),
        epot: ?Matrix(T),
        ekin: ?Matrix(T),
        etot: ?Matrix(T),
        norm: ?Matrix(T),
        wfn:  ?Matrix(T),
        // zig fmt: on

        index: usize = 0,

        pub fn init(ndim: usize, nstate: usize, npoint: usize, iters: usize, write: Write, gpa: Allocator) !@This() {
            const store_wfn = write.wavefunction != null;

            // zig fmt: off
            const store_epot = write.potential_energy != null or write.total_energy != null;
            const store_ekin = write.kinetic_energy   != null or write.total_energy != null;
            // zig fmt: on

            const wfn_nrow = std.math.pow(usize, npoint, ndim);

            return .{
                // zig fmt: off
                .pos  = if (write.position     != null) try Matrix(T).init(iters, ndim,   gpa) else null,
                .mom  = if (write.momentum     != null) try Matrix(T).init(iters, ndim,   gpa) else null,
                .pop  = if (write.population   != null) try Matrix(T).init(iters, nstate, gpa) else null,
                .norm = if (write.norm         != null) try Matrix(T).init(iters, 1,      gpa) else null,
                .etot = if (write.total_energy != null) try Matrix(T).init(iters, 1,      gpa) else null,
                // zig fmt: on

                .wfn = if (store_wfn) try Matrix(T).init(wfn_nrow, 2 * nstate * iters, gpa) else null,

                .ekin = if (store_ekin) try Matrix(T).init(iters, 1, gpa) else null,
                .epot = if (store_epot) try Matrix(T).init(iters, 1, gpa) else null,
            };
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            // zig fmt: off
            if (self.pos)  |*pos |  pos.deinit(gpa);
            if (self.mom)  |*mom |  mom.deinit(gpa);
            if (self.pop)  |*pop |  pop.deinit(gpa);
            if (self.epot) |*epot| epot.deinit(gpa);
            if (self.ekin) |*ekin| ekin.deinit(gpa);
            if (self.etot) |*etot| etot.deinit(gpa);
            if (self.norm) |*norm| norm.deinit(gpa);
            if (self.wfn)  |*wfn |  wfn.deinit(gpa);
            // zig fmt: on
        }

        pub fn append(self: *@This(), wfn: Wavefunction(T), obs: Observables(T)) void {
            const step_idx = self.index;

            if (self.wfn) |*storage| for (0..wfn.W.ncol()) |j| for (0..wfn.W.nrow()) |i| {
                storage.ptr(j, 2 * (step_idx * wfn.W.nrow() + i) + 0).* = wfn.W.at(i, j).re;
                storage.ptr(j, 2 * (step_idx * wfn.W.nrow() + i) + 1).* = wfn.W.at(i, j).im;
            };

            if (self.pos) |*pos| if (obs.pos) |v| {
                for (0..v.length()) |j| pos.ptr(step_idx, j).* = v.at(j);
            };

            if (self.mom) |*mom| if (obs.mom) |v| {
                for (0..v.length()) |j| mom.ptr(step_idx, j).* = v.at(j);
            };

            if (self.pop) |*pop| if (obs.pop) |v| {
                for (0..v.length()) |j| pop.ptr(step_idx, j).* = v.at(j);
            };

            // zig fmt: off
            if (self.epot) |*epot| {if (obs.epot) |v| epot.ptr(step_idx, 0).* = v;}
            if (self.ekin) |*ekin| {if (obs.ekin) |v| ekin.ptr(step_idx, 0).* = v;}
            if (self.norm) |*norm| {if (obs.norm) |v| norm.ptr(step_idx, 0).* = v;}
            // zig fmt: on

            if (self.etot) |*etot| {
                etot.ptr(step_idx, 0).* = obs.ekin.? + obs.epot.?;
            }

            self.index += 1;
        }

        pub fn exportAndDeinit(self: *@This(), io: std.Io, dt: f64, grid: Grid(T), write: Write, gpa: Allocator) !void {
            defer self.deinit(gpa);

            const end = dt * @as(T, @floatFromInt(self.index - 1));

            // zig fmt: off
            if (write.position        ) |path| try writeMatrixLspace(T, io, path, self.pos.?,  0, end);
            if (write.momentum        ) |path| try writeMatrixLspace(T, io, path, self.mom.?,  0, end);
            if (write.population      ) |path| try writeMatrixLspace(T, io, path, self.pop.?,  0, end);
            if (write.potential_energy) |path| try writeMatrixLspace(T, io, path, self.epot.?, 0, end);
            if (write.kinetic_energy  ) |path| try writeMatrixLspace(T, io, path, self.ekin.?, 0, end);
            if (write.norm            ) |path| try writeMatrixLspace(T, io, path, self.norm.?, 0, end);
            if (write.total_energy    ) |path| try writeMatrixLspace(T, io, path, self.etot.?, 0, end);
            // zig fmt: on

            if (write.wavefunction) |path| try writeMatrixHjoin(T, io, path, grid.r, self.wfn.?);
        }
    };
}

// LOGGERS =============================================================================================================

fn printHeader(io: std.Io, eigs: usize, ndim: usize, nstate: usize, neig: usize) !void {
    if (neig > 1) {
        try printf(io, "\nITP OF STATE {d}/{d}", .{ eigs + 1, neig });
    }

    const fmt = "\n{[0]s:8} {[1]s:12} {[2]s:12} {[3]s:12} {[4]s:[5]} {[6]s:[7]} {[8]s:[9]} {[10]s:10} {[11]s:4}\n";

    const tuple = .{
        "ITER",

        "EKIN (Eh)",
        "EPOT (Eh)",
        "ETOT (Eh)",

        // zig fmt: off
        "POS (a0)",    12 *   ndim,
        "MOM (hb/a0)", 12 *   ndim,
        "POP (-)",     11 * nstate,
        // zig fmt: on

        "NORM (-)",
        "TIME",
    };

    try printf(io, fmt, tuple);
}

fn printFinalEnergies(comptime T: type, io: std.Io, obs: std.ArrayList(Observables(T))) !void {
    try std.Io.File.stdout().writeStreamingAll(io, "\n");

    for (obs.items, 0..) |e, i| {
        const ekin = e.ekin orelse std.math.nan(T);
        const epot = e.epot orelse std.math.nan(T);

        const etot = ekin + epot;

        try printf(io, "FINAL ENERGY OF VIBRONIC STATE {d:02}: {d:13.8} Eh\n", .{ i, etot });
    }
}

fn printFinalPop(comptime T: type, io: std.Io, obs: Observables(T)) !void {
    if (obs.pop) |pop| {
        try std.Io.File.stdout().writeStreamingAll(io, "\n");

        for (0..pop.length()) |i| {
            try printf(io, "FINAL POPULATION OF ELECTRONIC STATE {d:02}: {d:.8}\n", .{ i, pop.at(i) });
        }
    }
}

fn printIteration(comptime T: type, io: std.Io, obs: Observables(T), i: usize, timer: *std.Io.Timestamp) !void {
    const ekin = obs.ekin orelse std.math.nan(T);
    const epot = obs.epot orelse std.math.nan(T);

    const etot = ekin + epot;

    try printf(io, "{d:8} {d:12.6} {d:12.6} {d:12.6} ", .{ i, ekin, epot, etot });

    if (obs.pos) |pos| {
        try printf(io, "[", .{});

        for (0..pos.length()) |j| {
            try printf(io, "{d:10.4}{s}", .{ pos.at(j), if (j == pos.length() - 1) "" else ", " });
        }

        try printf(io, "] ", .{});
    }

    if (obs.mom) |mom| {
        try printf(io, "[", .{});

        for (0..mom.length()) |j| {
            try printf(io, "{d:10.4}{s}", .{ mom.at(j), if (j == mom.length() - 1) "" else ", " });
        }

        try printf(io, "] ", .{});
    }

    if (obs.pop) |pop| {
        try printf(io, "[", .{});

        for (0..pop.length()) |j| {
            try printf(io, "{d:9.4}{s}", .{ pop.at(j), if (j == pop.length() - 1) "" else ", " });
        }

        try printf(io, "] ", .{});
    }

    if (obs.norm) |norm| {
        try printf(io, "{e:10.3} ", .{norm});
    }

    try printf(io, "{f}\n", .{timer.untilNow(io, .real)});

    timer.* = std.Io.Timestamp.now(io, .real);
}

// QUANTUM SYSTEM ======================================================================================================

fn QuantumSystem(comptime T: type) type {
    return struct {
        // zig fmt: off
        grid: Grid(T), ham: Hamiltonian(T), pot: Potential(T), wfn: Wavefunction(T),
        // zig fmt: on

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.grid.deinit(gpa);

            self.ham.deinit(gpa);
            self.pot.deinit(gpa);
            self.wfn.deinit(gpa);
        }
    };
}

// RUN =================================================================================================================

fn SimulationState(comptime T: type) type {
    return struct {
        // zig fmt: off
        qsys: QuantumSystem(T), prop: Propagator(T), orthw: std.ArrayList(Wavefunction(T)),
        // zig fmt: on

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.qsys.deinit(gpa);
            self.prop.deinit(gpa);

            for (self.orthw.items) |*e| e.deinit(gpa);

            self.orthw.deinit(gpa);
        }
    };
}

fn SolveContext(comptime T: type) type {
    return struct { opt: Options, sim: *SimulationState(T), eigs: usize, log: bool };
}

fn init(comptime T: type, opt: Options, gpa: Allocator) !SimulationState(T) {
    const pot = Potential(T).init(opt.potential);

    const dt = if (opt.imaginary) |_| Complex(T).init(0, -opt.time_step) else Complex(T).init(opt.time_step, 0);

    const plan_mode = switch (opt.fft.plan) {
        // zig fmt: off
        .estimate   => fftw  .FFTW_ESTIMATE,
        .measure    => fftw   .FFTW_MEASURE,
        .patient    => fftw   .FFTW_PATIENT,
        .exhaustive => fftw.FFTW_EXHAUSTIVE,
        // zig fmt: on
    };

    // zig fmt: off
    const grid = try Grid(T)        .init(opt.grid.bounds, opt.grid.npoint,                             gpa);
    const wfn  = try Wavefunction(T).init(pot.ndim(),      pot.nstate(),    opt.grid.npoint, plan_mode, gpa);
    const ham  = try Hamiltonian(T) .init(grid,            pot,             opt.mass,                   gpa);
    const prop = try Propagator(T)  .init(ham,             dt,                                          gpa);
    // zig fmt: on

    return .{ .qsys = .{ .grid = grid, .ham = ham, .wfn = wfn, .pot = pot }, .prop = prop, .orthw = .empty };
}

fn solve(comptime T: type, io: std.Io, ctx: SolveContext(T), gpa: Allocator, arena: Allocator) !Observables(T) {
    // zig fmt: off
    const ndim   = ctx.sim.qsys.grid.r.ncol();
    const nstate =  ctx.sim.qsys.wfn.W.nrow();
    const iters  =         ctx.opt.iterations;
    // zig fmt: on

    const neig = if (ctx.opt.imaginary) |imag| imag.nstate else 1;

    if (ctx.log) try printHeader(io, ctx.eigs, ctx.sim.qsys.grid.r.ncol(), ctx.sim.qsys.wfn.W.nrow(), neig);

    ctx.sim.qsys.wfn.setGaussian(ctx.opt.initial_conditions, ctx.sim.qsys.grid);

    var hist = try History(T).init(ndim, nstate, ctx.opt.grid.npoint, iters + 1, ctx.opt.write, gpa);

    if (ctx.opt.initial_conditions.adiabatic) {
        try ctx.sim.qsys.wfn.toDia(ctx.sim.qsys.ham, gpa);
    }

    var timer = std.Io.Timestamp.now(io, .real);

    for (0..ctx.opt.iterations + 1) |i| {
        const time = (@as(T, @floatFromInt(i)) - 0.5) * ctx.opt.time_step;

        if (i > 0 and ctx.sim.qsys.pot.isTd()) {
            try ctx.sim.qsys.ham.update(ctx.sim.qsys.grid, ctx.sim.qsys.pot, time);

            ctx.sim.prop.update(ctx.sim.qsys.ham);
        }

        if (i > 0) try ctx.sim.prop.step(&ctx.sim.qsys.wfn, gpa);

        if (ctx.opt.imaginary != null) for (ctx.sim.orthw.items) |e| {
            const overlap = e.overlap(ctx.sim.qsys.wfn, ctx.sim.qsys.grid);

            for (0..ctx.sim.qsys.wfn.W.data.len) |k| {
                ctx.sim.qsys.wfn.W.data[k] = ctx.sim.qsys.wfn.W.data[k].sub(overlap.mul(e.W.data[k]));
            }
        };

        if (ctx.opt.imaginary != null) ctx.sim.qsys.wfn.normalize(ctx.sim.qsys.grid);

        const is_log_step = ctx.log and ((i % ctx.opt.log_interval == 0) or (i == ctx.opt.iterations));

        if (ctx.sim.qsys.pot.isTd()) {
            const t = @as(T, @floatFromInt(i)) * ctx.opt.time_step;

            try ctx.sim.qsys.ham.update(ctx.sim.qsys.grid, ctx.sim.qsys.pot, t);
        }

        var obs = try Observables(T).init(&ctx.sim.qsys, ctx.opt.write, ctx.opt.adiabatic, is_log_step, gpa);

        defer obs.deinit(gpa);

        if (ctx.opt.write.wavefunction != null and ctx.opt.adiabatic) {
            try ctx.sim.qsys.wfn.toAdia(ctx.sim.qsys.ham, gpa);
        }

        hist.append(ctx.sim.qsys.wfn, obs);

        if (ctx.opt.write.wavefunction != null and ctx.opt.adiabatic) {
            try ctx.sim.qsys.wfn.toDia(ctx.sim.qsys.ham, gpa);
        }

        if (is_log_step) {
            try printIteration(T, io, obs, i, &timer);
        }
    }

    try hist.exportAndDeinit(io, ctx.opt.time_step, ctx.sim.qsys.grid, ctx.opt.write, gpa);

    return try Observables(T).init(&ctx.sim.qsys, ctx.opt.write, ctx.opt.adiabatic, true, arena);
}

pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator, arena: Allocator) !std.ArrayList(Observables(T)) {
    var output: std.ArrayList(Observables(T)) = .empty;

    if (log) try std.Io.File.stdout().writeStreamingAll(io, "\nQUANTUM DYNAMICS INIT: ");

    var timer = std.Io.Timestamp.now(io, .real);

    var sim = try init(T, opt, gpa);

    defer sim.deinit(gpa);

    if (log) try printf(io, "{f}\n", .{timer.untilNow(io, .real)});

    for (0..if (opt.imaginary) |imag| imag.nstate else 1) |i| {
        const obs = try solve(T, io, .{ .opt = opt, .sim = &sim, .eigs = i, .log = log }, gpa, arena);

        if (i < if (opt.imaginary) |imag| imag.nstate else 0) {
            try sim.orthw.append(gpa, try sim.qsys.wfn.clone(gpa));
        }

        try output.append(arena, obs);

        if (log) {
            try printFinalPop(T, io, obs);
        }
    }

    if (log) try printFinalEnergies(T, io, output);

    return output;
}
