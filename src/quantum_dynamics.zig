const std = @import("std");

const fftw = @cImport(@cInclude("fftw3.h"));

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;

const FftPlan = @import("fftw.zig").FftPlan;
const Matrix = @import("tensor.zig").Matrix;
const Potential = @import("potential.zig").Potential;
const PotentialOptions = @import("potential.zig").Options;
const Vector = @import("tensor.zig").Vector;

const eighBatch = @import("openblas.zig").eighBatch;
const printf = @import("read_write.zig").printf;
const writeMatrixHjoin = @import("read_write.zig").writeMatrixHjoin;
const writeMatrixLspace = @import("read_write.zig").writeMatrixLspace;

// OPTIONS =============================================================================================================

const InitialConditions = struct {
    position: []const f64,
    momentum: []const f64,
    gamma: []const f64,

    state: u32 = 0,
    adiabatic: bool = false,
};

const Imaginary = struct {
    nstate: u32 = 1,
};

const Write = struct {
    kinetic_energy: ?[]const u8 = null,
    momentum: ?[]const u8 = null,
    norm: ?[]const u8 = null,
    population: ?[]const u8 = null,
    position: ?[]const u8 = null,
    potential_energy: ?[]const u8 = null,
    total_energy: ?[]const u8 = null,
    wavefunction: ?[]const u8 = null,
};

pub const Options = struct {
    const FftOpt = struct {
        plan: enum { estimate, measure, patient, exhaustive } = .measure,
    };
    const GridOpt = struct {
        bounds: []const [2]f64,
        npoint: u32,
    };

    grid: GridOpt,
    initial_conditions: InitialConditions,
    potential: PotentialOptions,

    time_step: f64,
    iterations: u32,
    mass: f64,

    fft: FftOpt = .{},
    imaginary: ?Imaginary = null,
    write: Write = .{},

    adiabatic: bool = false,
    log_interval: u32 = 1,
};

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
            errdefer r.deinit(gpa);

            var k = try Matrix(T).init(ncol, bounds.len, gpa);
            errdefer k.deinit(gpa);

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

        ffft: FftPlan(Complex(T)),
        ifft: FftPlan(Complex(T)),

        pub fn init(ndim: usize, nstate: usize, npoint: usize, plan_mode: u32, gpa: Allocator) !@This() {
            var W = try Matrix(Complex(T)).init(nstate, std.math.pow(usize, npoint, ndim), gpa);
            errdefer W.deinit(gpa);

            const shape = try gpa.alloc(i32, ndim);
            defer gpa.free(shape);

            for (shape) |*e| {
                e.* = @as(i32, @intCast(npoint));
            }

            const ffft = try FftPlan(Complex(T)).init(W.rowSlice(0), shape, -1, plan_mode);
            errdefer ffft.deinit();

            const ifft = try FftPlan(Complex(T)).init(W.rowSlice(0), shape, 1, plan_mode);
            errdefer ifft.deinit();

            return .{ .W = W, .ffft = ffft, .ifft = ifft };
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.W.deinit(gpa);

            self.ffft.deinit();
            self.ifft.deinit();
        }

        pub fn clone(self: @This(), gpa: Allocator) !@This() {
            const ffft = try self.ffft.clone();
            errdefer ffft.deinit();

            const ifft = try self.ifft.clone();
            errdefer ifft.deinit();

            const W = try self.W.clone(gpa);
            errdefer W.deinit(gpa);

            return .{ .W = W, .ffft = ffft, .ifft = ifft };
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

                if (comptime sign == -1) self.ffft.execute(slice) else self.ifft.execute(slice);
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

// HAMILTONIAN =========================================================================================================

fn Hamiltonian(comptime T: type) type {
    return struct {
        V: Matrix(T),
        W: Matrix(T),
        U: Matrix(T),
        K: Vector(T),

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
            var R = try Matrix(Complex(T)).init(ham.V.nrow(), ham.V.ncol(), gpa);
            errdefer R.deinit(gpa);

            var K = try Vector(Complex(T)).initZero(ham.K.length(), gpa);
            errdefer K.deinit(gpa);

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
            try self.applyR(wfn, gpa);

            try self.applyK(wfn);

            try self.applyR(wfn, gpa);
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

        pub fn init(sim: *SimulationState(T), write: Write, adia: bool, log: bool, gpa: Allocator) !@This() {
            var obs = @This(){};
            errdefer obs.deinit(gpa);

            const calc_ekin, const calc_epot = .{ write.kinetic_energy != null, write.potential_energy != null };

            const calc_norm = write.norm != null;

            var calc = .{
                .pos = log or write.position != null,
                .mom = log or write.momentum != null,

                .pop = log or write.population != null,

                .norm = log or calc_norm,
                .ekin = log or calc_ekin,
                .epot = log or calc_epot,
            };

            calc.ekin = calc.ekin or write.total_energy != null;
            calc.epot = calc.epot or write.total_energy != null;

            if (calc.pos) {
                obs.pos = try sim.wfn.pos(sim.wfn_kpgrids, gpa);
            }

            if (calc.norm) {
                obs.norm = sim.wfn.norm(sim.wfn_kpgrids);
            }

            if (calc.epot) {
                obs.epot = sim.wfn.epot(sim.hams, sim.wfn_kpgrids);
            }

            if (calc.pop) {
                if (adia == true) {
                    obs.pop = try sim.wfn.popAdia(sim.hams, sim.wfn_kpgrids, gpa);
                }

                if (adia != true) {
                    obs.pop = try sim.wfn.pop(sim.wfn_kpgrids, gpa);
                }
            }

            const needs_fft = calc.mom or calc.ekin;

            if (needs_fft) {
                try sim.wfn.fft(-1);

                if (calc.mom) {
                    obs.mom = try sim.wfn.mom(sim.wfn_kpgrids, gpa);
                }

                if (calc.ekin) {
                    obs.ekin = sim.wfn.ekin(sim.hams, sim.wfn_kpgrids);
                }

                try sim.wfn.fft(1);
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
        pos: ?Matrix(T) = null,
        mom: ?Matrix(T) = null,
        pop: ?Matrix(T) = null,

        epot: ?Matrix(T) = null,
        ekin: ?Matrix(T) = null,
        etot: ?Matrix(T) = null,
        norm: ?Matrix(T) = null,

        wfn: ?Matrix(T) = null,

        index: usize = 0,

        pub fn init(ndim: usize, nstate: usize, npoint: usize, iters: usize, write: Write, gpa: Allocator) !@This() {
            var hist = @This(){};
            errdefer hist.deinit(gpa);

            var store_ekin, var store_epot = .{ write.kinetic_energy != null, write.potential_energy != null };

            const store_wfn = write.wavefunction != null;

            store_epot = store_epot or write.total_energy != null;
            store_ekin = store_ekin or write.total_energy != null;

            const store_etot = write.total_energy != null;

            const wfn_nrow = std.math.pow(usize, npoint, ndim);

            if (write.position != null) hist.pos = try Matrix(T).init(iters, ndim, gpa);
            if (write.momentum != null) hist.mom = try Matrix(T).init(iters, ndim, gpa);

            if (write.population != null) {
                hist.pop = try Matrix(T).init(iters, nstate, gpa);
            }

            if (write.norm != null) {
                hist.norm = try Matrix(T).init(iters, 1, gpa);
            }

            if (store_wfn) {
                hist.wfn = try Matrix(T).init(wfn_nrow, 2 * nstate * iters, gpa);
            }

            if (store_etot) hist.etot = try Matrix(T).init(iters, 1, gpa);
            if (store_ekin) hist.ekin = try Matrix(T).init(iters, 1, gpa);
            if (store_epot) hist.epot = try Matrix(T).init(iters, 1, gpa);

            return hist;
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            if (self.pos) |*pos| pos.deinit(gpa);
            if (self.mom) |*mom| mom.deinit(gpa);
            if (self.pop) |*pop| pop.deinit(gpa);

            if (self.epot) |*epot| epot.deinit(gpa);
            if (self.ekin) |*ekin| ekin.deinit(gpa);
            if (self.etot) |*etot| etot.deinit(gpa);
            if (self.norm) |*norm| norm.deinit(gpa);

            if (self.wfn) |*wfn| wfn.deinit(gpa);
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

            if (self.epot) |*epot| {
                if (obs.epot) |v| epot.ptr(step_idx, 0).* = v;
            }

            if (self.ekin) |*ekin| {
                if (obs.ekin) |v| ekin.ptr(step_idx, 0).* = v;
            }

            if (self.norm) |*norm| {
                if (obs.norm) |v| norm.ptr(step_idx, 0).* = v;
            }

            if (self.etot) |*etot| {
                etot.ptr(step_idx, 0).* = obs.ekin.? + obs.epot.?;
            }

            self.index += 1;
        }

        pub fn exportWrite(self: *@This(), io: std.Io, dt: f64, grid: Grid(T), write: Write) !void {
            const end = dt * @as(T, @floatFromInt(self.index - 1));

            if (write.position) |path| {
                try writeMatrixLspace(T, io, path, self.pos.?, 0, end);
            }

            if (write.momentum) |path| {
                try writeMatrixLspace(T, io, path, self.mom.?, 0, end);
            }

            if (write.population) |path| {
                try writeMatrixLspace(T, io, path, self.pop.?, 0, end);
            }

            if (write.potential_energy) |path| {
                try writeMatrixLspace(T, io, path, self.epot.?, 0, end);
            }

            if (write.kinetic_energy) |path| {
                try writeMatrixLspace(T, io, path, self.ekin.?, 0, end);
            }

            if (write.norm) |path| {
                try writeMatrixLspace(T, io, path, self.norm.?, 0, end);
            }

            if (write.total_energy) |path| {
                try writeMatrixLspace(T, io, path, self.etot.?, 0, end);
            }

            if (write.wavefunction) |path| {
                try writeMatrixHjoin(T, io, path, grid.r, self.wfn.?);
            }
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

        "POS (a0)",
        12 * ndim,

        "MOM (hb/a0)",
        12 * ndim,

        "POP (-)",
        11 * nstate,

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

// RUN =================================================================================================================

fn SimulationState(comptime T: type) type {
    return struct {
        wfn_kpgrids: Grid(T),
        hams: Hamiltonian(T),
        epoten: Potential(T),
        wfn: Wavefunction(T),
        propg: Propagator(T),

        orthw: std.ArrayList(Wavefunction(T)),

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            for (self.orthw.items) |*e| e.deinit(gpa);

            inline for (@typeInfo(@This()).@"struct".fields) |field| {
                @field(self, field.name).deinit(gpa);
            }
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
        .estimate => fftw.FFTW_ESTIMATE,
        .measure => fftw.FFTW_MEASURE,
        .patient => fftw.FFTW_PATIENT,
        .exhaustive => fftw.FFTW_EXHAUSTIVE,
    };

    var grid = try Grid(T).init(opt.grid.bounds, opt.grid.npoint, gpa);
    errdefer grid.deinit(gpa);

    var wfn = try Wavefunction(T).init(pot.ndim(), pot.nstate(), opt.grid.npoint, plan_mode, gpa);
    errdefer wfn.deinit(gpa);

    var ham = try Hamiltonian(T).init(grid, pot, opt.mass, gpa);
    errdefer ham.deinit(gpa);

    const prop = try Propagator(T).init(ham, dt, gpa);
    errdefer prop.deinit(gpa);

    return .{ .wfn_kpgrids = grid, .hams = ham, .wfn = wfn, .epoten = pot, .propg = prop, .orthw = .empty };
}

fn solve(comptime T: type, io: std.Io, ctx: SolveContext(T), gpa: Allocator, arena: Allocator) !Observables(T) {
    const ndim, const nstate = .{ ctx.sim.epoten.ndim(), ctx.sim.epoten.nstate() };

    const neig = if (ctx.opt.imaginary) |imag| imag.nstate else 1;

    if (ctx.log) try printHeader(io, ctx.eigs, ndim, nstate, neig);

    ctx.sim.wfn.setGaussian(ctx.opt.initial_conditions, ctx.sim.wfn_kpgrids);

    var hist = try History(T).init(ndim, nstate, ctx.opt.grid.npoint, ctx.opt.iterations + 1, ctx.opt.write, gpa);
    defer hist.deinit(gpa);

    if (ctx.opt.initial_conditions.adiabatic) {
        try ctx.sim.wfn.toDia(ctx.sim.hams, gpa);
    }

    var timer = std.Io.Timestamp.now(io, .real);

    for (0..ctx.opt.iterations + 1) |i| {
        const time = (@as(T, @floatFromInt(i)) - 0.5) * ctx.opt.time_step;

        if (i > 0 and ctx.sim.epoten.isTd()) {
            try ctx.sim.hams.update(ctx.sim.wfn_kpgrids, ctx.sim.epoten, time);

            ctx.sim.propg.update(ctx.sim.hams);
        }

        if (i > 0) try ctx.sim.propg.step(&ctx.sim.wfn, gpa);

        if (ctx.opt.imaginary != null) for (ctx.sim.orthw.items) |e| {
            const overlap = e.overlap(ctx.sim.wfn, ctx.sim.wfn_kpgrids);

            for (0..ctx.sim.wfn.W.data.len) |k| {
                ctx.sim.wfn.W.data[k] = ctx.sim.wfn.W.data[k].sub(overlap.mul(e.W.data[k]));
            }
        };

        if (ctx.opt.imaginary != null) ctx.sim.wfn.normalize(ctx.sim.wfn_kpgrids);

        const is_log_step = ctx.log and ((i % ctx.opt.log_interval == 0) or (i == ctx.opt.iterations));

        if (ctx.sim.epoten.isTd()) {
            const t = @as(T, @floatFromInt(i)) * ctx.opt.time_step;

            try ctx.sim.hams.update(ctx.sim.wfn_kpgrids, ctx.sim.epoten, t);
        }

        var obs = try Observables(T).init(ctx.sim, ctx.opt.write, ctx.opt.adiabatic, is_log_step, gpa);
        defer obs.deinit(gpa);

        if (ctx.opt.write.wavefunction != null and ctx.opt.adiabatic) {
            try ctx.sim.wfn.toAdia(ctx.sim.hams, gpa);
        }

        hist.append(ctx.sim.wfn, obs);

        if (ctx.opt.write.wavefunction != null and ctx.opt.adiabatic) {
            try ctx.sim.wfn.toDia(ctx.sim.hams, gpa);
        }

        if (is_log_step) {
            try printIteration(T, io, obs, i, &timer);
        }
    }

    try hist.exportWrite(io, ctx.opt.time_step, ctx.sim.wfn_kpgrids, ctx.opt.write);

    return try Observables(T).init(ctx.sim, ctx.opt.write, ctx.opt.adiabatic, true, arena);
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
            var cloned = try sim.wfn.clone(gpa);
            errdefer cloned.deinit(gpa);

            try sim.orthw.append(gpa, cloned);
        }

        try output.append(arena, obs);

        if (log) {
            try printFinalPop(T, io, obs);
        }
    }

    if (log) try printFinalEnergies(T, io, output);

    return output;
}
