const std = @import("std");

const fftw = @import("cimport.zig").fftw;

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;

const Grid = @import("wavepacket.zig").Grid;
const Hamiltonian = @import("wavepacket.zig").Hamiltonian;
const InitialConditions = @import("wavepacket.zig").InitialConditions;
const Matrix = @import("tensor.zig").Matrix;
const Potential = @import("potential.zig").Potential;
const PotentialOptions = @import("potential.zig").Options;
const Vector = @import("tensor.zig").Vector;
const Wavefunction = @import("wavepacket.zig").Wavefunction;

const calcSpectrum = @import("spectral_analysis.zig").calcSpectrum;
const printf = @import("read_write.zig").printf;
const writeMatrixHjoin = @import("read_write.zig").writeMatrixHjoin;
const writeMatrixLspace = @import("read_write.zig").writeMatrixLspace;

pub const Options = struct {
    initial_conditions: InitialConditions,
    potential: PotentialOptions,

    time_step: f64,
    iterations: u32,
    mass: f64,

    write: Write = .{},

    adiabatic: bool = false,
    log_interval: u32 = 1,

    absorbing_potential: ?struct {
        track_population: bool = false,
        bounds: []const [2]f64,
        exponent: f64 = 0.001,
        stop_norm: ?f64 = null,
    } = null,

    fft: struct {
        plan: enum { estimate, measure, patient, exhaustive } = .measure,
    } = .{},

    grid: struct {
        bounds: []const [2]f64,
        npoint: u32,
    },

    imaginary: ?struct {
        nstate: u32 = 1,
    } = null,

    spectrum: ?struct {
        padding: u32 = 0,
        threshold: f64 = 1e-6,

        write: struct {
            acf: ?[]const u8 = null,
            spectrum: ?[]const u8 = null,
        } = .{},
    } = null,
};

pub fn Result(comptime T: type) type {
    return struct {
        observables: std.ArrayList(Observables(T)),

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            for (0..self.observables.items.len) |i| {
                self.observables.items[i].deinit(gpa);
            }

            self.observables.deinit(gpa);
        }
    };
}

const Write = struct {
    acf: ?[]const u8 = null,
    kinetic_energy: ?[]const u8 = null,
    momentum: ?[]const u8 = null,
    norm: ?[]const u8 = null,
    population: ?[]const u8 = null,
    position: ?[]const u8 = null,
    potential_energy: ?[]const u8 = null,
    total_energy: ?[]const u8 = null,
    wavefunction: ?[]const u8 = null,
};

fn History(comptime T: type) type {
    return struct {
        acf: ?Vector(Complex(T)) = null,

        pos: ?Matrix(T) = null,
        mom: ?Matrix(T) = null,
        pop: ?Matrix(T) = null,

        epot: ?Matrix(T) = null,
        ekin: ?Matrix(T) = null,
        etot: ?Matrix(T) = null,
        norm: ?Matrix(T) = null,

        wfn: ?Matrix(T) = null,

        index: usize = 0,

        pub fn init(ndim: usize, nstate: usize, npoint: usize, iters: usize, write: Write, spectrum: bool, gpa: Allocator) !@This() {
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

            if (spectrum) {
                hist.acf = try Vector(Complex(T)).init(iters, gpa);
            }

            return hist;
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            if (self.acf) |*acf| acf.deinit(gpa);
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

            if (self.acf) |*acf| {
                acf.ptr(step_idx).* = obs.acf.?;
            }

            self.index += 1;
        }

        pub fn exportWrite(self: *@This(), io: std.Io, dt: T, grid: Grid(T), write: Write, spectrum: anytype, nstate: usize, gpa: Allocator) !void {
            const end = dt * @as(T, @floatFromInt(self.index - 1));

            if (write.acf) |path| {
                try writeMatrixLspace(Complex(T), io, path, self.acf.?.takeRows(self.index).asMatrix(), 0, end);
            }

            if (write.position) |path| {
                try writeMatrixLspace(T, io, path, self.pos.?.takeRows(self.index), 0, end);
            }

            if (write.momentum) |path| {
                try writeMatrixLspace(T, io, path, self.mom.?.takeRows(self.index), 0, end);
            }

            if (write.population) |path| {
                try writeMatrixLspace(T, io, path, self.pop.?.takeRows(self.index), 0, end);
            }

            if (write.potential_energy) |path| {
                try writeMatrixLspace(T, io, path, self.epot.?.takeRows(self.index), 0, end);
            }

            if (write.kinetic_energy) |path| {
                try writeMatrixLspace(T, io, path, self.ekin.?.takeRows(self.index), 0, end);
            }

            if (write.norm) |path| {
                try writeMatrixLspace(T, io, path, self.norm.?.takeRows(self.index), 0, end);
            }

            if (write.total_energy) |path| {
                try writeMatrixLspace(T, io, path, self.etot.?.takeRows(self.index), 0, end);
            }

            if (write.wavefunction) |path| {
                try writeMatrixHjoin(T, io, path, grid.r, null, self.wfn.?, 2 * nstate * self.index);
            }

            if (spectrum) |spec| {
                var acfpd, var sigma = try calcSpectrum(T, self.acf.?.takeRows(self.index), dt, spec.padding, spec.threshold, gpa);

                defer {
                    acfpd.deinit(gpa);
                    sigma.deinit(gpa);
                }

                if (spec.write.acf) |path| {
                    const pad = dt * @as(T, @floatFromInt(spec.padding));

                    try writeMatrixLspace(Complex(T), io, path, acfpd.asMatrix(), 0, end + pad);
                }

                if (spec.write.spectrum) |path| {
                    const nyquist, const nt = .{ std.math.pi / dt, @as(T, @floatFromInt(sigma.length())) };

                    try writeMatrixLspace(T, io, path, sigma.asMatrix(), -nyquist, nyquist * (nt - 2) / nt);
                }
            }
        }
    };
}

fn Observables(comptime T: type) type {
    return struct {
        pos: ?Vector(T) = null,
        mom: ?Vector(T) = null,
        pop: ?Vector(T) = null,

        acf: ?Complex(T) = null,

        epot: ?T = null,
        ekin: ?T = null,
        norm: ?T = null,

        pub fn init(sim: *SimulationState(T), wfn0: ?Wavefunction(T), write: Write, adia: bool, log: bool, gpa: Allocator) !@This() {
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

            if (wfn0 != null) {
                obs.acf = wfn0.?.overlap(sim.wfn, sim.wfn_kpgrids);
            }

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

                if (obs.pop) |*pop| for (0..pop.length()) |j| {
                    pop.ptr(j).* += sim.pop_apabs.at(j);
                };
            }

            const needs_fft = calc.mom or calc.ekin;

            if (needs_fft) {
                sim.wfn.fft(-1);

                defer {
                    sim.wfn.fft(1);
                }

                if (calc.mom) {
                    obs.mom = try sim.wfn.mom(sim.wfn_kpgrids, gpa);
                }

                if (calc.ekin) {
                    obs.ekin = sim.wfn.ekin(sim.hams, sim.wfn_kpgrids);
                }
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

fn Propagator(comptime T: type) type {
    return struct {
        R: Matrix(Complex(T)),
        K: Vector(Complex(T)),
        cap_weight: Vector(T),

        dt: Complex(T),

        pub fn init(grid: Grid(T), ham: Hamiltonian(T), capopt: anytype, dt: Complex(T), gpa: Allocator) !@This() {
            var R = try Matrix(Complex(T)).init(ham.V.nrow(), ham.V.ncol(), gpa);
            errdefer R.deinit(gpa);

            var K = try Vector(Complex(T)).initZero(ham.K.length(), gpa);
            errdefer K.deinit(gpa);

            var cap_weight = try Vector(T).initZero(ham.K.length(), gpa);
            errdefer cap_weight.deinit(gpa);

            for (0..ham.K.data.len) |i| {
                K.data[i] = std.math.complex.exp(Complex(T).init(0, -ham.K.data[i]).mul(dt));
            }

            var prop = @This(){ .R = R, .K = K, .cap_weight = cap_weight, .dt = dt };

            prop.update(grid, ham, capopt);

            return prop;
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.R.deinit(gpa);
            self.K.deinit(gpa);

            self.cap_weight.deinit(gpa);
        }

        pub fn step(self: @This(), sim: *SimulationState(T), adia: bool, track_pop: bool, gpa: Allocator) !void {
            if (track_pop) self.accumAbsorbed(sim, adia);

            try self.applyR(&sim.wfn, gpa);

            try self.applyK(&sim.wfn);

            if (track_pop) self.accumAbsorbed(sim, adia);

            try self.applyR(&sim.wfn, gpa);
        }

        pub fn update(self: *@This(), grid: Grid(T), ham: Hamiltonian(T), capopt: anytype) void {
            const nstate = std.math.sqrt(ham.V.ncol());

            for (0..self.R.nrow()) |i| {
                var cap_sum: T = 0;

                if (capopt) |cap| for (0..grid.r.ncol()) |j| {
                    const r_val = grid.r.at(i, j);

                    const min_bound = cap.bounds[j][0];
                    const max_bound = cap.bounds[j][1];

                    const dist = @max(0, @max(min_bound - r_val, r_val - max_bound));

                    if (dist > 0) {
                        cap_sum += std.math.exp(cap.exponent * dist) - 1;
                    }
                };

                const cap_decay = std.math.exp(-0.5 * cap_sum * self.dt.magnitude());

                self.cap_weight.ptr(i).* = 1 - cap_decay * cap_decay;

                for (0..nstate) |k| for (0..nstate) |j| {
                    var sum = Complex(T).init(0, 0);

                    for (0..nstate) |m| {
                        const U_km = ham.U.at(i, k * nstate + m);
                        const U_jm = ham.U.at(i, j * nstate + m);

                        const phase = std.math.complex.exp(Complex(T).init(0, -0.5 * ham.W.at(i, m)).mul(self.dt));

                        sum = sum.add(Complex(T).init(phase.re * U_km * U_jm, phase.im * U_km * U_jm));
                    }

                    self.R.ptr(i, k * nstate + j).* = sum.mul(Complex(T).init(cap_decay, 0));
                };
            }
        }

        fn accumAbsorbed(self: @This(), sim: *SimulationState(T), adia: bool) void {
            for (0..sim.wfn.W.nrow()) |i| {
                var sum: T = 0;

                for (0..sim.wfn.W.ncol()) |j| {
                    const w = self.cap_weight.at(j);

                    if (adia) {
                        var adia_re: T = 0;
                        var adia_im: T = 0;

                        for (0..sim.wfn.W.nrow()) |k| {
                            const u = sim.hams.U.at(j, k * sim.wfn.W.nrow() + i);

                            adia_re += sim.wfn.W.at(k, j).re * u;
                            adia_im += sim.wfn.W.at(k, j).im * u;
                        }

                        sum += w * (adia_re * adia_re + adia_im * adia_im);
                    }

                    if (!adia) {
                        sum += w * sim.wfn.W.at(i, j).squaredMagnitude();
                    }
                }

                sim.pop_apabs.ptr(i).* += sum * sim.wfn_kpgrids.dr;
            }
        }

        fn applyK(self: @This(), wfn: *Wavefunction(T)) !void {
            wfn.fft(-1);

            defer {
                wfn.fft(1);
            }

            for (0..wfn.W.nrow()) |i| for (0..wfn.W.ncol()) |j| {
                wfn.W.ptr(i, j).* = self.K.at(j).mul(wfn.W.at(i, j));
            };
        }

        fn applyR(self: @This(), wfn: *Wavefunction(T), gpa: Allocator) !void {
            var temp = try gpa.alloc(Complex(T), wfn.W.nrow());
            defer gpa.free(temp);

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
        }
    };
}

fn SimulationState(comptime T: type) type {
    return struct {
        wfn_kpgrids: Grid(T),
        hams: Hamiltonian(T),
        epoten: Potential(T),
        wfn: Wavefunction(T),
        propg: Propagator(T),
        pop_apabs: Vector(T),

        orthw: std.ArrayList(Wavefunction(T)),

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            for (0..self.orthw.items.len) |i| self.orthw.items[i].deinit(gpa);

            inline for (@typeInfo(@This()).@"struct".fields) |field| {
                @field(self, field.name).deinit(gpa);
            }
        }
    };
}

fn SolveContext(comptime T: type) type {
    return struct { opt: Options, sim: *SimulationState(T), eigs: usize, log: bool };
}

pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator) !Result(T) {
    try checkInvalidInput(opt);

    var result: Result(T) = .{ .observables = .empty };
    errdefer result.deinit(gpa);

    if (log) try std.Io.File.stdout().writeStreamingAll(io, "\nQUANTUM DYNAMICS INIT: ");

    var timer = std.Io.Timestamp.now(io, .real);

    var sim = try init(T, io, opt, gpa);
    defer sim.deinit(gpa);

    if (log) try printf(io, "{f}\n", .{timer.untilNow(io, .real)});

    for (0..if (opt.imaginary) |imag| imag.nstate else 1) |i| {
        var obs = try solve(T, io, .{ .opt = opt, .sim = &sim, .eigs = i, .log = log }, gpa);
        errdefer obs.deinit(gpa);

        if (i < if (opt.imaginary) |imag| imag.nstate else 0) {
            var cloned = try sim.wfn.clone(gpa);
            errdefer cloned.deinit(gpa);

            try sim.orthw.append(gpa, cloned);
        }

        try result.observables.append(gpa, obs);

        if (log) {
            try printFinalPop(T, io, obs);
        }
    }

    if (log) try printFinalEnergies(T, io, result.observables);

    return result;
}

fn checkInvalidInput(opt: Options) !void {
    if (opt.time_step <= 0) {
        std.log.err("TIME STEP MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.mass <= 0) {
        std.log.err("MASS MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.log_interval == 0) {
        std.log.err("LOG INTERVAL MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.grid.npoint == 0) {
        std.log.err("GRID POINT COUNT MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.grid.bounds.len == 0) {
        std.log.err("GRID BOUNDS MUST NOT BE EMPTY", .{});

        return error.InvalidInput;
    }

    for (opt.grid.bounds) |b| if (b[0] >= b[1]) {
        std.log.err("GRID BOUNDS MIN MUST BE LESS THAN MAX", .{});

        return error.InvalidInput;
    };

    if (opt.initial_conditions.position.len == 0) {
        std.log.err("INITIAL POSITION VECTOR MUST NOT BE EMPTY", .{});

        return error.InvalidInput;
    }

    if (opt.initial_conditions.momentum.len != opt.initial_conditions.position.len) {
        std.log.err("INITIAL MOMENTUM AND POSITION VECTORS MUST HAVE THE SAME LENGTH", .{});

        return error.InvalidInput;
    }

    if (opt.initial_conditions.gamma.len != opt.initial_conditions.position.len) {
        std.log.err("INITIAL GAMMA VECTOR MUST HAVE THE SAME LENGTH AS POSITION VECTOR", .{});

        return error.InvalidInput;
    }

    if (opt.initial_conditions.position.len != opt.grid.bounds.len) {
        std.log.err("INITIAL CONDITIONS DIMENSION DOES NOT MATCH GRID DIMENSION", .{});

        return error.InvalidInput;
    }

    if (opt.absorbing_potential) |cap| {
        if (cap.bounds.len != opt.grid.bounds.len) {
            std.log.err("ABSORBING POTENTIAL BOUNDS DIMENSION DOES NOT MATCH GRID DIMENSION", .{});

            return error.InvalidInput;
        }

        for (cap.bounds, opt.grid.bounds) |cb, gb| if (cb[0] < gb[0] or cb[1] > gb[1] or cb[0] >= cb[1]) {
            std.log.err("ABSORBING POTENTIAL BOUNDS MUST LIE WITHIN GRID BOUNDS AND CB[0] < CB[1]", .{});

            return error.InvalidInput;
        };

        if (cap.exponent <= 0) {
            std.log.err("ABSORBING POTENTIAL EXPONENT MUST BE GREATER THAN 0", .{});

            return error.InvalidInput;
        }

        if (cap.stop_norm) |sn| if (sn <= 0 or sn >= 1) {
            std.log.err("ABSORBING POTENTIAL STOP_NORM MUST BE BETWEEN 0 AND 1", .{});

            return error.InvalidInput;
        };
    }

    if (opt.imaginary) |imag| if (imag.nstate == 0) {
        std.log.err("ITP TARGET STATE COUNT MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    };

    if (opt.spectrum) |spec| if (spec.threshold <= 0 or spec.threshold >= 1) {
        std.log.err("SPECTRUM THRESHOLD MUST BE BETWEEN 0 AND 1", .{});

        return error.InvalidInput;
    };
}

fn init(comptime T: type, io: std.Io, opt: Options, gpa: Allocator) !SimulationState(T) {
    const pot = try Potential(T).init(io, opt.potential, gpa);

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

    var prop = try Propagator(T).init(grid, ham, opt.absorbing_potential, dt, gpa);
    errdefer prop.deinit(gpa);

    var pop_apabs = try Vector(T).initZero(pot.nstate(), gpa);
    errdefer pop_apabs.deinit(gpa);

    return .{ .wfn_kpgrids = grid, .hams = ham, .wfn = wfn, .epoten = pot, .propg = prop, .pop_apabs = pop_apabs, .orthw = .empty };
}

fn printFinalEnergies(comptime T: type, io: std.Io, obs: std.ArrayList(Observables(T))) !void {
    try std.Io.File.stdout().writeStreamingAll(io, "\n");

    for (0..obs.items.len) |i| {
        const ekin = obs.items[i].ekin orelse std.math.nan(T);
        const epot = obs.items[i].epot orelse std.math.nan(T);

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

fn printHeader(io: std.Io, eigs: usize, ndim: usize, nstate: usize, neig: usize) !void {
    if (neig > 1) {
        try printf(io, "\nITP OF STATE {d}/{d}", .{ eigs + 1, neig });
    }

    if (neig == 1) {
        try printf(io, "\nREAL-TIME PROPAGATION", .{});
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

fn solve(comptime T: type, io: std.Io, ctx: SolveContext(T), gpa: Allocator) !Observables(T) {
    const ndim, const nstate = .{ ctx.sim.epoten.ndim(), ctx.sim.epoten.nstate() };

    const neig = if (ctx.opt.imaginary) |imag| imag.nstate else 1;

    if (ctx.log) try printHeader(io, ctx.eigs, ndim, nstate, neig);

    ctx.sim.wfn.setGaussian(ctx.opt.initial_conditions, ctx.sim.wfn_kpgrids);

    ctx.sim.pop_apabs.zero();

    const calc_acf = ctx.opt.spectrum != null or ctx.opt.write.acf != null;

    var hist = try History(T).init(ndim, nstate, ctx.opt.grid.npoint, ctx.opt.iterations + 1, ctx.opt.write, calc_acf, gpa);
    defer hist.deinit(gpa);

    if (ctx.opt.initial_conditions.adiabatic) {
        try ctx.sim.wfn.toDia(ctx.sim.hams, gpa);
    }

    const track_cap_pop = ctx.opt.absorbing_potential != null and ctx.opt.absorbing_potential.?.track_population;

    var wfn0 = if (calc_acf) try ctx.sim.wfn.clone(gpa) else null;
    defer if (wfn0) |*w| w.deinit(gpa);

    var timer = std.Io.Timestamp.now(io, .real);

    for (0..ctx.opt.iterations + 1) |i| {
        const time = (@as(T, @floatFromInt(i)) - 0.5) * ctx.opt.time_step;

        if (i > 0 and ctx.sim.epoten.isTd()) {
            try ctx.sim.hams.update(ctx.sim.wfn_kpgrids, ctx.sim.epoten, time, gpa);

            ctx.sim.propg.update(ctx.sim.wfn_kpgrids, ctx.sim.hams, ctx.opt.absorbing_potential);
        }

        if (i > 0) try ctx.sim.propg.step(ctx.sim, ctx.opt.adiabatic, track_cap_pop, gpa);

        if (ctx.opt.imaginary != null) for (0..ctx.sim.orthw.items.len) |j| {
            const overlap = ctx.sim.orthw.items[j].overlap(ctx.sim.wfn, ctx.sim.wfn_kpgrids);

            for (0..ctx.sim.wfn.W.data.len) |k| {
                ctx.sim.wfn.W.data[k] = ctx.sim.wfn.W.data[k].sub(overlap.mul(ctx.sim.orthw.items[j].W.data[k]));
            }
        };

        if (ctx.opt.imaginary != null) ctx.sim.wfn.normalize(ctx.sim.wfn_kpgrids);

        const is_log_step = ctx.log and ((i % ctx.opt.log_interval == 0) or (i == ctx.opt.iterations));

        if (ctx.sim.epoten.isTd()) {
            const t = @as(T, @floatFromInt(i)) * ctx.opt.time_step;

            try ctx.sim.hams.update(ctx.sim.wfn_kpgrids, ctx.sim.epoten, t, gpa);
        }

        var obs = try Observables(T).init(ctx.sim, wfn0, ctx.opt.write, ctx.opt.adiabatic, is_log_step, gpa);
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

        if (ctx.opt.absorbing_potential) |ap| if (ap.stop_norm) |stop_norm| {
            const norm = ctx.sim.wfn.norm(ctx.sim.wfn_kpgrids);

            if (norm < stop_norm) {
                var stop_obs = try Observables(T).init(ctx.sim, wfn0, ctx.opt.write, ctx.opt.adiabatic, true, gpa);
                defer stop_obs.deinit(gpa);

                if (!is_log_step) try printIteration(T, io, stop_obs, i, &timer);

                break;
            }
        };
    }

    try hist.exportWrite(io, ctx.opt.time_step, ctx.sim.wfn_kpgrids, ctx.opt.write, ctx.opt.spectrum, nstate, gpa);

    return try Observables(T).init(ctx.sim, wfn0, ctx.opt.write, ctx.opt.adiabatic, true, gpa);
}
