//! Grid-based quantum dynamics propagation of molecular wavefunctions on coupled potential energy surfaces.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;

const partial_wave = @import("partial_wave.zig");
const quantum_observables = @import("quantum_observables.zig");
const split_operator = @import("split_operator.zig");

const FluxAnalysis = @import("flux_analysis.zig").FluxAnalysis;
const Grid = @import("wavepacket.zig").Grid;
const Hamiltonian = @import("wavepacket.zig").Hamiltonian;
const InitialConditions = @import("wavepacket.zig").InitialConditions;
const Matrix = @import("tensor.zig").Matrix;
const Potential = @import("potential.zig").Potential;
const PotentialOptions = @import("potential.zig").Options;
const Vector = @import("tensor.zig").Vector;
const Wavefunction = @import("wavepacket.zig").Wavefunction;
const validation = @import("validation.zig");

const calcSpectrum = @import("spectral_analysis.zig").calcSpectrum;
const geru = @import("linear_algebra.zig").geru;
const printf = @import("read_write.zig").printf;
const writeMatrixHjoin = @import("read_write.zig").writeMatrixHjoin;
const writeMatrixLinspace = @import("read_write.zig").writeMatrixLinspace;
const writeMatrixLogspace = @import("read_write.zig").writeMatrixLogspace;

pub const estimateMaxJ = partial_wave.estimateMaxJ;
pub const injectAngularFname = partial_wave.injectAngularFname;
pub const Observables = quantum_observables.Observables;
pub const PartialWaveContext = partial_wave.PartialWaveContext;
pub const PartialWaveOptions = partial_wave.PartialWaveOptions;
pub const PartialWaveWrite = partial_wave.PartialWaveWrite;
pub const Propagator = split_operator.SplitOperator;
pub const SplitOperator = split_operator.SplitOperator;

/// Configuration options for split-operator quantum dynamics wavepacket propagation on a grid.
pub const Options = struct {
    absorbing_potential: ?AbsorbingPotential = null,
    adiabatic: bool = false,
    flux_analysis: ?FluxAnalysisOptions = null,
    grid: GridOptions,
    imaginary: ?Imaginary = null,
    initial_conditions: InitialConditions,
    iterations: u32,
    j_quantum_number: u32 = 0,
    log_interval: u32 = 1,
    mass: []const f64,
    memory: Memory = .{},
    partial_waves: ?PartialWaveOptions = null,
    potential: PotentialOptions,
    spectrum: ?SpectrumOptions = null,
    time_step: f64,
    write: Write = .{},
};

/// Negative imaginary potential parameters for absorbing outgoing boundary wavepackets.
const AbsorbingPotential = struct {
    bounds: []const [2]f64,
    exponent: f64 = 0.001,
    stop_norm: ?f64 = null,
    track_population: bool = false,
};

/// Quantum flux analysis settings for computing energy-resolved transmission cross sections.
const FluxAnalysisOptions = struct {
    asymptotic_offset: f64 = 0,
    e_max: f64,
    e_min: f64,
    flux_bounds: []const [2]f64,
    logarithmic: bool = false,
    npoint: usize,
    write: FluxAnalysisWrite = .{},
};

/// Output file paths for recording reaction flux cross sections.
const FluxAnalysisWrite = struct {
    cross_section: ?[]const u8 = null,
};

/// Coordinate grid bounds and point discretization for wavefunction representation.
const GridOptions = struct {
    bounds: []const [2]f64,
    cylindrical: bool = false,
    npoint: u32,
};

/// Settings for imaginary time propagation to find the lowest energy eigenstates.
const Imaginary = struct {
    nstate: u32 = 1,
};

/// Memory caching toggles for grid coordinates, potential surfaces, and propagators.
const Memory = struct {
    grid: bool = true,
    potential: bool = true,
    propagator: bool = true,
};

/// Spectral analysis parameters for Fourier transforming wavefunction autocorrelation functions.
const SpectrumOptions = struct {
    padding: u32 = 0,
    threshold: f64 = 1e-6,
    write: SpectrumWrite = .{},
};

/// Output file paths for recording autocorrelation functions and absorption spectra.
const SpectrumWrite = struct {
    acf: ?[]const u8 = null,
    spectrum: ?[]const u8 = null,
};

/// File paths for exporting time-dependent wavefunctions and expectation values to disk.
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

/// Stores accumulated quantum observables and cross sections from wavepacket propagation.
pub fn Result(comptime T: type) type {
    return struct {
        observables: std.ArrayList(Observables(T)),

        cross_section: ?Matrix(T) = null,

        /// Deallocates memory associated with the quantum dynamics results.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            for (0..self.observables.items.len) |i| {
                self.observables.items[i].deinit(gpa);
            }

            self.observables.deinit(gpa);

            if (self.cross_section) |*cs| {
                cs.deinit(gpa);
            }
        }
    };
}

/// Accumulates time-dependent wavefunctions, autocorrelation functions, and expectation values.
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

        flux_acc: ?Matrix(Complex(T)) = null,
        flux_phs: ?Vector(Complex(T)) = null,
        wfn_init: ?Matrix(Complex(T)) = null,
        e_thrs: T = 0,

        index: usize = 0,

        /// Allocates memory for wavefunctions and observables history arrays.
        pub fn init(ndim: usize, nstate: usize, npoint: usize, iters: usize, write: Write, spectrum: bool, opt: Options, gpa: Allocator) !@This() {
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

            if (opt.flux_analysis) |flux_opt| {
                const n = flux_opt.npoint;

                hist.flux_acc = try Matrix(Complex(T)).initZero(n * nstate, wfn_nrow, gpa);

                hist.flux_phs = try Vector(Complex(T)).init(n, gpa);

                hist.wfn_init = try Matrix(Complex(T)).init(nstate, wfn_nrow, gpa);
            }

            return hist;
        }

        /// Deallocates the wavefunction and observable history matrices.
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

            if (self.flux_acc) |*flux_acc| flux_acc.deinit(gpa);
            if (self.flux_phs) |*flux_phs| flux_phs.deinit(gpa);
            if (self.wfn_init) |*wfn_init| wfn_init.deinit(gpa);
        }

        /// Appends the current wavefunction and observables to the propagation history.
        pub fn append(self: *@This(), wfn: Wavefunction(T), obs: Observables(T), opt: Options) void {
            const step_idx, const t = .{ self.index, @as(T, @floatFromInt(self.index)) * opt.time_step };

            if (self.wfn) |*storage| for (0..wfn.W.ncol()) |j| for (0..wfn.W.nrow()) |i| {
                storage.ptr(j, 2 * (step_idx * wfn.W.nrow() + i) + 0).* = wfn.W.at(i, j).re;
                storage.ptr(j, 2 * (step_idx * wfn.W.nrow() + i) + 1).* = wfn.W.at(i, j).im;
            };

            if (self.flux_acc) |*accum| {
                const flux_opt, const phases = .{ opt.flux_analysis.?, self.flux_phs.? };

                const n_energies, const n_states, const n_grid = .{ phases.length(), wfn.W.nrow(), wfn.W.ncol() };

                for (0..n_energies) |ei| {
                    const E = FluxAnalysis(T).getCollisionEnergy(flux_opt, ei) + self.e_thrs;

                    phases.data[ei] = Complex(T).init(std.math.cos(E * t), std.math.sin(E * t));
                }

                for (0..n_states) |f| {
                    const offset = f * n_energies * n_grid;

                    const slice = accum.data[offset .. offset + n_energies * n_grid];

                    var a_mat = Matrix(Complex(T)).fromSlice(n_energies, n_grid, slice);

                    geru(Complex(T), &a_mat, phases, wfn.W.row(f), Complex(T).init(1, 0));
                }
            }

            if (step_idx == 0 and self.wfn_init != null) for (0..wfn.W.nrow()) |i| for (0..wfn.W.ncol()) |j| {
                self.wfn_init.?.ptr(i, j).* = wfn.W.at(i, j);
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

        /// Writes the accumulated history and calculated spectra to output files.
        pub fn exportWrite(self: *@This(), io: std.Io, grid: Grid(T), opt: Options, pot: Potential(T), gpa: Allocator) !?Matrix(T) {
            const dt, const end = .{ opt.time_step, opt.time_step * @as(T, @floatFromInt(self.index - 1)) };

            if (opt.write.acf) |path| {
                try writeMatrixLinspace(Complex(T), io, path, self.acf.?.takeRows(self.index).asMatrix(), 0, end);
            }

            if (opt.write.position) |path| {
                try writeMatrixLinspace(T, io, path, self.pos.?.takeRows(self.index), 0, end);
            }

            if (opt.write.momentum) |path| {
                try writeMatrixLinspace(T, io, path, self.mom.?.takeRows(self.index), 0, end);
            }

            if (opt.write.population) |path| {
                try writeMatrixLinspace(T, io, path, self.pop.?.takeRows(self.index), 0, end);
            }

            if (opt.write.potential_energy) |path| {
                try writeMatrixLinspace(T, io, path, self.epot.?.takeRows(self.index), 0, end);
            }

            if (opt.write.kinetic_energy) |path| {
                try writeMatrixLinspace(T, io, path, self.ekin.?.takeRows(self.index), 0, end);
            }

            if (opt.write.norm) |path| {
                try writeMatrixLinspace(T, io, path, self.norm.?.takeRows(self.index), 0, end);
            }

            if (opt.write.total_energy) |path| {
                try writeMatrixLinspace(T, io, path, self.etot.?.takeRows(self.index), 0, end);
            }

            if (opt.write.wavefunction) |path| {
                if (grid.r) |r| {
                    try writeMatrixHjoin(T, io, path, r, null, self.wfn.?, 2 * pot.nstate() * self.index);
                }

                if (grid.r == null) {
                    var temp_r = try Matrix(T).init(grid.nrow(), grid.ncol(), gpa);
                    defer temp_r.deinit(gpa);

                    for (0..grid.nrow()) |i| for (0..grid.ncol()) |j| {
                        temp_r.ptr(i, j).* = grid.getR(i, j);
                    };

                    try writeMatrixHjoin(T, io, path, temp_r, null, self.wfn.?, 2 * pot.nstate() * self.index);
                }
            }

            if (opt.spectrum) |spec| {
                const rows = self.acf.?.takeRows(self.index);

                var acfpd, var sigma = try calcSpectrum(T, rows, dt, spec.padding, spec.threshold, gpa);

                defer {
                    acfpd.deinit(gpa);
                    sigma.deinit(gpa);
                }

                if (spec.write.acf) |path| {
                    const pad = dt * @as(T, @floatFromInt(spec.padding));

                    try writeMatrixLinspace(Complex(T), io, path, acfpd.asMatrix(), 0, end + pad);
                }

                if (spec.write.spectrum) |path| {
                    const nyquist, const nt = .{ std.math.pi / dt, @as(T, @floatFromInt(sigma.length())) };

                    try writeMatrixLinspace(T, io, path, sigma.asMatrix(), -nyquist, nyquist * (nt - 2) / nt);
                }
            }

            var result_sigma: ?Matrix(T) = null;

            if (opt.flux_analysis) |flux_opt| {
                var fa = try FluxAnalysis(T).init(opt, gpa);
                defer fa.deinit(gpa);

                var sigma = try fa.run(grid, self.wfn_init.?, self.flux_acc.?, gpa);
                errdefer sigma.deinit(gpa);

                if (flux_opt.write.cross_section) |path| {
                    const e_min: T, const e_max: T = .{ @floatCast(flux_opt.e_min), @floatCast(flux_opt.e_max) };

                    if (flux_opt.logarithmic) {
                        try writeMatrixLogspace(T, io, path, sigma, e_min, e_max);
                    } else {
                        try writeMatrixLinspace(T, io, path, sigma, e_min, e_max);
                    }
                }

                result_sigma = sigma;
            }

            return result_sigma;
        }
    };
}

/// Worker routine wrapper for concurrent partial wavepacket propagation tasks.
fn PartialWaveWorker(comptime T: type) type {
    return struct {
        /// Executes partial wave propagation worker thread.
        pub fn run(ctx: *PartialWaveContext(T), io: std.Io, opt: Options, log: bool, gpa: Allocator, thread_sigma: ?*Matrix(T)) void {
            ctx.worker(io, opt, log, gpa, thread_sigma, init, solve);
        }
    };
}

/// Bundles the grids, Hamiltonian, wavefunction, propagator, and boundaries for simulation.
fn SimulationState(comptime T: type) type {
    return struct {
        wfn_kpgrids: Grid(T),
        hams: Hamiltonian(T),
        epoten: Potential(T),
        wfn: Wavefunction(T),
        propg: Propagator(T),
        pop_apabs: Vector(T),

        orthw: std.ArrayList(Wavefunction(T)),

        /// Deallocates all resources held in the quantum simulation state.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            for (0..self.orthw.items.len) |i| self.orthw.items[i].deinit(gpa);

            inline for (@typeInfo(@This()).@"struct".field_names) |name| {
                @field(self, name).deinit(gpa);
            }
        }
    };
}

/// Bundles simulation options, state, target eigenvalue index, and logging flags.
fn SolveContext(comptime T: type) type {
    return struct { opt: Options, sim: *SimulationState(T), eigs: usize, log: bool };
}

/// Executes the grid-based split-operator wavepacket propagation simulation.
pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator) !Result(T) {
    try validation.checkInvalidInput(opt);

    if (opt.partial_waves != null) {
        return try runPartialWaves(T, io, opt, log, gpa);
    }

    var result: Result(T) = .{ .observables = .empty };
    errdefer result.deinit(gpa);

    if (log) try std.Io.File.stdout().writeStreamingAll(io, "\nQUANTUM DYNAMICS INIT: ");

    var timer = std.Io.Timestamp.now(io, .real);

    var sim = try init(T, io, opt, gpa);
    defer sim.deinit(gpa);

    if (log) try printf(io, "{f}\n", .{timer.untilNow(io, .real)});

    for (0..if (opt.imaginary) |imag| imag.nstate else 1) |i| {
        var res = try solve(T, io, .{ .opt = opt, .sim = &sim, .eigs = i, .log = log }, gpa);
        errdefer res.deinit(gpa);

        if (i < if (opt.imaginary) |imag| imag.nstate else 0) {
            var cloned = try sim.wfn.clone(gpa);
            errdefer cloned.deinit(gpa);

            try sim.orthw.append(gpa, cloned);
        }

        var obs = res.observables.swapRemove(0);
        errdefer obs.deinit(gpa);

        res.observables.deinit(gpa);

        if (res.cross_section) |cs| {
            if (result.cross_section) |*old| old.deinit(gpa);

            result.cross_section = cs;
        }

        try result.observables.append(gpa, obs);

        if (log) {
            try printFinalPop(T, io, obs);
        }
    }

    if (log) try printFinalEnergies(T, io, result.observables);

    return result;
}

/// Initializes the grid, wavefunction, Hamiltonian, and Fourier plans.
fn init(comptime T: type, io: std.Io, opt: Options, gpa: Allocator) !SimulationState(T) {
    var pot = try Potential(T).init(io, opt.potential, gpa);
    errdefer pot.deinit(gpa);

    if (opt.initial_conditions.state >= pot.nstate()) {
        std.log.err("REQUESTED ELECTRONIC STATE EXCEEDS NUMBER OF AVAILABLE STATES", .{});

        return error.InvalidInput;
    }

    const dt = if (opt.imaginary) |_| Complex(T).init(0, -opt.time_step) else Complex(T).init(opt.time_step, 0);

    var grid = try Grid(T).init(opt.grid.bounds, opt.grid.npoint, opt.grid.cylindrical, !opt.memory.grid, gpa);
    errdefer grid.deinit(gpa);

    var wfn = try Wavefunction(T).init(pot.ndim(), pot.nstate(), opt.grid.npoint, gpa);
    errdefer wfn.deinit(gpa);

    const mass = try gpa.alloc(T, opt.mass.len);
    defer gpa.free(mass);

    for (opt.mass, 0..) |m, i| {
        mass[i] = @floatCast(m);
    }

    var ham = try Hamiltonian(T).init(grid, pot, mass, opt.j_quantum_number, !opt.memory.potential, gpa);
    errdefer ham.deinit(gpa);

    var prop = try Propagator(T).init(grid, ham, pot, opt.absorbing_potential, dt, !opt.memory.propagator, gpa);
    errdefer prop.deinit(gpa);

    var pop_apabs = try Vector(T).initZero(pot.nstate(), gpa);
    errdefer pop_apabs.deinit(gpa);

    var ss: SimulationState(T) = undefined;

    ss.epoten = pot;
    ss.propg = prop;

    ss.wfn, ss.wfn_kpgrids, ss.hams, ss.pop_apabs, ss.orthw = .{ wfn, grid, ham, pop_apabs, .empty };

    return ss;
}

/// Prints final energies of target vibronic states to stdout.
fn printFinalEnergies(comptime T: type, io: std.Io, obs: std.ArrayList(Observables(T))) !void {
    try std.Io.File.stdout().writeStreamingAll(io, "\n");

    for (0..obs.items.len) |i| {
        const ekin = obs.items[i].ekin orelse std.math.nan(T);
        const epot = obs.items[i].epot orelse std.math.nan(T);

        const etot = ekin + epot;

        try printf(io, "FINAL ENERGY OF VIBRONIC STATE {d:02}: {d:13.8} Eh\n", .{ i, etot });
    }
}

/// Prints the final electronic state populations to stdout.
fn printFinalPop(comptime T: type, io: std.Io, obs: Observables(T)) !void {
    if (obs.pop) |pop| {
        try std.Io.File.stdout().writeStreamingAll(io, "\n");

        for (0..pop.length()) |i| {
            try printf(io, "FINAL POPULATION OF ELECTRONIC STATE {d:02}: {d:.8}\n", .{ i, pop.at(i) });
        }
    }
}

/// Prints column headers for real-time quantum dynamics logging.
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

/// Prints the current time step's quantum observables and wall time.
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

/// Computes total reactive cross sections by summing partial wave contributions across angular momentum states.
fn runPartialWaves(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator) !Result(T) {
    const pw = opt.partial_waves.?;
    const fa = opt.flux_analysis.?;

    const j_max = if (pw.max) |m| m else try estimateMaxJ(T, io, opt, gpa);

    const nthreads = @min(pw.threads orelse 1, j_max / pw.step + 1);

    if (nthreads < 1) {
        std.log.err("PARTIAL WAVES THREAD COUNT MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (log) {
        const params = .{ 0, j_max, pw.step, nthreads };

        try printf(io, "\nPARTIAL-WAVES RUNS FROM J={d} TO J={d} WITH STEP {d} USING {d} THREADS\n", params);

        const cols = .{ "INDEX", "J", "VCENT (Eh)", "PEAK P (-)", "INTEG P (Eh)", "TIME" };

        try printf(io, "{s:6} {s:7} {s:12} {s:12} {s:12} {s}\n", cols);
    }

    const thread_sigmas = if (pw.write.cross_section != null or opt.flux_analysis != null) try gpa.alloc(Matrix(T), nthreads) else null;

    defer if (thread_sigmas) |sigmas| {
        for (sigmas) |*s| if (s.data.len > 0) s.deinit(gpa);

        gpa.free(sigmas);
    };

    if (thread_sigmas) |sigmas| for (sigmas) |*s| {
        s.* = .{ .data = &.{}, .shape = .{ 0, 0 } };
    };

    var ctx = PartialWaveContext(T){ .j_max = j_max, .next_j = .init(0) };

    if (nthreads == 1) {
        ctx.worker(io, opt, log, gpa, if (thread_sigmas) |sigmas| &sigmas[0] else null, init, solve);
    }

    if (nthreads > 1) {
        var threads = try gpa.alloc(std.Thread, nthreads);
        defer gpa.free(threads);

        for (0..nthreads) |t| {
            const sigma = if (thread_sigmas) |sigmas| &sigmas[t] else null;

            threads[t] = try std.Thread.spawn(.{}, PartialWaveWorker(T).run, .{ &ctx, io, opt, log, gpa, sigma });
        }

        for (threads) |t| t.join();
    }

    var total: ?Matrix(T) = null;
    errdefer if (total) |*t| t.deinit(gpa);

    if (thread_sigmas) |sigmas| for (sigmas) |s| if (s.data.len > 0) {
        if (total == null) total = try Matrix(T).initZero(s.nrow(), s.ncol(), gpa);

        for (0..s.data.len) |i| total.?.data[i] += s.data[i];
    };

    if (total) |t| if (pw.write.cross_section) |cs| {
        const e_min: T, const e_max: T = .{ @floatCast(fa.e_min), @floatCast(fa.e_max) };

        if (fa.logarithmic) {
            try writeMatrixLogspace(T, io, cs, t, e_min, e_max);
        }

        if (!fa.logarithmic) {
            try writeMatrixLinspace(T, io, cs, t, e_min, e_max);
        }
    };

    return Result(T){ .observables = .empty, .cross_section = total };
}

/// Propagates the wavepacket over the specified iterations using split-operator steps.
fn solve(comptime T: type, io: std.Io, ctx: SolveContext(T), gpa: Allocator) !Result(T) {
    const opt, var e_thrs: T = .{ ctx.opt, 0 };

    if (opt.flux_analysis) |fa| {
        const asymp: T = @floatCast(fa.asymptotic_offset);

        e_thrs = (try FluxAnalysis(T).getThreshold(opt, ctx.sim.epoten, gpa)) - asymp;
    }

    const ndim, const nstate, const npoint = .{ ctx.sim.epoten.ndim(), ctx.sim.epoten.nstate(), opt.grid.npoint };

    const neig = if (opt.imaginary) |imag| imag.nstate else 1;

    if (ctx.log) try printHeader(io, ctx.eigs, ndim, nstate, neig);

    ctx.sim.wfn.setGaussian(opt.initial_conditions, ctx.sim.wfn_kpgrids);

    ctx.sim.pop_apabs.zero();

    const c_acf = opt.spectrum != null or opt.write.acf != null;

    var hist = try History(T).init(ndim, nstate, npoint, opt.iterations + 1, opt.write, c_acf, opt, gpa);
    defer hist.deinit(gpa);

    hist.e_thrs = e_thrs;

    if (opt.initial_conditions.adiabatic) {
        try ctx.sim.wfn.toDia(ctx.sim.hams, ctx.sim.wfn_kpgrids, ctx.sim.epoten, 0, gpa);
    }

    const track_cap_pop = if (opt.absorbing_potential) |ap| ap.track_population else false;

    var wfn0 = if (c_acf) try ctx.sim.wfn.clone(gpa) else null;
    defer if (wfn0) |*w| w.deinit(gpa);

    var timer = std.Io.Timestamp.now(io, .real);

    for (0..opt.iterations + 1) |i| {
        const time = (@as(T, @floatFromInt(i)) - 0.5) * opt.time_step;

        const ap, const obst = .{ opt.absorbing_potential, @as(T, @floatFromInt(i)) * opt.time_step };

        if (i > 0 and ctx.sim.epoten.isTd()) {
            try ctx.sim.hams.update(ctx.sim.wfn_kpgrids, ctx.sim.epoten, time, gpa);

            try ctx.sim.propg.update(ctx.sim.wfn_kpgrids, ctx.sim.hams, ctx.sim.epoten, time, ap);
        }

        if (i > 0) try ctx.sim.propg.step(ctx.sim, ap, opt.adiabatic, track_cap_pop, time, gpa);

        if (opt.imaginary != null) for (0..ctx.sim.orthw.items.len) |j| {
            const overlap = ctx.sim.orthw.items[j].overlap(ctx.sim.wfn, ctx.sim.wfn_kpgrids);

            for (0..ctx.sim.wfn.W.data.len) |k| {
                ctx.sim.wfn.W.data[k] = ctx.sim.wfn.W.data[k].sub(overlap.mul(ctx.sim.orthw.items[j].W.data[k]));
            }
        };

        if (opt.imaginary != null) ctx.sim.wfn.normalize(ctx.sim.wfn_kpgrids);

        const is_log_step = ctx.log and ((i % opt.log_interval == 0) or (i == opt.iterations));

        if (ctx.sim.epoten.isTd()) {
            try ctx.sim.hams.update(ctx.sim.wfn_kpgrids, ctx.sim.epoten, obst, gpa);
        }

        var obs = try Observables(T).init(ctx.sim, wfn0, opt.write, opt.adiabatic, is_log_step, obst, gpa);
        defer obs.deinit(gpa);

        const need_a_wfn = (opt.write.wavefunction != null or opt.flux_analysis != null) and opt.adiabatic;

        if (need_a_wfn) {
            try ctx.sim.wfn.toAdia(ctx.sim.hams, ctx.sim.wfn_kpgrids, ctx.sim.epoten, obst, gpa);
        }

        hist.append(ctx.sim.wfn, obs, opt);

        if (need_a_wfn) {
            try ctx.sim.wfn.toDia(ctx.sim.hams, ctx.sim.wfn_kpgrids, ctx.sim.epoten, obst, gpa);
        }

        if (is_log_step) {
            try printIteration(T, io, obs, i, &timer);
        }

        if (opt.absorbing_potential) |cap| if (cap.stop_norm) |stop_norm| {
            const norm = ctx.sim.wfn.norm(ctx.sim.wfn_kpgrids);

            if (norm < stop_norm) {
                var stop_o = try Observables(T).init(ctx.sim, wfn0, opt.write, opt.adiabatic, true, obst, gpa);
                defer stop_o.deinit(gpa);

                if (ctx.log and !is_log_step) try printIteration(T, io, stop_o, i, &timer);

                break;
            }
        };
    }

    const maybe_sigma = try hist.exportWrite(io, ctx.sim.wfn_kpgrids, opt, ctx.sim.epoten, gpa);

    var result: Result(T) = .{ .observables = .empty, .cross_section = maybe_sigma };
    errdefer result.deinit(gpa);

    const time = @as(T, @floatFromInt(opt.iterations)) * opt.time_step;

    var obs = try Observables(T).init(ctx.sim, wfn0, opt.write, opt.adiabatic, true, time, gpa);
    errdefer obs.deinit(gpa);

    try result.observables.append(gpa, obs);

    return result;
}
