//! Classical molecular dynamics trajectory simulation on adiabatic or diabatic potential energy surfaces.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;

const classical_forces = @import("classical_forces.zig");
const ensemble_mod = @import("ensemble.zig");
const initial_conditions_mod = @import("initial_conditions.zig");

const Ehrenfest = @import("ehrenfest_dynamics.zig").Ehrenfest;
const EhrenfestOptions = @import("ehrenfest_dynamics.zig").Options;
const Matrix = @import("tensor.zig").Matrix;
const MethodOptions = @import("potential.zig").MethodOptions;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Potential = @import("potential.zig").Potential;
const PotentialOptions = @import("potential.zig").Options;
const SurfaceHopping = @import("surface_hopping.zig").SurfaceHopping;
const SurfaceHoppingOptions = @import("surface_hopping.zig").Options;
const Thermostat = @import("thermostat.zig").Thermostat;
const ThermostatOptions = @import("thermostat.zig").Options;
const Vector = @import("tensor.zig").Vector;
const validation = @import("validation.zig");

const getMass = @import("constants.zig").getMass;
const getTrajectoryPath = @import("read_write.zig").getTrajectoryPath;
const printf = @import("read_write.zig").printf;
const writeMatrixLinspace = @import("read_write.zig").writeMatrixLinspace;
const writeXyzFrame = @import("read_write.zig").writeXyzFrame;
const writeXyzMomentumFrame = @import("read_write.zig").writeXyzMomentumFrame;
const writeXyzMomentumTrajectory = @import("read_write.zig").writeXyzMomentumTrajectory;
const writeXyzTrajectory = @import("read_write.zig").writeXyzTrajectory;

const AMU2AU = @import("constants.zig").AMU2AU;
const AU2K = @import("constants.zig").AU2K;

pub const DistributionOptions = initial_conditions_mod.DistributionOptions;
pub const Ensemble = ensemble_mod.Ensemble;
pub const GradientBuffer = classical_forces.GradientBuffer;
pub const InitialConditions = initial_conditions_mod.InitialConditions;
pub const MaxwellBoltzmannOptions = initial_conditions_mod.MaxwellBoltzmannOptions;
pub const WignerOptions = initial_conditions_mod.WignerOptions;

const sampleWignerAbInitio = initial_conditions_mod.sampleWignerAbInitio;

/// Configuration options for the classical molecular dynamics simulation.
pub const Options = struct {
    adiabatic: bool = true,
    initial_conditions: InitialConditions = .{},
    iterations: u32,
    log_interval: u32 = 1,
    mass: ?[]const f64 = null,
    nonadiabatic: ?NonadiabaticOptions = null,
    potential: PotentialOptions,
    stream: Stream = .{},
    thermostat: ?ThermostatOptions = null,
    time_step: f64,
    trajectories: u32,
    write: Write = .{},
};

/// Tagged union for multi-state non-adiabatic trajectory propagation methods.
pub const NonadiabaticOptions = union(enum) {
    ehrenfest: EhrenfestOptions,
    surface_hopping: SurfaceHoppingOptions,
};

/// Output file destinations for streaming trajectory observables to disk during dynamics.
const Stream = struct {
    kinetic_energy: ?[]const u8 = null,
    momentum: ?[]const u8 = null,
    population: ?[]const u8 = null,
    position: ?[]const u8 = null,
    potential_energy: ?[]const u8 = null,
    state_coefficient: ?[]const u8 = null,
    state_potential_energy: ?[]const u8 = null,
    temperature: ?[]const u8 = null,
    total_energy: ?[]const u8 = null,
};

/// Output paths for recording trajectory observables to disk during dynamics.
const Write = struct {
    kinetic_energy: ?[]const u8 = null,
    momentum: ?[]const u8 = null,
    population: ?[]const u8 = null,
    position: ?[]const u8 = null,
    potential_energy: ?[]const u8 = null,
    state_coefficient: ?[]const u8 = null,
    state_potential_energy: ?[]const u8 = null,
    temperature: ?[]const u8 = null,
    total_energy: ?[]const u8 = null,
};

/// Container structure holding the results and observables computed during classical propagation.
pub fn Result(comptime T: type) type {
    return struct {
        observables: Observables(T),

        /// Deallocates memory associated with the trajectory simulation result observables.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.observables.deinit(gpa);
        }
    };
}

/// Manages real-time streaming output files for multi-trajectory observable logging during dynamics.
const StreamWriter = struct {
    files_ekin: ?[]std.Io.File = null,
    files_epot: ?[]std.Io.File = null,
    files_etot: ?[]std.Io.File = null,
    files_temp: ?[]std.Io.File = null,

    files_mom: ?[]std.Io.File = null,
    files_pop: ?[]std.Io.File = null,
    files_pos: ?[]std.Io.File = null,

    files_state_coef: ?[]std.Io.File = null,
    files_state_epot: ?[]std.Io.File = null,

    /// Initializes and opens trajectory output files on disk for configured properties.
    pub fn init(io: std.Io, stream: Stream, atoms: []const i32, ndim: usize, ntraj: usize, nstep: usize, nstate: usize, gpa: Allocator) !?@This() {
        const any_stream = inline for (@typeInfo(Stream).@"struct".field_names) |name| {
            if (@field(stream, name) != null) break true;
        } else false;

        if (!any_stream) return null;

        var self = @This(){};
        errdefer self.deinit(io, gpa);

        if (stream.kinetic_energy) |path| {
            self.files_ekin = try initMatrixFiles(io, path, ntraj, nstep, 2, gpa);
        }

        if (stream.momentum) |path| {
            self.files_mom = if (atoms.len > 0)
                try initFiles(io, path, ntraj, gpa)
            else
                try initMatrixFiles(io, path, ntraj, nstep, ndim + 1, gpa);
        }

        if (stream.population) |path| {
            self.files_pop = try initMatrixFiles(io, path, ntraj, nstep, nstate + 1, gpa);
        }

        if (stream.position) |path| {
            self.files_pos = if (atoms.len > 0)
                try initFiles(io, path, ntraj, gpa)
            else
                try initMatrixFiles(io, path, ntraj, nstep, ndim + 1, gpa);
        }

        if (stream.potential_energy) |path| {
            self.files_epot = try initMatrixFiles(io, path, ntraj, nstep, 2, gpa);
        }

        if (stream.state_coefficient) |path| {
            self.files_state_coef = try initMatrixFiles(io, path, ntraj, nstep, nstate + 1, gpa);
        }

        if (stream.state_potential_energy) |path| {
            self.files_state_epot = try initMatrixFiles(io, path, ntraj, nstep, nstate + 1, gpa);
        }

        if (stream.temperature) |path| {
            self.files_temp = try initMatrixFiles(io, path, ntraj, nstep, 2, gpa);
        }

        if (stream.total_energy) |path| {
            self.files_etot = try initMatrixFiles(io, path, ntraj, nstep, 2, gpa);
        }

        return self;
    }

    /// Closes all trajectory output files and deallocates handles.
    pub fn deinit(self: *@This(), io: std.Io, gpa: Allocator) void {
        inline for (@typeInfo(@This()).@"struct".field_names) |name| {
            if (@field(self, name)) |files| {
                for (files) |*file| file.close(io);

                gpa.free(files);
            }
        }
    }

    /// Appends the current step's energies and temperature to trajectory files.
    pub fn writeEnergy(self: @This(), io: std.Io, comptime T: type, ensemble: Ensemble(T), gb: GradientBuffer(T), nstate: usize, coefs: ?Matrix(Complex(T)), time: T) !void {
        const no_e = self.files_ekin == null and self.files_epot == null;
        if (no_e and self.files_etot == null and self.files_temp == null) return;

        const ntraj = ensemble.r.nrow();

        var buffer: [65536]u8 = undefined;

        for (0..ntraj) |k| {
            var ekin: T = 0;

            if (self.files_ekin != null or self.files_etot != null or self.files_temp != null) {
                for (0..ensemble.p.ncol()) |j| {
                    ekin += ensemble.p.at(k, j) * ensemble.p.at(k, j) / (2 * ensemble.m[j]);
                }
            }

            var epot: T = 0;

            if (self.files_epot != null or self.files_etot != null) {
                if (coefs) |coef| {
                    if (gb.adia) {
                        for (0..nstate) |l| {
                            epot += coef.at(k, l).squaredMagnitude() * gb.W.at(k, l);
                        }
                    }

                    if (!gb.adia) {
                        for (0..nstate) |l| for (0..nstate) |m| {
                            const rho_lm = coef.at(k, l).conjugate().mul(coef.at(k, m)).re;

                            epot += rho_lm * gb.V.at(k, l * nstate + m);
                        };
                    }
                }

                if (coefs == null and gb.adia) {
                    epot = gb.W.at(k, ensemble.s.at(k));
                }

                if (coefs == null and !gb.adia) {
                    const sk = ensemble.s.at(k);

                    epot = gb.V.at(k, sk * nstate + sk);
                }
            }

            if (self.files_ekin) |files| {
                var writer = files[k].writerStreaming(io, &buffer);

                try writer.interface.print("{d:20.14} {d:20.14}\n", .{ time, ekin });

                try writer.interface.flush();
            }

            if (self.files_epot) |files| {
                var writer = files[k].writerStreaming(io, &buffer);

                try writer.interface.print("{d:20.14} {d:20.14}\n", .{ time, epot });

                try writer.interface.flush();
            }

            if (self.files_temp) |files| {
                const ndof_t = @as(T, @floatFromInt(ensemble.ndof));

                const temp_val = if (ensemble.ndof == 0) 0 else (2 * ekin / ndof_t) * AU2K;

                var writer = files[k].writerStreaming(io, &buffer);

                try writer.interface.print("{d:20.14} {d:20.14}\n", .{ time, temp_val });

                try writer.interface.flush();
            }

            if (self.files_etot) |files| {
                var writer = files[k].writerStreaming(io, &buffer);

                try writer.interface.print("{d:20.14} {d:20.14}\n", .{ time, ekin + epot });

                try writer.interface.flush();
            }
        }
    }

    /// Appends the current step's coordinates to trajectory files in XYZ or space-separated format.
    pub fn writeFrame(self: @This(), io: std.Io, comptime T: type, atoms: []const i32, r: Matrix(T), s: ?Vector(usize), time: T) !void {
        if (self.files_pos) |files| {
            var buffer: [65536]u8 = undefined;

            for (0..files.len) |k| {
                var writer = files[k].writerStreaming(io, &buffer);

                if (atoms.len > 0) {
                    const state = if (s) |st| st.at(k) else null;

                    try writeXyzFrame(T, &writer, atoms, r.rowSlice(k), time, state);
                }

                if (atoms.len == 0) {
                    try writer.interface.print("{d:20.14}", .{time});

                    for (0..r.ncol()) |j| {
                        try writer.interface.print(" {d:20.14}", .{r.at(k, j)});
                    }

                    try writer.interface.print("\n", .{});
                }

                try writer.interface.flush();
            }
        }
    }

    /// Appends the current step's momentum vectors to trajectory files in XYZ or space-separated format.
    pub fn writeMomentumFrame(self: @This(), io: std.Io, comptime T: type, atoms: []const i32, p: Matrix(T), s: ?Vector(usize), time: T) !void {
        if (self.files_mom) |files| {
            var buffer: [65536]u8 = undefined;

            for (0..files.len) |k| {
                var writer = files[k].writerStreaming(io, &buffer);

                if (atoms.len > 0) {
                    const state = if (s) |st| st.at(k) else null;

                    try writeXyzMomentumFrame(T, &writer, atoms, p.rowSlice(k), time, state);
                }

                if (atoms.len == 0) {
                    try writer.interface.print("{d:20.14}", .{time});

                    for (0..p.ncol()) |j| {
                        try writer.interface.print(" {d:20.14}", .{p.at(k, j)});
                    }

                    try writer.interface.print("\n", .{});
                }

                try writer.interface.flush();
            }
        }
    }

    /// Appends the current step's electronic state populations to trajectory files.
    pub fn writePopulation(self: @This(), io: std.Io, comptime T: type, ensemble: Ensemble(T), nstate: usize, coefs: ?Matrix(Complex(T)), time: T) !void {
        if (self.files_pop) |files| {
            var buffer: [65536]u8 = undefined;

            for (0..files.len) |k| {
                var writer = files[k].writerStreaming(io, &buffer);

                try writer.interface.print("{d:20.14}", .{time});

                for (0..nstate) |m| {
                    var pop_m: T = 0;

                    if (coefs) |coef| {
                        pop_m = coef.at(k, m).squaredMagnitude();
                    }

                    if (coefs == null) {
                        pop_m = if (m == ensemble.s.at(k)) 1 else 0;
                    }

                    try writer.interface.print(" {d:20.14}", .{pop_m});
                }

                try writer.interface.print("\n", .{});

                try writer.interface.flush();
            }
        }
    }

    /// Appends the current step's electronic state coefficient magnitudes to trajectory files.
    pub fn writeStateCoefficient(self: @This(), io: std.Io, comptime T: type, nstate: usize, coefs: ?Matrix(Complex(T)), time: T) !void {
        if (self.files_state_coef) |files| {
            var buffer: [65536]u8 = undefined;

            for (0..files.len) |k| {
                var writer = files[k].writerStreaming(io, &buffer);

                try writer.interface.print("{d:20.14}", .{time});

                for (0..nstate) |m| {
                    const c2 = if (coefs) |coef| coef.at(k, m).squaredMagnitude() else 0;

                    try writer.interface.print(" {d:20.14}", .{c2});
                }

                try writer.interface.print("\n", .{});

                try writer.interface.flush();
            }
        }
    }

    /// Appends the current step's electronic state potential energies to trajectory files.
    pub fn writeStatePotentialEnergy(self: @This(), io: std.Io, comptime T: type, gb: GradientBuffer(T), nstate: usize, time: T) !void {
        if (self.files_state_epot) |files| {
            var buffer: [65536]u8 = undefined;

            for (0..files.len) |k| {
                var writer = files[k].writerStreaming(io, &buffer);

                try writer.interface.print("{d:20.14}", .{time});

                for (0..nstate) |s| {
                    const epot = if (gb.adia) gb.W.at(k, s) else gb.V.at(k, s * nstate + s);

                    try writer.interface.print(" {d:20.14}", .{epot});
                }

                try writer.interface.print("\n", .{});

                try writer.interface.flush();
            }
        }
    }

    /// Opens output files for each trajectory in the ensemble.
    fn initFiles(io: std.Io, path: []const u8, ntraj: usize, gpa: Allocator) ![]std.Io.File {
        const files = try gpa.alloc(std.Io.File, ntraj);
        errdefer gpa.free(files);

        var count: usize = 0;

        errdefer for (0..count) |k| files[k].close(io);

        while (count < ntraj) : (count += 1) {
            const fname = try getTrajectoryPath(path, count, ntraj, gpa);
            defer if (ntraj > 1) gpa.free(fname);

            files[count] = try std.Io.Dir.cwd().createFile(io, fname, .{});
        }

        return files;
    }

    /// Opens trajectory files and writes the matrix dimensions header to each stream.
    fn initMatrixFiles(io: std.Io, path: []const u8, ntraj: usize, nstep: usize, ncol: usize, gpa: Allocator) ![]std.Io.File {
        const files = try initFiles(io, path, ntraj, gpa);

        errdefer {
            for (files) |*file| file.close(io);

            gpa.free(files);
        }

        var buffer: [1024]u8 = undefined;

        for (files) |file| {
            var writer = file.writerStreaming(io, &buffer);

            try writer.interface.print("{d} {d}\n", .{ nstep, ncol });

            try writer.interface.flush();
        }

        return files;
    }
};

/// Accumulates time-dependent expectation values and dynamics trajectories.
fn History(comptime T: type) type {
    return struct {
        pos: ?Matrix(T) = null,
        mom: ?Matrix(T) = null,
        pop: ?Matrix(T) = null,

        state_coef: ?Matrix(T) = null,
        state_epot: ?Matrix(T) = null,

        epot: ?Matrix(T) = null,
        ekin: ?Matrix(T) = null,
        temp: ?Matrix(T) = null,
        etot: ?Matrix(T) = null,

        index: usize = 0,

        /// Allocates memory for storing dynamics history of positions, momenta, and populations.
        pub fn init(ndim: usize, nstate: usize, iters: usize, write: Write, gpa: Allocator) !@This() {
            var hist = @This(){};
            errdefer hist.deinit(gpa);

            var store_ekin, var store_epot = .{ write.kinetic_energy != null, write.potential_energy != null };

            store_epot = store_epot or write.total_energy != null;
            store_ekin = store_ekin or write.total_energy != null;

            const store_etot = write.total_energy != null;

            if (write.position != null) hist.pos = try Matrix(T).init(iters, ndim, gpa);
            if (write.momentum != null) hist.mom = try Matrix(T).init(iters, ndim, gpa);

            if (write.population != null) {
                hist.pop = try Matrix(T).init(iters, nstate, gpa);
            }

            if (write.state_coefficient != null) {
                hist.state_coef = try Matrix(T).init(iters, nstate, gpa);
            }

            if (write.state_potential_energy != null) {
                hist.state_epot = try Matrix(T).init(iters, nstate, gpa);
            }

            if (store_ekin) hist.ekin = try Matrix(T).init(iters, 1, gpa);
            if (store_epot) hist.epot = try Matrix(T).init(iters, 1, gpa);
            if (store_etot) hist.etot = try Matrix(T).init(iters, 1, gpa);

            if (write.temperature != null) hist.temp = try Matrix(T).init(iters, 1, gpa);

            return hist;
        }

        /// Deallocates history storage arrays.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            if (self.pos) |*pos| pos.deinit(gpa);
            if (self.mom) |*mom| mom.deinit(gpa);
            if (self.pop) |*pop| pop.deinit(gpa);

            if (self.state_coef) |*state_coef| state_coef.deinit(gpa);
            if (self.state_epot) |*state_epot| state_epot.deinit(gpa);

            if (self.epot) |*epot| epot.deinit(gpa);
            if (self.ekin) |*ekin| ekin.deinit(gpa);
            if (self.temp) |*temp| temp.deinit(gpa);
            if (self.etot) |*etot| etot.deinit(gpa);
        }

        /// Appends the current step's observables to the simulation history.
        pub fn append(self: *@This(), obs: Observables(T)) void {
            const step_idx = self.index;

            if (self.pos) |*pos| if (obs.pos) |v| {
                for (0..v.length()) |j| pos.ptr(step_idx, j).* = v.at(j);
            };

            if (self.mom) |*mom| if (obs.mom) |v| {
                for (0..v.length()) |j| mom.ptr(step_idx, j).* = v.at(j);
            };

            if (self.pop) |*pop| if (obs.pop) |v| {
                for (0..v.length()) |j| pop.ptr(step_idx, j).* = v.at(j);
            };

            if (self.state_coef) |*state_coef| if (obs.state_coef) |v| {
                for (0..v.length()) |j| state_coef.ptr(step_idx, j).* = v.at(j);
            };

            if (self.state_epot) |*state_epot| if (obs.state_epot) |v| {
                for (0..v.length()) |j| state_epot.ptr(step_idx, j).* = v.at(j);
            };

            if (self.epot) |*epot| {
                epot.ptr(step_idx, 0).* = obs.epot.?;
            }

            if (self.ekin) |*ekin| {
                ekin.ptr(step_idx, 0).* = obs.ekin.?;
            }

            if (self.temp) |*temp| {
                temp.ptr(step_idx, 0).* = obs.temp.?;
            }

            if (self.etot) |*etot| {
                etot.ptr(step_idx, 0).* = obs.ekin.? + obs.epot.?;
            }

            self.index += 1;
        }

        /// Writes the accumulated history of dynamics observables to output files.
        pub fn exportWrite(self: *@This(), io: std.Io, atoms: []const i32, dt: f64, write: Write) !void {
            const end = dt * @as(T, @floatFromInt(self.index - 1));

            if (write.position) |path| {
                if (atoms.len > 0) {
                    try writeXyzTrajectory(T, io, path, atoms, self.pos.?.takeRows(self.index), @floatCast(dt));
                }

                if (atoms.len == 0) {
                    try writeMatrixLinspace(T, io, path, self.pos.?.takeRows(self.index), 0, end);
                }
            }

            if (write.momentum) |path| {
                if (atoms.len > 0) {
                    try writeXyzMomentumTrajectory(T, io, path, atoms, self.mom.?.takeRows(self.index), @floatCast(dt));
                }

                if (atoms.len == 0) {
                    try writeMatrixLinspace(T, io, path, self.mom.?.takeRows(self.index), 0, end);
                }
            }

            if (write.population) |path| {
                try writeMatrixLinspace(T, io, path, self.pop.?.takeRows(self.index), 0, end);
            }

            if (write.potential_energy) |path| {
                try writeMatrixLinspace(T, io, path, self.epot.?.takeRows(self.index), 0, end);
            }

            if (write.state_coefficient) |path| {
                try writeMatrixLinspace(T, io, path, self.state_coef.?.takeRows(self.index), 0, end);
            }

            if (write.state_potential_energy) |path| {
                try writeMatrixLinspace(T, io, path, self.state_epot.?.takeRows(self.index), 0, end);
            }

            if (write.kinetic_energy) |path| {
                try writeMatrixLinspace(T, io, path, self.ekin.?.takeRows(self.index), 0, end);
            }

            if (write.temperature) |path| {
                try writeMatrixLinspace(T, io, path, self.temp.?.takeRows(self.index), 0, end);
            }

            if (write.total_energy) |path| {
                try writeMatrixLinspace(T, io, path, self.etot.?.takeRows(self.index), 0, end);
            }
        }
    };
}

/// Holds physical observables evaluated at a specific time step of dynamics.
fn Observables(comptime T: type) type {
    return struct {
        pos: ?Vector(T) = null,
        mom: ?Vector(T) = null,
        pop: ?Vector(T) = null,

        state_coef: ?Vector(T) = null,
        state_epot: ?Vector(T) = null,

        epot: ?T = null,
        ekin: ?T = null,
        temp: ?T = null,

        /// Computes the physical observables from the current simulation state.
        pub fn init(sim: SimulationState(T), write: Write, log: bool, has_thermo: bool, gpa: Allocator) !@This() {
            var obs = @This(){};
            errdefer obs.deinit(gpa);

            const calc_ekin, const calc_epot = .{ write.kinetic_energy != null, write.potential_energy != null };

            var calc = .{
                .pos = log or write.position != null,
                .mom = log or write.momentum != null,

                .pop = log or write.population != null,

                .ekin = log or calc_ekin,
                .epot = log or calc_epot,

                .temp = (log and has_thermo) or write.temperature != null,
            };

            calc.ekin = calc.ekin or write.total_energy != null;
            calc.epot = calc.epot or write.total_energy != null;

            if (calc.mom) obs.mom = try sim.ensemble.mom(gpa);
            if (calc.pos) obs.pos = try sim.ensemble.pos(gpa);

            const nonadia = sim.propag.nonadia_dynamic;

            const coefs = if (nonadia) |n| (if (n == .ehrenfest) n.ehrenfest.coefics else null) else null;

            if (calc.pop) obs.pop = try sim.ensemble.pop(sim.elpoten.nstate(), coefs, gpa);

            if (write.state_coefficient != null or (log and sim.propag.coefs() != null)) {
                if (sim.propag.coefs()) |c_coefs| {
                    obs.state_coef = try sim.ensemble.stateCoef(sim.elpoten.nstate(), c_coefs, gpa);
                }
            }

            if (write.state_potential_energy != null) {
                obs.state_epot = try sim.gb.epotState(sim.elpoten.nstate(), gpa);
            }

            if (calc.ekin) {
                obs.ekin = sim.ensemble.ekin();
            }

            if (calc.epot) {
                obs.epot = sim.gb.epot(sim.ensemble, sim.elpoten.nstate(), coefs);
            }

            if (calc.temp) {
                obs.temp = sim.ensemble.temp();
            }

            return obs;
        }

        /// Deallocates arrays stored in the observables struct.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            if (self.pos) |*pos| pos.deinit(gpa);
            if (self.mom) |*mom| mom.deinit(gpa);
            if (self.pop) |*pop| pop.deinit(gpa);

            if (self.state_coef) |*state_coef| state_coef.deinit(gpa);
            if (self.state_epot) |*state_epot| state_epot.deinit(gpa);
        }
    };
}

/// Numerical integrator for classical trajectories and surface hopping updates.
fn Propagator(comptime T: type) type {
    return struct {
        /// Trajectory propagation method utilizing adiabatic/diabatic surface hopping or mean-field Ehrenfest dynamics.
        pub const Namd = union(enum) {
            surface_hopping: SurfaceHopping(T),
            ehrenfest: Ehrenfest(T),
        };

        nonadia_dynamic: ?Namd = null,
        thermo: ?Thermostat(T) = null,

        dt: T,

        /// Initializes the time propagator and optional surface hopping solver.
        pub fn init(opt: Options, nstate: usize, gpa: Allocator) !@This() {
            const istate = opt.initial_conditions.state;

            var nonadia_dynamic: ?Namd = null;

            if (opt.nonadiabatic) |naopt| if (naopt == .surface_hopping) {
                const sh_opt, const trajs = .{ naopt.surface_hopping, opt.trajectories };

                const sh = try SurfaceHopping(T).init(sh_opt, nstate, trajs, istate, opt.adiabatic, gpa);

                nonadia_dynamic = .{ .surface_hopping = sh };
            };

            if (opt.nonadiabatic) |naopt| if (naopt == .ehrenfest) {
                const eh = try Ehrenfest(T).init(naopt.ehrenfest, nstate, opt.trajectories, istate, opt.adiabatic, gpa);

                nonadia_dynamic = .{ .ehrenfest = eh };
            };

            const thermo = if (opt.thermostat) |topt| Thermostat(T).init(topt, @floatCast(opt.time_step)) else null;

            return .{ .dt = @floatCast(opt.time_step), .nonadia_dynamic = nonadia_dynamic, .thermo = thermo };
        }

        /// Deallocates surface hopping resources.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            if (self.nonadia_dynamic) |*namd| switch (namd.*) {
                inline else => |*n| n.deinit(gpa),
            };
        }

        /// Returns electronic wavefunction coefficients for nonadiabatic dynamics if active.
        pub fn coefs(self: @This()) ?Matrix(Complex(T)) {
            if (self.nonadia_dynamic) |n| switch (n) {
                .surface_hopping => |sh| return sh.coefs(),
                .ehrenfest => |eh| return eh.coefs(),
            };

            return null;
        }

        /// Integrates equations of motion using the Velocity Verlet algorithm.
        pub fn step(self: *@This(), ens: *Ensemble(T), gb: *GradientBuffer(T), pot: Potential(T), time: T) !void {
            if (self.thermo) |*thermo| {
                thermo.apply(&ens.p, ens.m, ens.ndof, 0.5);
            }

            for (0..ens.r.nrow()) |i| for (0..ens.r.ncol()) |j| {
                ens.p.ptr(i, j).* += 0.5 * ens.m[j] * ens.a.at(i, j) * self.dt;

                ens.r.ptr(i, j).* += (ens.p.at(i, j) / ens.m[j]) * self.dt;
            };

            const is_ehrenfest = if (self.nonadia_dynamic) |nd| nd == .ehrenfest else false;

            try gb.update(ens.r, pot, time, if (is_ehrenfest) null else ens.s.data);

            if (self.nonadia_dynamic) |*n| if (n.* == .ehrenfest) {
                try n.ehrenfest.step(ens, gb.V, gb.W, gb.U, gb.grad_V, pot.isTd(), self.dt);
            };

            const nonadia = self.nonadia_dynamic;

            const coefics = if (nonadia) |*n| (if (n.* == .ehrenfest) &n.ehrenfest.coefics else null) else null;

            gb.apply(ens, pot, coefics);

            for (0..ens.r.nrow()) |i| for (0..ens.r.ncol()) |j| {
                ens.p.ptr(i, j).* += 0.5 * ens.m[j] * ens.a.at(i, j) * self.dt;
            };

            if (self.thermo) |*thermo| {
                thermo.apply(&ens.p, ens.m, ens.ndof, 0.5);
            }

            if (self.nonadia_dynamic) |*n| if (n.* == .surface_hopping) {
                if (try n.surface_hopping.hop(ens, gb.V, gb.W, gb.U, gb.grad_V, pot.isTd(), self.dt)) {
                    for (0..ens.s.length()) |i| if (ens.s.at(i) != n.surface_hopping.targets[i]) {
                        switch (pot) {
                            .ab_initio => |ab| {
                                try ab.updateTrajectoryGradient(&gb.V, &gb.grad_V, ens.r, i, ens.s.at(i));

                                for (0..pot.nstate()) |k| {
                                    gb.W.ptr(i, k).* = gb.V.at(i, k * pot.nstate() + k);
                                }
                            },

                            inline else => {},
                        }
                    };

                    gb.apply(ens, pot, null);
                }
            };
        }
    };
}

/// Bundles the active ensemble, potential, propagator, and gradient buffers.
fn SimulationState(comptime T: type) type {
    return struct {
        ensemble: Ensemble(T),
        gb: GradientBuffer(T),
        elpoten: Potential(T),
        propag: Propagator(T),

        /// Deallocates all resources held within the simulation state.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            inline for (@typeInfo(@This()).@"struct".field_names) |name| {
                @field(self, name).deinit(gpa);
            }
        }
    };
}

/// Bundles the simulation options, state, and logging flag for the solver.
fn SolveContext(comptime T: type) type {
    return struct { opt: Options, sim: *SimulationState(T), log: bool };
}

/// Runs the classical dynamics simulation and returns the accumulated observables.
pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator) !Result(T) {
    try validation.checkInvalidInput(opt);

    var timer = std.Io.Timestamp.now(io, .real);

    var sim = try init(T, io, opt, gpa);
    defer sim.deinit(gpa);

    if (log) try printf(io, "\nCLASSICAL DYNAMICS INIT: {f}\n", .{timer.untilNow(io, .real)});

    var obs = try solve(T, io, .{ .opt = opt, .sim = &sim, .log = log }, gpa, gpa);
    errdefer obs.deinit(gpa);

    if (log) {
        try printFinalPop(T, io, obs, opt);
    }

    return .{ .observables = obs };
}

/// Initializes the simulation state, potential, and initial ensemble.
fn init(comptime T: type, io: std.Io, opt: Options, gpa: Allocator) !SimulationState(T) {
    var pot = try Potential(T).init(io, opt.potential, gpa);
    errdefer pot.deinit(gpa);

    if (opt.initial_conditions.state >= pot.nstate()) {
        std.log.err("REQUESTED ELECTRONIC STATE EXCEEDS NUMBER OF AVAILABLE STATES", .{});

        return error.InvalidInput;
    }

    const mass = try gpa.alloc(T, pot.ndim());
    defer gpa.free(mass);

    switch (pot) {
        .ab_initio => |ab| for (0..ab.msys.atoms.len) |i| {
            const m_au = try getMass(T, ab.msys.atoms[i]) * AMU2AU;

            mass[i * 3 + 0] = m_au;
            mass[i * 3 + 1] = m_au;
            mass[i * 3 + 2] = m_au;
        },

        inline else => for (opt.mass.?, 0..) |m, i| {
            mass[i] = @floatCast(m);
        },
    }

    var ensemble = try Ensemble(T).init(pot.ndim(), opt.trajectories, mass, gpa);
    errdefer ensemble.deinit(gpa);

    var gb = try GradientBuffer(T).init(pot.ndim(), pot.nstate(), opt.trajectories, opt.adiabatic, gpa);
    errdefer gb.deinit(gpa);

    var prop = try Propagator(T).init(opt, pot.nstate(), gpa);
    errdefer prop.deinit(gpa);

    const default_pos = switch (pot) {
        .ab_initio => |ab| ab.msys.coors,

        inline else => null,
    };

    const is_ab_initio_wigner = if (opt.initial_conditions.distribution) |dist| switch (dist) {
        .wigner => |wig| wig.method != null or opt.potential == .ab_initio,
        inline else => false,
    } else false;

    if (is_ab_initio_wigner) {
        try sampleWignerAbInitio(T, io, &ensemble, opt, gpa);
    }

    if (!is_ab_initio_wigner) {
        ensemble.setGaussian(opt.initial_conditions, default_pos);
    }

    if (pot == .ab_initio) {
        ensemble.removeTranslation();

        _ = ensemble.removeRotation();

        const temp = if (opt.initial_conditions.distribution) |dist| switch (dist) {
            .maxwell_boltzmann => |mb| mb.temperature,
            inline else => null,
        } else null;

        if (temp) |t_val| {
            ensemble.rescaleTemperature(@floatCast(t_val));
        }
    }

    const is_ehrenfest = if (opt.nonadiabatic) |na| na == .ehrenfest else false;

    try gb.update(ensemble.r, pot, 0, if (is_ehrenfest) null else ensemble.s.data);

    if (prop.nonadia_dynamic) |*n| if (n.* == .surface_hopping) {
        n.surface_hopping.update(if (opt.adiabatic) gb.W else gb.V, gb.U);
    };

    if (prop.nonadia_dynamic) |*n| if (n.* == .ehrenfest) {
        n.ehrenfest.update(if (opt.adiabatic) gb.W else gb.V, gb.U);
    };

    const coefs = if (prop.nonadia_dynamic) |*n| (if (n.* == .ehrenfest) &n.ehrenfest.coefics else null) else null;

    gb.apply(&ensemble, pot, coefs);

    return .{ .ensemble = ensemble, .elpoten = pot, .propag = prop, .gb = gb };
}

/// Prints the final population fractions of electronic states to stdout.
fn printFinalPop(comptime T: type, io: std.Io, obs: Observables(T), opt: Options) !void {
    if (opt.nonadiabatic == null) return;

    if (obs.pop) |pop| {
        if (pop.length() <= 1) return;

        try std.Io.File.stdout().writeStreamingAll(io, "\n");

        const na = opt.nonadiabatic.?;

        const is_fssh = na == .surface_hopping and na.surface_hopping == .fewest_switches;

        for (0..pop.length()) |i| {
            if (is_fssh and obs.state_coef != null) {
                const fmt = "FINAL POPULATION OF ELECTRONIC STATE {d:02}: {d:.8} (COEF: {d:.8})\n";

                try printf(io, fmt, .{ i, pop.at(i), obs.state_coef.?.at(i) });
            }

            if (!is_fssh or obs.state_coef == null) {
                try printf(io, "FINAL POPULATION OF ELECTRONIC STATE {d:02}: {d:.8}\n", .{ i, pop.at(i) });
            }
        }
    }
}

/// Prints the column headers for the real-time dynamics logging output.
fn printHeader(io: std.Io, ndim: usize, nstate: usize, has_thermo: bool) !void {
    try std.Io.File.stdout().writeStreamingAll(io, "\nREAL-TIME PROPAGATION");

    const col_width = @as(usize, 12) * @min(ndim, @as(usize, 3)) + (if (ndim > 3) @as(usize, 5) else @as(usize, 0));
    const pop_width = @as(usize, 11) * @min(nstate, @as(usize, 3)) + (if (nstate > 3) @as(usize, 5) else @as(usize, 0));

    if (has_thermo) {
        const fmt = "\n{[0]s:8} {[1]s:12} {[2]s:12} {[3]s:12} {[4]s:12} {[5]s:[6]} {[7]s:[8]} {[9]s:[10]} {[11]s:4}\n";

        const tuple = .{
            "ITER",

            "EKIN (Eh)",
            "EPOT (Eh)",
            "ETOT (Eh)",
            "TEMP (K)",

            "POS (a0)",
            col_width,

            "MOM (hb/a0)",
            col_width,

            "POP (-)",
            pop_width,

            "TIME",
        };

        try printf(io, fmt, tuple);
    }

    if (!has_thermo) {
        const fmt = "\n{[0]s:8} {[1]s:12} {[2]s:12} {[3]s:12} {[4]s:[5]} {[6]s:[7]} {[8]s:[9]} {[10]s:4}\n";

        const tuple = .{
            "ITER",

            "EKIN (Eh)",
            "EPOT (Eh)",
            "ETOT (Eh)",

            "POS (a0)",
            col_width,

            "MOM (hb/a0)",
            col_width,

            "POP (-)",
            pop_width,

            "TIME",
        };

        try printf(io, fmt, tuple);
    }
}

/// Prints the current iteration step's physical observables and elapsed time.
fn printIteration(comptime T: type, io: std.Io, obs: Observables(T), i: usize, has_thermo: bool, timer: *std.Io.Timestamp) !void {
    const ekin = obs.ekin orelse std.math.nan(T);
    const epot = obs.epot orelse std.math.nan(T);

    const etot = ekin + epot;

    try printf(io, "{d:8} {d:12.6} {d:12.6} {d:12.6} ", .{ i, ekin, epot, etot });

    if (has_thermo) if (obs.temp) |temp| {
        try printf(io, "{d:12.2} ", .{temp});
    };

    if (obs.pos) |pos_vec| {
        try printf(io, "[", .{});

        const n_show = @min(pos_vec.length(), 3);

        for (0..n_show) |j| {
            try printf(io, "{d:10.4}{s}", .{ pos_vec.at(j), if (j == pos_vec.length() - 1) "" else ", " });
        }

        if (pos_vec.length() > 3) {
            try printf(io, "...", .{});
        }

        try printf(io, "] ", .{});
    }

    if (obs.mom) |mom_vec| {
        try printf(io, "[", .{});

        const n_show = @min(mom_vec.length(), 3);

        for (0..n_show) |j| {
            try printf(io, "{d:10.4}{s}", .{ mom_vec.at(j), if (j == mom_vec.length() - 1) "" else ", " });
        }

        if (mom_vec.length() > 3) {
            try printf(io, "...", .{});
        }

        try printf(io, "] ", .{});
    }

    if (obs.pop) |pop_vec| {
        try printf(io, "[", .{});

        const n_show = @min(pop_vec.length(), 3);

        for (0..n_show) |j| {
            try printf(io, "{d:9.4}{s}", .{ pop_vec.at(j), if (j == pop_vec.length() - 1) "" else ", " });
        }

        if (pop_vec.length() > 3) {
            try printf(io, "...", .{});
        }

        try printf(io, "] ", .{});
    }

    try printf(io, "{f}\n", .{timer.untilNow(io, .real)});

    timer.* = std.Io.Timestamp.now(io, .real);
}

/// Propagates the classical equations of motion over the specified number of time steps.
fn solve(comptime T: type, io: std.Io, ctx: SolveContext(T), gpa: Allocator, _: Allocator) !Observables(T) {
    const ndim, const nstate = .{ ctx.sim.elpoten.ndim(), ctx.sim.elpoten.nstate() };

    const has_thermo = ctx.opt.thermostat != null;

    if (ctx.log) try printHeader(io, ndim, nstate, has_thermo);

    const atoms = switch (ctx.sim.elpoten) {
        .ab_initio => |ab| ab.msys.atoms,

        inline else => &.{},
    };

    const nrow, const iters = .{ ctx.sim.ensemble.r.nrow(), ctx.opt.iterations + 1 };

    var stream_writer = try StreamWriter.init(io, ctx.opt.stream, atoms, ndim, nrow, iters, nstate, gpa);
    defer if (stream_writer) |*sw| sw.deinit(io, gpa);

    var hist = try History(T).init(ndim, nstate, ctx.opt.iterations + 1, ctx.opt.write, gpa);
    defer hist.deinit(gpa);

    var timer = std.Io.Timestamp.now(io, .real);

    for (0..ctx.opt.iterations + 1) |i| {
        const time = @as(T, @floatFromInt(i)) * ctx.opt.time_step;

        if (i > 0) {
            try ctx.sim.propag.step(&ctx.sim.ensemble, &ctx.sim.gb, ctx.sim.elpoten, time);
        }

        if (stream_writer) |sw| {
            const prop = ctx.sim.propag;

            const coefs = if (prop.nonadia_dynamic) |n| (if (n == .ehrenfest) n.ehrenfest.coefics else null) else null;

            const is_sh = if (prop.nonadia_dynamic) |nd| nd == .surface_hopping else false;

            const s = if (is_sh) ctx.sim.ensemble.s else null;

            try sw.writeFrame(io, T, atoms, ctx.sim.ensemble.r, s, time);

            try sw.writeMomentumFrame(io, T, atoms, ctx.sim.ensemble.p, s, time);

            try sw.writeStatePotentialEnergy(io, T, ctx.sim.gb, nstate, time);

            try sw.writeEnergy(io, T, ctx.sim.ensemble, ctx.sim.gb, nstate, coefs, time);

            try sw.writePopulation(io, T, ctx.sim.ensemble, nstate, coefs, time);

            try sw.writeStateCoefficient(io, T, nstate, prop.coefs(), time);
        }

        const is_log_step = ctx.log and ((i % ctx.opt.log_interval == 0) or (i == ctx.opt.iterations));

        var obs = try Observables(T).init(ctx.sim.*, ctx.opt.write, is_log_step, has_thermo, gpa);
        defer obs.deinit(gpa);

        hist.append(obs);

        if (is_log_step) {
            try printIteration(T, io, obs, i, has_thermo, &timer);
        }
    }

    try hist.exportWrite(io, atoms, ctx.opt.time_step, ctx.opt.write);

    return try Observables(T).init(ctx.sim.*, ctx.opt.write, true, has_thermo, gpa);
}
