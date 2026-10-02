//! Classical molecular dynamics trajectory simulation on adiabatic or diabatic potential energy surfaces.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;

const AbInitioWrite = @import("potential.zig").AbInitioWrite;
const Ehrenfest = @import("ehrenfest.zig").Ehrenfest;
const EhrenfestOptions = @import("ehrenfest.zig").Options;
const Matrix = @import("tensor.zig").Matrix;
const Potential = @import("potential.zig").Potential;
const PotentialOptions = @import("potential.zig").Options;
const ScalarDual = @import("dual.zig").ScalarDual;
const SurfaceHopping = @import("surface_hopping.zig").SurfaceHopping;
const SurfaceHoppingOptions = @import("surface_hopping.zig").Options;
const Thermostat = @import("thermostat.zig").Thermostat;
const ThermostatOptions = @import("thermostat.zig").Options;
const Vector = @import("tensor.zig").Vector;

const eighBatch = @import("linear_algebra.zig").eighBatch;
const eighSlice = @import("linear_algebra.zig").eighSlice;
const getMass = @import("constant.zig").getMass;
const getTrajectoryPath = @import("read_write.zig").getTrajectoryPath;
const norm = @import("linear_algebra.zig").norm;
const printf = @import("read_write.zig").printf;
const writeMatrixLspace = @import("read_write.zig").writeMatrixLspace;
const writeXyzFrame = @import("read_write.zig").writeXyzFrame;

const AMU2AU = @import("constant.zig").AMU2AU;
const AU2K = @import("constant.zig").AU2K;

/// Configuration options for the classical molecular dynamics simulation.
pub const Options = struct {
    adiabatic: bool = true,
    initial_conditions: InitialConditions = .{},
    iterations: u32,
    log_interval: u32 = 1,
    mass: ?[]const f64 = null,
    nonadiabatic: ?NonadiabaticOptions = null,
    potential: PotentialOptions,
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

/// Initial phase space parameters and Gaussian width for trajectory sampling.
const InitialConditions = struct {
    gamma: ?[]const f64 = null,
    momentum: ?[]const f64 = null,
    position: ?[]const f64 = null,
    seed: u32 = 1,
    state: u32 = 0,
    temperature: ?f64 = null,
};

/// Output paths for recording trajectory observables to disk during dynamics.
const Write = struct {
    kinetic_energy: ?[]const u8 = null,
    momentum: ?[]const u8 = null,
    population: ?[]const u8 = null,
    position: ?[]const u8 = null,
    potential_energy: ?[]const u8 = null,
    state_potential_energy: ?[]const u8 = null,
    temperature: ?[]const u8 = null,
    total_energy: ?[]const u8 = null,
};

/// Generates a representation of a classical trajectory ensemble with positions, momenta, and active states.
pub fn Ensemble(comptime T: type) type {
    return struct {
        r: Matrix(T),
        p: Matrix(T),
        a: Matrix(T),
        m: []const T,

        s: Vector(usize),

        /// Initializes the ensemble trajectories with allocated memory for positions, momenta, and forces.
        pub fn init(ndim: usize, ntraj: usize, mass: []const T, gpa: Allocator) !@This() {
            var r = try Matrix(T).init(ntraj, ndim, gpa);
            errdefer r.deinit(gpa);

            var p = try Matrix(T).init(ntraj, ndim, gpa);
            errdefer p.deinit(gpa);

            var a = try Matrix(T).init(ntraj, ndim, gpa);
            errdefer a.deinit(gpa);

            var s = try Vector(usize).init(ntraj, gpa);
            errdefer s.deinit(gpa);

            const m = try gpa.alloc(T, mass.len);
            errdefer gpa.free(m);

            @memcpy(m, mass);

            return .{ .r = r, .p = p, .a = a, .s = s, .m = m };
        }

        /// Deallocates the positions, momenta, and force vectors of the trajectory ensemble.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.r.deinit(gpa);
            self.p.deinit(gpa);
            self.a.deinit(gpa);
            self.s.deinit(gpa);

            gpa.free(self.m);
        }

        /// Calculates the average classical kinetic energy of the ensemble trajectories.
        pub fn ekin(self: @This()) T {
            var sum: T = 0;

            for (0..self.p.nrow()) |i| for (0..self.p.ncol()) |j| {
                sum += self.p.at(i, j) * self.p.at(i, j) / (2 * self.m[j]);
            };

            return sum / @as(T, @floatFromInt(self.p.nrow()));
        }

        /// Computes the mean momentum vector averaged over all trajectories in the ensemble.
        pub fn mom(self: @This(), gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(self.p.ncol(), gpa);

            for (0..self.p.nrow()) |i| for (0..self.p.ncol()) |j| {
                value.ptr(j).* += self.p.at(i, j);
            };

            value.divs(@floatFromInt(self.p.nrow()));

            return value;
        }

        /// Computes the population fraction of trajectories occupying each electronic state.
        pub fn pop(self: @This(), nstate: usize, coefs: ?Matrix(Complex(T)), U: ?Matrix(T), adia: bool, gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(nstate, gpa);

            if (coefs) |coef| {
                if (adia) for (0..self.r.nrow()) |i| for (0..nstate) |m| {
                    var a_im = Complex(T).init(0, 0);

                    for (0..nstate) |k| {
                        const u_km, const c_ik = .{ U.?.at(i, k * nstate + m), coef.at(i, k) };

                        a_im = a_im.add(Complex(T).init(c_ik.re * u_km, c_ik.im * u_km));
                    }

                    value.ptr(m).* += a_im.squaredMagnitude();
                };

                if (!adia) for (0..coef.nrow()) |i| for (0..nstate) |k| {
                    value.ptr(k).* += coef.at(i, k).squaredMagnitude();
                };

                value.divs(@floatFromInt(self.r.nrow()));
            }

            if (coefs == null) {
                for (0..self.s.length()) |i| {
                    value.ptr(self.s.at(i)).* += 1;
                }

                value.divs(@floatFromInt(self.s.length()));
            }

            return value;
        }

        /// Computes the mean position vector averaged over all trajectories in the ensemble.
        pub fn pos(self: @This(), gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(self.r.ncol(), gpa);

            for (0..self.r.nrow()) |i| for (0..self.r.ncol()) |j| {
                value.ptr(j).* += self.r.at(i, j);
            };

            value.divs(@floatFromInt(self.r.nrow()));

            return value;
        }

        /// Samples initial positions and momenta from a Wigner-like Gaussian distribution or reference geometry.
        pub fn setGaussian(self: *@This(), ic: InitialConditions, default_pos: ?[]const T) void {
            var split_mix = std.Random.SplitMix64.init(ic.seed);

            var rng = std.Random.DefaultPrng.init(split_mix.next());

            for (0..self.s.length()) |i| {
                self.s.ptr(i).* = ic.state;
            }

            const random = rng.random();

            for (0..self.r.nrow()) |i| for (0..self.r.ncol()) |j| {
                const r0 = if (ic.position) |r_pos| r_pos[j] else default_pos.?[j];

                if (ic.gamma) |gamma| {
                    const stdev = 1 / std.math.sqrt(2 * gamma[j]);

                    self.r.ptr(i, j).* = r0 + stdev * random.floatNorm(T);
                }

                if (ic.gamma == null) {
                    self.r.ptr(i, j).* = r0;
                }
            };

            for (0..self.p.nrow()) |i| for (0..self.p.ncol()) |j| {
                const p0 = if (ic.momentum) |p_mom| p_mom[j] else 0;

                if (ic.temperature) |t_val| {
                    const stdev = std.math.sqrt(self.m[j] * @as(T, @floatCast(t_val)) / AU2K);

                    self.p.ptr(i, j).* = p0 + stdev * random.floatNorm(T);

                    continue;
                }

                if (ic.gamma) |gamma| {
                    const stdev = std.math.sqrt(gamma[j] / 2);

                    self.p.ptr(i, j).* = p0 + stdev * random.floatNorm(T);

                    continue;
                }

                self.p.ptr(i, j).* = p0;
            };
        }

        /// Calculates the instantaneous kinetic temperature of the ensemble in Kelvin.
        pub fn temp(self: @This()) T {
            return (2 * self.ekin() / @as(T, @floatFromInt(self.p.ncol()))) * AU2K;
        }
    };
}

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

/// Manages streaming output files for multi-trajectory ab initio property logging.
const AbInitioStreamWriter = struct {
    files_pos: ?[]std.Io.File = null,

    /// Initializes and opens trajectory output files on disk for configured properties.
    pub fn init(io: std.Io, write: AbInitioWrite, ntraj: usize, gpa: Allocator) !@This() {
        var self = @This(){};
        errdefer self.deinit(io, gpa);

        if (write.position) |path| {
            self.files_pos = try initFiles(io, path, ntraj, gpa);
        }

        return self;
    }

    /// Closes all trajectory output files and deallocates handles.
    pub fn deinit(self: *@This(), io: std.Io, gpa: Allocator) void {
        if (self.files_pos) |files| {
            for (files) |*file| file.close(io);

            gpa.free(files);
        }
    }

    /// Appends the current step's nuclear coordinates to trajectory files.
    pub fn writeFrame(self: @This(), io: std.Io, comptime T: type, atoms: []const i32, r: Matrix(T), time: T) !void {
        if (self.files_pos) |files| {
            var buffer: [65536]u8 = undefined;

            for (0..files.len) |k| {
                var writer = files[k].writerStreaming(io, &buffer);

                try writeXyzFrame(T, &writer, atoms, r.rowSlice(k), time);

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
};

/// Helper struct managing memory for potential energy gradients and wavefunctions.
fn GradientBuffer(comptime T: type) type {
    return struct {
        adia: bool,

        r_dual: Matrix(ScalarDual(T)),
        V_dual: Matrix(ScalarDual(T)),

        V: Matrix(T),
        W: Matrix(T),
        U: Matrix(T),

        grad_V: Matrix(T),

        /// Allocates memory for dual numbers and gradients needed in force evaluations.
        pub fn init(ndim: usize, nstate: usize, ntraj: usize, adia: bool, gpa: Allocator) !@This() {
            var r_dual = try Matrix(ScalarDual(T)).init(ntraj, ndim, gpa);
            errdefer r_dual.deinit(gpa);

            var V_dual = try Matrix(ScalarDual(T)).init(ntraj, nstate * nstate, gpa);
            errdefer V_dual.deinit(gpa);

            var V = try Matrix(T).init(ntraj, nstate * nstate, gpa);
            errdefer V.deinit(gpa);

            var W = try Matrix(T).init(ntraj, nstate, gpa);
            errdefer W.deinit(gpa);

            var U = try Matrix(T).init(ntraj, nstate * nstate, gpa);
            errdefer U.deinit(gpa);

            const grad_V = try Matrix(T).init(ntraj, ndim * nstate * nstate, gpa);
            errdefer grad_V.deinit(gpa);

            return .{ .r_dual = r_dual, .V_dual = V_dual, .V = V, .W = W, .U = U, .grad_V = grad_V, .adia = adia };
        }

        /// Deallocates gradient matrices, eigenvalues, and eigenvectors.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.r_dual.deinit(gpa);
            self.V_dual.deinit(gpa);

            self.V.deinit(gpa);
            self.W.deinit(gpa);
            self.U.deinit(gpa);

            self.grad_V.deinit(gpa);
        }

        /// Computes and applies nuclear forces to trajectories based on potential gradients.
        pub fn apply(self: *@This(), ensemble: *Ensemble(T), pot: Potential(T), coefs: ?*const Matrix(Complex(T))) void {
            const nstate = pot.nstate();

            for (0..ensemble.r.nrow()) |i| for (0..ensemble.r.ncol()) |j| {
                if (coefs) |coef| {
                    var der: T = 0;

                    for (0..nstate) |k| for (0..nstate) |l| {
                        const rho_kl = coef.at(i, k).conjugate().mul(coef.at(i, l)).re;

                        der += rho_kl * self.grad_V.at(i, j * nstate * nstate + k * nstate + l);
                    };

                    ensemble.a.ptr(i, j).* = -der / ensemble.m[j];
                }

                if (coefs == null and self.adia) {
                    var der: T = 0;

                    for (0..nstate) |k| {
                        const U_ka = self.U.at(i, k * nstate + ensemble.s.at(i));

                        for (0..nstate) |l| {
                            const U_la = self.U.at(i, l * nstate + ensemble.s.at(i));

                            der += U_ka * U_la * self.grad_V.at(i, j * nstate * nstate + k * nstate + l);
                        }
                    }

                    ensemble.a.ptr(i, j).* = -der / ensemble.m[j];
                }

                if (coefs == null and !self.adia) {
                    const k = ensemble.s.at(i) * nstate + ensemble.s.at(i);

                    ensemble.a.ptr(i, j).* = -self.grad_V.at(i, j * nstate * nstate + k) / ensemble.m[j];
                }
            };
        }

        /// Calculates average potential energy from cached gradient buffer matrices.
        pub fn epot(self: @This(), ensemble: Ensemble(T), nstate: usize, coefs: ?Matrix(Complex(T))) T {
            var sum: T = 0;

            if (coefs) |coef| for (0..ensemble.r.nrow()) |i| {
                var traj_epot: T = 0;

                for (0..nstate) |k| for (0..nstate) |l| {
                    const rho_kl = coef.at(i, k).conjugate().mul(coef.at(i, l)).re;

                    traj_epot += rho_kl * self.V.at(i, k * nstate + l);
                };

                sum += traj_epot;
            };

            if (coefs == null and self.adia) for (0..ensemble.r.nrow()) |i| {
                sum += self.W.at(i, ensemble.s.at(i));
            };

            if (coefs == null and !self.adia) for (0..ensemble.r.nrow()) |i| {
                sum += self.V.at(i, ensemble.s.at(i) * nstate + ensemble.s.at(i));
            };

            return sum / @as(T, @floatFromInt(ensemble.r.nrow()));
        }

        /// Computes the ensemble-averaged potential energy for each electronic state.
        pub fn epotState(self: @This(), nstate: usize, gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(nstate, gpa);

            for (0..self.V.nrow()) |i| for (0..nstate) |k| {
                const v_k = if (self.adia) self.W.at(i, k) else self.V.at(i, k * nstate + k);

                value.ptr(k).* += v_k;
            };

            value.divs(@floatFromInt(self.V.nrow()));

            return value;
        }

        /// Evaluates the potential and its gradients using automatic differentiation or direct ab initio evaluation.
        pub fn update(self: *@This(), r: Matrix(T), pot: Potential(T), time: T, states: []const usize) !void {
            switch (pot) {
                .ab_initio => |ab| {
                    try ab.updateGradients(T, &self.V, &self.grad_V, r, states);

                    if (self.adia) for (0..self.V.nrow()) |i| for (0..pot.nstate()) |k| {
                        self.W.ptr(i, k).* = self.V.at(i, k * pot.nstate() + k);

                        for (0..pot.nstate()) |l| {
                            self.U.ptr(i, k * pot.nstate() + l).* = if (k == l) 1 else 0;
                        }
                    };
                },

                inline else => {
                    try pot.evalBatch(T, &self.V, r, time);

                    if (self.adia) {
                        try eighBatch(T, &self.W, &self.U, self.V);
                    }

                    for (0..r.nrow()) |i| for (0..r.ncol()) |j| {
                        for (0..r.ncol()) |k| {
                            self.r_dual.ptr(i, k).* = ScalarDual(T).init(r.at(i, k), if (k == j) 1 else 0);
                        }

                        try pot.eval(ScalarDual(T), self.V_dual.rowSlice(i), self.r_dual.rowSlice(i), ScalarDual(T).init(time, 0));

                        for (0..pot.nstate() * pot.nstate()) |k| {
                            self.grad_V.ptr(i, j * pot.nstate() * pot.nstate() + k).* = self.V_dual.at(i, k).der;
                        }
                    };
                },
            }
        }
    };
}

/// Accumulates time-dependent expectation values and dynamics trajectories.
fn History(comptime T: type) type {
    return struct {
        pos: ?Matrix(T) = null,
        mom: ?Matrix(T) = null,
        pop: ?Matrix(T) = null,

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
        pub fn exportWrite(self: *@This(), io: std.Io, dt: f64, write: Write) !void {
            const end = dt * @as(T, @floatFromInt(self.index - 1));

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

            if (write.state_potential_energy) |path| {
                try writeMatrixLspace(T, io, path, self.state_epot.?.takeRows(self.index), 0, end);
            }

            if (write.kinetic_energy) |path| {
                try writeMatrixLspace(T, io, path, self.ekin.?.takeRows(self.index), 0, end);
            }

            if (write.temperature) |path| {
                try writeMatrixLspace(T, io, path, self.temp.?.takeRows(self.index), 0, end);
            }

            if (write.total_energy) |path| {
                try writeMatrixLspace(T, io, path, self.etot.?.takeRows(self.index), 0, end);
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

            const coefs = if (sim.propag.nonadia_dynamic) |n| (if (n == .ehrenfest) n.ehrenfest.coefics else null) else null;

            if (calc.pop) obs.pop = try sim.ensemble.pop(sim.elpoten.nstate(), coefs, sim.gb.U, sim.gb.adia, gpa);

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
                const eh = try Ehrenfest(T).init(naopt.ehrenfest, nstate, opt.trajectories, gpa);

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

        /// Integrates equations of motion using the Velocity Verlet algorithm.
        pub fn step(self: *@This(), ens: *Ensemble(T), gb: *GradientBuffer(T), pot: Potential(T), time: T) !void {
            for (0..ens.r.nrow()) |i| for (0..ens.r.ncol()) |j| {
                ens.p.ptr(i, j).* += 0.5 * ens.m[j] * ens.a.at(i, j) * self.dt;

                ens.r.ptr(i, j).* += 0.5 * (ens.p.at(i, j) / ens.m[j]) * self.dt;
            };

            if (self.thermo) |*thermo| {
                thermo.apply(&ens.p, ens.m);
            }

            for (0..ens.r.nrow()) |i| for (0..ens.r.ncol()) |j| {
                ens.r.ptr(i, j).* += 0.5 * (ens.p.at(i, j) / ens.m[j]) * self.dt;
            };

            try gb.update(ens.r, pot, time, ens.s.data);

            if (self.nonadia_dynamic) |*n| if (n.* == .ehrenfest) {
                try n.ehrenfest.step(gb.V, self.dt);
            };

            const coefics = if (self.nonadia_dynamic) |*n| (if (n.* == .ehrenfest) &n.ehrenfest.coefics else null) else null;

            gb.apply(ens, pot, coefics);

            for (0..ens.r.nrow()) |i| for (0..ens.r.ncol()) |j| {
                ens.p.ptr(i, j).* += 0.5 * ens.m[j] * ens.a.at(i, j) * self.dt;
            };

            if (self.nonadia_dynamic) |*n| if (n.* == .surface_hopping) {
                if (try n.surface_hopping.hop(ens, gb.V, gb.W, gb.U, self.dt)) {
                    for (0..ens.s.length()) |i| if (ens.s.at(i) != n.surface_hopping.targets[i]) {
                        switch (pot) {
                            .ab_initio => |ab| {
                                try ab.updateTrajectoryGradient(T, &gb.V, &gb.grad_V, ens.r, i, ens.s.at(i));
                            },

                            inline else => {},
                        }
                    };
                }

                gb.apply(ens, pot, null);
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
            inline for (@typeInfo(@This()).@"struct".fields) |field| {
                @field(self, field.name).deinit(gpa);
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
    try checkInvalidInput(opt);

    if (log) try std.Io.File.stdout().writeStreamingAll(io, "\nCLASSICAL DYNAMICS INIT: ");

    var timer = std.Io.Timestamp.now(io, .real);

    var sim = try init(T, io, opt, gpa);
    defer sim.deinit(gpa);

    if (log) try printf(io, "{f}\n", .{timer.untilNow(io, .real)});

    var obs = try solve(T, io, .{ .opt = opt, .sim = &sim, .log = log }, gpa, gpa);
    errdefer obs.deinit(gpa);

    if (log) {
        try printFinalPop(T, io, obs);
    }

    return .{ .observables = obs };
}

/// Validates simulation options to ensure positive mass, time step, and consistent dimensions.
fn checkInvalidInput(opt: Options) !void {
    if (opt.time_step <= 0) {
        std.log.err("TIME STEP MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.trajectories == 0) {
        std.log.err("NUMBER OF TRAJECTORIES MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.log_interval == 0) {
        std.log.err("LOG INTERVAL MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.mass) |mass_vec| {
        if (opt.potential == .ab_initio) {
            std.log.err("MASS VECTOR MUST NOT BE PROVIDED FOR AB INITIO POTENTIAL", .{});

            return error.InvalidInput;
        }

        for (mass_vec) |m| if (m <= 0) {
            std.log.err("MASS MUST BE GREATER THAN 0", .{});

            return error.InvalidInput;
        };

        if (opt.initial_conditions.position) |pos| if (mass_vec.len != pos.len) {
            std.log.err("MASS VECTOR MUST HAVE THE SAME LENGTH AS POSITION VECTOR", .{});

            return error.InvalidInput;
        };
    }

    if (opt.mass == null and opt.potential != .ab_initio) {
        std.log.err("MASS VECTOR MUST BE PROVIDED FOR MODEL POTENTIALS", .{});

        return error.InvalidInput;
    }

    if (opt.initial_conditions.position) |pos| {
        if (pos.len == 0) {
            std.log.err("INITIAL POSITION VECTOR MUST NOT BE EMPTY", .{});

            return error.InvalidInput;
        }

        if (opt.initial_conditions.momentum) |mom| if (mom.len != pos.len) {
            std.log.err("INITIAL MOMENTUM AND POSITION VECTORS MUST HAVE THE SAME LENGTH", .{});

            return error.InvalidInput;
        };

        if (opt.initial_conditions.gamma) |gam| if (gam.len != pos.len) {
            std.log.err("INITIAL GAMMA VECTOR MUST HAVE THE SAME LENGTH AS POSITION VECTOR", .{});

            return error.InvalidInput;
        };
    }

    if (opt.initial_conditions.position == null and opt.potential != .ab_initio) {
        std.log.err("INITIAL POSITION VECTOR MUST BE PROVIDED FOR MODEL POTENTIALS", .{});

        return error.InvalidInput;
    }

    if (opt.potential == .ab_initio and opt.nonadiabatic != null) {
        switch (opt.nonadiabatic.?) {
            .ehrenfest => {
                std.log.err("EHRENFEST DYNAMICS IS NOT SUPPORTED FOR AB INITIO POTENTIAL WITHOUT NACVS", .{});

                return error.InvalidInput;
            },
            .surface_hopping => |sh| switch (sh) {
                .fewest_switches => {
                    std.log.err("FSSH IS NOT SUPPORTED FOR AB INITIO POTENTIAL WITHOUT NACVS", .{});

                    return error.InvalidInput;
                },
                inline else => {},
            },
        }
    }

    if (opt.thermostat) |topt| switch (topt) {
        .berendsen => |bopt| {
            if (bopt.temperature < 0) {
                std.log.err("TEMPERATURE MUST BE NON-NEGATIVE", .{});

                return error.InvalidInput;
            }

            if (bopt.tau <= 0) {
                std.log.err("COUPLING CONSTANT TAU MUST BE POSITIVE", .{});

                return error.InvalidInput;
            }
        },
        .langevin => |lopt| {
            if (lopt.temperature < 0) {
                std.log.err("TEMPERATURE MUST BE NON-NEGATIVE", .{});

                return error.InvalidInput;
            }

            if (lopt.gamma < 0) {
                std.log.err("FRICTION COEFFICIENT MUST BE NON-NEGATIVE", .{});

                return error.InvalidInput;
            }
        },
    };
}

/// Initializes the simulation state, potential, and initial ensemble.
fn init(comptime T: type, io: std.Io, opt: Options, gpa: Allocator) !SimulationState(T) {
    var pot = try Potential(T).init(io, opt.potential, gpa);
    errdefer pot.deinit(gpa);

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

    ensemble.setGaussian(opt.initial_conditions, default_pos);

    try gb.update(ensemble.r, pot, 0, ensemble.s.data);

    if (prop.nonadia_dynamic) |*n| if (n.* == .surface_hopping) {
        n.surface_hopping.update(if (opt.adiabatic) gb.W else gb.V, gb.U);
    };

    if (prop.nonadia_dynamic) |*n| if (n.* == .ehrenfest) {
        n.ehrenfest.setInitialState(opt.initial_conditions.state, opt.adiabatic, gb.U);
    };

    const coefs = if (prop.nonadia_dynamic) |*n| (if (n.* == .ehrenfest) &n.ehrenfest.coefics else null) else null;

    gb.apply(&ensemble, pot, coefs);

    return .{ .ensemble = ensemble, .elpoten = pot, .propag = prop, .gb = gb };
}

/// Prints the final population fractions of electronic states to stdout.
fn printFinalPop(comptime T: type, io: std.Io, obs: Observables(T)) !void {
    if (obs.pop) |pop| {
        try std.Io.File.stdout().writeStreamingAll(io, "\n");

        for (0..pop.length()) |i| {
            try printf(io, "FINAL POPULATION OF ELECTRONIC STATE {d:02}: {d:.8}\n", .{ i, pop.at(i) });
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

    var stream_writer: ?AbInitioStreamWriter = null;

    if (ctx.sim.elpoten == .ab_initio) {
        const ab = ctx.sim.elpoten.ab_initio;

        if (ab.options.write.any()) {
            stream_writer = try AbInitioStreamWriter.init(io, ab.options.write, ctx.sim.ensemble.r.nrow(), gpa);
        }
    }

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
            try sw.writeFrame(io, T, atoms, ctx.sim.ensemble.r, time);
        }

        const is_log_step = ctx.log and ((i % ctx.opt.log_interval == 0) or (i == ctx.opt.iterations));

        var obs = try Observables(T).init(ctx.sim.*, ctx.opt.write, is_log_step, has_thermo, gpa);
        defer obs.deinit(gpa);

        hist.append(obs);

        if (is_log_step) {
            try printIteration(T, io, obs, i, has_thermo, &timer);
        }
    }

    try hist.exportWrite(io, ctx.opt.time_step, ctx.opt.write);

    return try Observables(T).init(ctx.sim.*, ctx.opt.write, true, has_thermo, gpa);
}
