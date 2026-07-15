//! Classical molecular dynamics trajectory simulation on adiabatic or diabatic potential energy surfaces.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;

const Ehrenfest = @import("ehrenfest.zig").Ehrenfest;
const EhrenfestOptions = @import("ehrenfest.zig").Options;
const Matrix = @import("tensor.zig").Matrix;
const Potential = @import("potential.zig").Potential;
const PotentialOptions = @import("potential.zig").Options;
const ScalarDual = @import("dual.zig").ScalarDual;
const SurfaceHopping = @import("surface_hopping.zig").SurfaceHopping;
const SurfaceHoppingOptions = @import("surface_hopping.zig").Options;
const Vector = @import("tensor.zig").Vector;

const eighBatch = @import("linear_algebra.zig").eighBatch;
const eighSlice = @import("linear_algebra.zig").eighSlice;
const norm = @import("linear_algebra.zig").norm;
const printf = @import("read_write.zig").printf;
const writeMatrixLspace = @import("read_write.zig").writeMatrixLspace;

/// Tagged union for multi-state non-adiabatic trajectory propagation methods.
pub const NonadiabaticOptions = union(enum) {
    surface_hopping: SurfaceHoppingOptions,
    ehrenfest: EhrenfestOptions,
};

/// Configuration options for the classical molecular dynamics simulation.
pub const Options = struct {
    initial_conditions: InitialConditions,
    potential: PotentialOptions,

    time_step: f64,
    iterations: u32,
    mass: []const f64,
    trajectories: u32,

    nonadiabatic: ?NonadiabaticOptions = null,
    write: Write = .{},

    adiabatic: bool = true,
    log_interval: u32 = 1,
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

        /// Calculates the average potential energy of the ensemble based on active electronic states.
        pub fn epot(self: @This(), pot: Potential(T), time: T, adiabatic: bool, coefics: ?Matrix(Complex(T)), gpa: Allocator) !T {
            var sum: T = 0;

            const V_slice = try gpa.alloc(T, pot.nstate() * pot.nstate());
            defer gpa.free(V_slice);

            const U_slice = try gpa.alloc(T, pot.nstate() * pot.nstate());
            defer gpa.free(U_slice);

            const W_slice = try gpa.alloc(T, pot.nstate());
            defer gpa.free(W_slice);

            if (coefics) |coef| for (0..self.r.nrow()) |i| {
                pot.eval(T, V_slice, self.r.rowSlice(i), time);

                var traj_epot: T = 0;

                for (0..pot.nstate()) |k| for (0..pot.nstate()) |l| {
                    const rho_kl = coef.at(i, k).conjugate().mul(coef.at(i, l)).re;

                    traj_epot += rho_kl * V_slice[k * pot.nstate() + l];
                };

                sum += traj_epot;
            };

            if (coefics == null and adiabatic) for (0..self.r.nrow()) |i| {
                pot.eval(T, V_slice, self.r.rowSlice(i), time);

                try eighSlice(T, W_slice, U_slice, V_slice);

                sum += W_slice[self.s.at(i)];
            };

            if (coefics == null and !adiabatic) for (0..self.r.nrow()) |i| {
                pot.eval(T, V_slice, self.r.rowSlice(i), time);

                sum += V_slice[self.s.at(i) * pot.nstate() + self.s.at(i)];
            };

            return sum / @as(T, @floatFromInt(self.r.nrow()));
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
        pub fn pop(self: @This(), nstate: usize, coefics: ?Matrix(Complex(T)), U: ?Matrix(T), adia: bool, gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(nstate, gpa);

            if (coefics) |coef| {
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

            if (coefics == null) {
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

        /// Samples initial positions and momenta from a Wigner-like Gaussian distribution.
        pub fn setGaussian(self: *@This(), ic: InitialConditions) void {
            var split_mix = std.Random.SplitMix64.init(ic.seed);
            var rng = std.Random.DefaultPrng.init(split_mix.next());

            for (0..self.s.length()) |i| {
                self.s.ptr(i).* = ic.state;
            }

            const random = rng.random();

            for (0..self.r.nrow()) |i| for (0..self.r.ncol()) |j| {
                const stdev = 1 / std.math.sqrt(2 * ic.gamma[j]);

                self.r.ptr(i, j).* = ic.position[j] + stdev * random.floatNorm(T);
            };

            for (0..self.p.nrow()) |i| for (0..self.p.ncol()) |j| {
                const stdev = std.math.sqrt(ic.gamma[j] / 2);

                self.p.ptr(i, j).* = ic.momentum[j] + stdev * random.floatNorm(T);
            };
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

/// Initial phase space parameters and Gaussian width for trajectory sampling.
const InitialConditions = struct {
    position: []const f64,
    momentum: []const f64,
    gamma: []const f64,

    state: u32 = 0,
    seed: u32 = 1,
};

/// Output paths for recording trajectory observables to disk during dynamics.
const Write = struct {
    kinetic_energy: ?[]const u8 = null,
    momentum: ?[]const u8 = null,
    population: ?[]const u8 = null,
    position: ?[]const u8 = null,
    potential_energy: ?[]const u8 = null,
    total_energy: ?[]const u8 = null,
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
        pub fn apply(self: *@This(), ensemble: *Ensemble(T), pot: Potential(T), coefics: ?*const Matrix(Complex(T))) void {
            const nstate = pot.nstate();

            for (0..ensemble.r.nrow()) |i| for (0..ensemble.r.ncol()) |j| {
                if (coefics) |coef| {
                    var der: T = 0;

                    for (0..nstate) |k| for (0..nstate) |l| {
                        const rho_kl = coef.at(i, k).conjugate().mul(coef.at(i, l)).re;

                        der += rho_kl * self.grad_V.at(i, j * nstate * nstate + k * nstate + l);
                    };

                    ensemble.a.ptr(i, j).* = -der / ensemble.m[j];
                }

                if (coefics == null and self.adia) {
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

                if (coefics == null and !self.adia) {
                    const k = ensemble.s.at(i) * nstate + ensemble.s.at(i);

                    ensemble.a.ptr(i, j).* = -self.grad_V.at(i, j * nstate * nstate + k) / ensemble.m[j];
                }
            };
        }

        /// Evaluates the potential and its gradients using automatic differentiation.
        pub fn update(self: *@This(), r: Matrix(T), pot: Potential(T), time: T) !void {
            pot.evalBatch(T, &self.V, r, time);

            if (self.adia) {
                try eighBatch(T, &self.W, &self.U, self.V);
            }

            for (0..r.nrow()) |i| for (0..r.ncol()) |j| {
                for (0..r.ncol()) |k| {
                    self.r_dual.ptr(i, k).* = ScalarDual(T).init(r.at(i, k), if (k == j) 1 else 0);
                }

                pot.eval(ScalarDual(T), self.V_dual.rowSlice(i), self.r_dual.rowSlice(i), ScalarDual(T).init(time, 0));

                for (0..pot.nstate() * pot.nstate()) |k| {
                    self.grad_V.ptr(i, j * pot.nstate() * pot.nstate() + k).* = self.V_dual.at(i, k).der;
                }
            };
        }
    };
}

/// Accumulates time-dependent expectation values and dynamics trajectories.
fn History(comptime T: type) type {
    return struct {
        pos: ?Matrix(T) = null,
        mom: ?Matrix(T) = null,
        pop: ?Matrix(T) = null,

        epot: ?Matrix(T) = null,
        ekin: ?Matrix(T) = null,
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

            if (store_ekin) hist.ekin = try Matrix(T).init(iters, 1, gpa);
            if (store_epot) hist.epot = try Matrix(T).init(iters, 1, gpa);
            if (store_etot) hist.etot = try Matrix(T).init(iters, 1, gpa);

            return hist;
        }

        /// Deallocates history storage arrays.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            if (self.pos) |*pos| pos.deinit(gpa);
            if (self.mom) |*mom| mom.deinit(gpa);
            if (self.pop) |*pop| pop.deinit(gpa);

            if (self.epot) |*epot| epot.deinit(gpa);
            if (self.ekin) |*ekin| ekin.deinit(gpa);
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

            if (self.epot) |*epot| {
                epot.ptr(step_idx, 0).* = obs.epot.?;
            }

            if (self.ekin) |*ekin| {
                ekin.ptr(step_idx, 0).* = obs.ekin.?;
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

            if (write.kinetic_energy) |path| {
                try writeMatrixLspace(T, io, path, self.ekin.?.takeRows(self.index), 0, end);
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

        epot: ?T = null,
        ekin: ?T = null,

        /// Computes the physical observables from the current simulation state.
        pub fn init(sim: SimulationState(T), time: T, write: Write, log: bool, gpa: Allocator) !@This() {
            var obs = @This(){};
            errdefer obs.deinit(gpa);

            const calc_ekin, const calc_epot = .{ write.kinetic_energy != null, write.potential_energy != null };

            var calc = .{
                .pos = log or write.position != null,
                .mom = log or write.momentum != null,

                .pop = log or write.population != null,

                .ekin = log or calc_ekin,
                .epot = log or calc_epot,
            };

            calc.ekin = calc.ekin or write.total_energy != null;
            calc.epot = calc.epot or write.total_energy != null;

            if (calc.mom) obs.mom = try sim.ensemble.mom(gpa);
            if (calc.pos) obs.pos = try sim.ensemble.pos(gpa);

            const coefics = if (sim.propag.namd) |*n| (if (n.* == .ehrenfest) n.ehrenfest.coefics else null) else null;

            if (calc.pop) obs.pop = try sim.ensemble.pop(sim.elpoten.nstate(), coefics, sim.gb.U, sim.gb.adia, gpa);

            if (calc.ekin) {
                obs.ekin = sim.ensemble.ekin();
            }

            if (calc.epot) {
                obs.epot = try sim.ensemble.epot(sim.elpoten, time, sim.gb.adia, coefics, gpa);
            }

            return obs;
        }

        /// Deallocates arrays stored in the observables struct.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            if (self.pos) |*pos| pos.deinit(gpa);
            if (self.mom) |*mom| mom.deinit(gpa);
            if (self.pop) |*pop| pop.deinit(gpa);
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

        namd: ?Namd = null,

        dt: T,

        /// Initializes the time propagator and optional surface hopping solver.
        pub fn init(opt: Options, nstate: usize, gpa: Allocator) !@This() {
            const istate = opt.initial_conditions.state;

            var namd: ?Namd = null;

            if (opt.nonadiabatic) |naopt| if (naopt == .surface_hopping) {
                const sh = try SurfaceHopping(T).init(naopt.surface_hopping, nstate, opt.trajectories, istate, opt.adiabatic, gpa);

                namd = .{ .surface_hopping = sh };
            };

            if (opt.nonadiabatic) |naopt| if (naopt == .ehrenfest) {
                const eh = try Ehrenfest(T).init(naopt.ehrenfest, nstate, opt.trajectories, gpa);

                namd = .{ .ehrenfest = eh };
            };

            return .{ .dt = @floatCast(opt.time_step), .namd = namd };
        }

        /// Deallocates surface hopping resources.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            if (self.namd) |*namd| switch (namd.*) {
                inline else => |*n| n.deinit(gpa),
            };
        }

        /// Integrates equations of motion using the Velocity Verlet algorithm.
        pub fn step(self: *@This(), ens: *Ensemble(T), gb: *GradientBuffer(T), pot: Potential(T), time: T) !void {
            for (0..ens.r.nrow()) |i| for (0..ens.r.ncol()) |j| {
                ens.p.ptr(i, j).* += 0.5 * ens.m[j] * ens.a.at(i, j) * self.dt;

                ens.r.ptr(i, j).* += (ens.p.at(i, j) / ens.m[j]) * self.dt;
            };

            try gb.update(ens.r, pot, time);

            if (self.namd) |*n| if (n.* == .ehrenfest) {
                try n.ehrenfest.step(gb.V, self.dt);
            };

            const coefics = if (self.namd) |*n| (if (n.* == .ehrenfest) &n.ehrenfest.coefics else null) else null;

            gb.apply(ens, pot, coefics);

            for (0..ens.r.nrow()) |i| for (0..ens.r.ncol()) |j| {
                ens.p.ptr(i, j).* += 0.5 * ens.m[j] * ens.a.at(i, j) * self.dt;
            };

            if (self.namd) |*n| if (n.* == .surface_hopping) {
                try n.surface_hopping.hop(ens, gb.V, gb.W, gb.U, self.dt);

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

    for (opt.mass) |m| if (m <= 0) {
        std.log.err("MASS MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    };

    if (opt.mass.len != opt.initial_conditions.position.len) {
        std.log.err("MASS VECTOR MUST HAVE THE SAME LENGTH AS POSITION VECTOR", .{});

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
}

/// Initializes the simulation state, potential, and initial ensemble.
fn init(comptime T: type, io: std.Io, opt: Options, gpa: Allocator) !SimulationState(T) {
    var pot = try Potential(T).init(io, opt.potential, gpa);
    errdefer pot.deinit(gpa);

    const mass = try gpa.alloc(T, opt.mass.len);
    defer gpa.free(mass);

    for (opt.mass, 0..) |m, i| {
        mass[i] = @floatCast(m);
    }

    var ensemble = try Ensemble(T).init(pot.ndim(), opt.trajectories, mass, gpa);
    errdefer ensemble.deinit(gpa);

    var gb = try GradientBuffer(T).init(pot.ndim(), pot.nstate(), opt.trajectories, opt.adiabatic, gpa);
    errdefer gb.deinit(gpa);

    var prop = try Propagator(T).init(opt, pot.nstate(), gpa);
    errdefer prop.deinit(gpa);

    ensemble.setGaussian(opt.initial_conditions);

    try gb.update(ensemble.r, pot, 0);

    if (prop.namd) |*n| if (n.* == .surface_hopping) {
        n.surface_hopping.update(if (opt.adiabatic) gb.W else gb.V, gb.U);
    };

    if (prop.namd) |*n| if (n.* == .ehrenfest) {
        n.ehrenfest.setInitialState(opt.initial_conditions.state, opt.adiabatic, gb.U);
    };

    const coefics = if (prop.namd) |*n| (if (n.* == .ehrenfest) &n.ehrenfest.coefics else null) else null;

    gb.apply(&ensemble, pot, coefics);

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
fn printHeader(io: std.Io, ndim: usize, nstate: usize) !void {
    try std.Io.File.stdout().writeStreamingAll(io, "\nREAL-TIME PROPAGATION");

    const fmt = "\n{[0]s:8} {[1]s:12} {[2]s:12} {[3]s:12} {[4]s:[5]} {[6]s:[7]} {[8]s:[9]} {[10]s:4}\n";

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

        "TIME",
    };

    try printf(io, fmt, tuple);
}

/// Prints the current iteration step's physical observables and elapsed time.
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

    try printf(io, "{f}\n", .{timer.untilNow(io, .real)});

    timer.* = std.Io.Timestamp.now(io, .real);
}

/// Propagates the classical equations of motion over the specified number of time steps.
fn solve(comptime T: type, io: std.Io, ctx: SolveContext(T), gpa: Allocator, _: Allocator) !Observables(T) {
    const ndim, const nstate = .{ ctx.sim.elpoten.ndim(), ctx.sim.elpoten.nstate() };

    if (ctx.log) try printHeader(io, ndim, nstate);

    var hist = try History(T).init(ndim, nstate, ctx.opt.iterations + 1, ctx.opt.write, gpa);
    defer hist.deinit(gpa);

    var timer = std.Io.Timestamp.now(io, .real);

    for (0..ctx.opt.iterations + 1) |i| {
        const time = @as(T, @floatFromInt(i)) * ctx.opt.time_step;

        if (i > 0) {
            try ctx.sim.propag.step(&ctx.sim.ensemble, &ctx.sim.gb, ctx.sim.elpoten, time);
        }

        const is_log_step = ctx.log and ((i % ctx.opt.log_interval == 0) or (i == ctx.opt.iterations));

        var obs = try Observables(T).init(ctx.sim.*, time, ctx.opt.write, is_log_step, gpa);
        defer obs.deinit(gpa);

        hist.append(obs);

        if (is_log_step) {
            try printIteration(T, io, obs, i, &timer);
        }
    }

    try hist.exportWrite(io, ctx.opt.time_step, ctx.opt.write);

    const end_time = @as(T, @floatFromInt(ctx.opt.iterations)) * ctx.opt.time_step;

    return try Observables(T).init(ctx.sim.*, end_time, ctx.opt.write, true, gpa);
}
