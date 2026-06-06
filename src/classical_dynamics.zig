const std = @import("std");

const Allocator = std.mem.Allocator;

// zig fmt: off
const Matrix           = @import("tensor.zig"   ).    Matrix;
const Potential        = @import("potential.zig"). Potential;
const PotentialOptions = @import("potential.zig").   Options;
const ScalarDual       = @import("dual.zig"     ).ScalarDual;
const Vector           = @import("tensor.zig"   )    .Vector;
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

    state: u32 = 0,
    seed:  u32 = 1,
};
// zig fmt: on

// zig fmt: off
const Write = struct {
    kinetic_energy:   ?[]const u8 = null,
    momentum:         ?[]const u8 = null,
    population:       ?[]const u8 = null,
    position:         ?[]const u8 = null,
    potential_energy: ?[]const u8 = null,
    total_energy:     ?[]const u8 = null,
};
// zig fmt: on

// zig fmt: off
pub const Options = struct {
    initial_conditions: InitialConditions,
    potential:           PotentialOptions,
    time_step:                        f64,
    iterations:                       u32,
    mass:                             f64,
    trajectories:                     u32,

    write:        Write     =   .{},
    adiabatic:    bool      =  true,
    log_interval: u32       =     1,
};
// zig fmt: on

// ENSEMBLE ============================================================================================================

pub fn Ensemble(comptime T: type) type {
    return struct {
        r: Matrix(T),
        p: Matrix(T),
        a: Matrix(T),

        s: Vector(usize),

        // zig fmt: off
        mass:       T,
        nstate: usize,
        // zig fmt: on

        pub fn init(ndim: usize, nstate: usize, ntraj: usize, mass: T, gpa: Allocator) !@This() {
            return .{
                .r = try Matrix(T).init(ntraj, ndim, gpa),
                .p = try Matrix(T).init(ntraj, ndim, gpa),
                .a = try Matrix(T).init(ntraj, ndim, gpa),

                .s = try Vector(usize).init(ntraj, gpa),

                // zig fmt: off
                .mass   =   mass,
                .nstate = nstate,
                // zig fmt: on
            };
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.r.deinit(gpa);
            self.p.deinit(gpa);
            self.a.deinit(gpa);
            self.s.deinit(gpa);
        }

        pub fn pos(self: @This(), gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(self.r.ncol(), gpa);

            for (0..self.r.nrow()) |i| for (0..self.r.ncol()) |j| {
                value.ptr(j).* += self.r.at(i, j);
            };

            value.divs(@floatFromInt(self.r.nrow()));

            return value;
        }

        pub fn mom(self: @This(), gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(self.p.ncol(), gpa);

            for (0..self.p.nrow()) |i| for (0..self.p.ncol()) |j| {
                value.ptr(j).* += self.p.at(i, j);
            };

            value.divs(@floatFromInt(self.p.nrow()));

            return value;
        }

        pub fn pop(self: @This(), gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(self.nstate, gpa);

            for (0..self.s.length()) |i| {
                value.ptr(self.s.at(i)).* += 1.0;
            }

            value.divs(@floatFromInt(self.s.length()));

            return value;
        }

        pub fn ekin(self: @This()) T {
            var sum: T = 0;

            for (0..self.p.nrow()) |i| {
                var psq: T = 0;

                for (0..self.p.ncol()) |j| {
                    psq += self.p.at(i, j) * self.p.at(i, j);
                }

                sum += psq;
            }

            return sum / (2 * self.mass * @as(T, @floatFromInt(self.p.nrow())));
        }

        pub fn epot(self: @This(), pot: Potential(T), t: T, gpa: Allocator) !T {
            var sum: T = 0;

            const V_slice = try gpa.alloc(T, pot.nstate() * pot.nstate());

            for (0..self.r.nrow()) |i| {
                pot.eval(V_slice, self.r.rowSlice(i), t);

                sum += V_slice[self.s.at(i) * self.nstate + self.s.at(i)];
            }

            gpa.free(V_slice);

            return sum / @as(T, @floatFromInt(self.r.nrow()));
        }

        pub fn setGaussian(self: *@This(), ic: InitialConditions) void {
            var split_mix = std.Random.SplitMix64.init(ic.seed);

            for (0..self.s.length()) |i| {
                self.s.ptr(i).* = ic.state;
            }

            var rng = std.Random.DefaultPrng.init(split_mix.next());

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

// PROPAGATOR ==========================================================================================================

pub fn Propagator(comptime T: type) type {
    return struct {
        r_dual: Matrix(ScalarDual(T)),
        V_dual: Matrix(ScalarDual(T)),

        dt: T,

        pub fn init(ndim: usize, nstate: usize, ntraj: usize, dt: T, gpa: Allocator) !@This() {
            return .{
                // zig fmt: off
                .r_dual = try Matrix(ScalarDual(T)).init(ntraj, ndim,            gpa),
                .V_dual = try Matrix(ScalarDual(T)).init(ntraj, nstate * nstate, gpa),
                // zig fmt: on

                .dt = dt,
            };
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.r_dual.deinit(gpa);
            self.V_dual.deinit(gpa);
        }

        pub fn calcAcc(self: *@This(), ensemble: *Ensemble(T), pot: Potential(T), time: T) void {
            for (0..ensemble.r.ncol()) |j| {
                for (0..ensemble.r.nrow()) |i| for (0..ensemble.r.ncol()) |k| {
                    self.r_dual.ptr(i, k).* = ScalarDual(T).init(ensemble.r.ptr(i, k).*, if (k == j) 1 else 0);
                };

                pot.evalDualBatch(&self.V_dual, self.r_dual, time);

                for (0..ensemble.r.nrow()) |i| {
                    const k = ensemble.s.at(i) * pot.nstate() + ensemble.s.at(i);

                    ensemble.a.ptr(i, j).* = -self.V_dual.at(i, k).der / ensemble.mass;
                }
            }
        }

        pub fn step(self: *@This(), ensemble: *Ensemble(T), pot: Potential(T), time: T) void {
            for (0..ensemble.r.nrow()) |i| for (0..ensemble.r.ncol()) |j| {
                ensemble.p.ptr(i, j).* += 0.5 * ensemble.mass * ensemble.a.at(i, j) * self.dt;

                ensemble.r.ptr(i, j).* += (ensemble.p.at(i, j) / ensemble.mass) * self.dt;
            };

            self.calcAcc(ensemble, pot, time);

            for (0..ensemble.r.nrow()) |i| for (0..ensemble.r.ncol()) |j| {
                ensemble.p.ptr(i, j).* += 0.5 * ensemble.mass * ensemble.a.at(i, j) * self.dt;
            };
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

        pub fn init(csys: ClassicalSystem(T), t: T, write: Write, log: bool, gpa: Allocator) !@This() {
            var obs = @This(){};

            // zig fmt: off
            var calc = .{
                .pos  = log or write.position         != null,
                .epot = log or write.potential_energy != null,
                .pop  = log or write.population       != null,
                .mom  = log or write.momentum         != null,
                .ekin = log or write.kinetic_energy   != null,
            };
            // zig fmt: on

            calc.ekin = calc.ekin or write.total_energy != null;
            calc.epot = calc.epot or write.total_energy != null;

            if (calc.mom) obs.mom = try csys.ensemble.mom(gpa);
            if (calc.pos) obs.pos = try csys.ensemble.pos(gpa);
            if (calc.pop) obs.pop = try csys.ensemble.pop(gpa);

            // zig fmt: off
            if (calc.ekin) obs.ekin =     csys.ensemble.ekin(                );
            if (calc.epot) obs.epot = try csys.ensemble.epot(csys.pot, t, gpa);
            // zig fmt: on

            return obs;
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            if (self.pos) |*pos| pos.deinit(gpa);
            if (self.mom) |*mom| mom.deinit(gpa);
            if (self.pop) |*pop| pop.deinit(gpa);
        }
    };
}

// LOGGERS =============================================================================================================

fn printHeader(io: std.Io, ndim: usize, nstate: usize) !void {
    const fmt = "\n{[0]s:8} {[1]s:12} {[2]s:12} {[3]s:12} {[4]s:[5]} {[6]s:[7]} {[8]s:[9]} {[10]s:4}\n";

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

        "TIME",
    };

    try printf(io, fmt, tuple);
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

    try printf(io, "{f}\n", .{timer.untilNow(io, .real)});

    timer.* = std.Io.Timestamp.now(io, .real);
}

// CLASSICAL SYSTEM ====================================================================================================

fn ClassicalSystem(comptime T: type) type {
    return struct { ensemble: Ensemble(T), pot: Potential(T) };
}

// RUN =================================================================================================================

fn SimulationState(comptime T: type) type {
    return struct { csys: ClassicalSystem(T), prop: Propagator(T) };
}

fn SolveContext(comptime T: type) type {
    return struct { opt: Options, sim: *SimulationState(T), log: bool };
}

fn init(comptime T: type, opt: Options, arena: Allocator) !SimulationState(T) {
    const pot = Potential(T).init(opt.potential);

    // zig fmt: off
    var ensemble = try Ensemble  (T).init(pot.ndim(), pot.nstate(), opt.trajectories, opt.mass,      arena);
    var prop     = try Propagator(T).init(pot.ndim(), pot.nstate(), opt.trajectories, opt.time_step, arena);
    // zig fmt: on

    ensemble.setGaussian(opt.initial_conditions);

    prop.calcAcc(&ensemble, pot, 0);

    return .{ .csys = .{ .ensemble = ensemble, .pot = pot }, .prop = prop };
}

fn solve(comptime T: type, io: std.Io, ctx: SolveContext(T), gpa: Allocator, arena: Allocator) !void {
    if (ctx.log) try printHeader(io, ctx.sim.csys.pot.ndim(), ctx.sim.csys.pot.nstate());

    var timer = std.Io.Timestamp.now(io, .real);

    for (0..ctx.opt.iterations + 1) |i| {
        const time = (@as(T, @floatFromInt(i)) - 0.5) * ctx.opt.time_step;

        if (i > 0) ctx.sim.prop.step(&ctx.sim.csys.ensemble, ctx.sim.csys.pot, time);

        const is_log_step = ctx.log and ((i % ctx.opt.log_interval == 0) or (i == ctx.opt.iterations));

        var obs = try Observables(T).init(ctx.sim.csys, time + 0.5 * ctx.opt.time_step, ctx.opt.write, is_log_step, gpa);

        if (is_log_step) {
            try printIteration(T, io, obs, i, &timer);
        }

        obs.deinit(gpa);
    }

    _ = arena;
}

pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator, arena: Allocator) !void {
    if (log) try std.Io.File.stdout().writeStreamingAll(io, "\nCLASSICAL DYNAMICS INIT: ");

    var timer = std.Io.Timestamp.now(io, .real);

    var sim = try init(T, opt, arena);

    if (log) try printf(io, "{f}\n", .{timer.untilNow(io, .real)});

    try solve(T, io, .{ .opt = opt, .sim = &sim, .log = log }, gpa, arena);
}
