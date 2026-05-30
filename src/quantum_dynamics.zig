const std = @import("std");

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;

const Matrix = @import("tensor.zig").Matrix;
const Potential = @import("potential.zig").Potential;
const PotentialOptions = @import("potential.zig").Options;
const Vector = @import("tensor.zig").Vector;

const eighMany = @import("openblas.zig").eighMany;
const fftn = @import("fftw.zig").fftn;
const printf = @import("writer.zig").printf;
const writeMatrix = @import("writer.zig").writeMatrix;

// OPTIONS =============================================================================================================

const InitialConditions = struct {
    position: []const f64,
    momentum: []const f64,
    gamma: []const f64,
    state: u32 = 0,
    adiabatic: bool = false,
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
    adiabatic: bool = false,
    grid: struct {
        bounds: []const [2]f64,
        npoint: u32,
    },
    imaginary: ?struct {
        nstate: u32 = 1,
    } = null,
    initial_conditions: InitialConditions,
    iterations: u32,
    log_interval: u32 = 1,
    mass: f64,
    potential: PotentialOptions,
    time_step: f64,
    write: Write = .{},
};

// GRID ================================================================================================================

fn Grid(comptime T: type) type {
    return struct {
        shape: []i32,
        r: Matrix(T),
        k: Matrix(T),
        dr: T,
        dk: T,

        pub fn init(bounds: []const [2]T, npoint: u32, gpa: Allocator) !@This() {
            const nrow = std.math.pow(usize, npoint, bounds.len);

            var r = try Matrix(T).init(nrow, bounds.len, gpa);
            var k = try Matrix(T).init(nrow, bounds.len, gpa);

            var dr: T = 1;

            for (0..bounds.len) |i| {
                const min = bounds[i][0];
                const max = bounds[i][1];

                dr *= (max - min) / @as(T, @floatFromInt(npoint));
            }

            const dk = dr / @as(T, @floatFromInt(nrow));

            for (0..nrow) |i| {
                var j = i;

                for (0..bounds.len) |l| {
                    const m = bounds.len - l - 1;

                    const n = @as(T, @floatFromInt(j % npoint));

                    const min = bounds[m][0];
                    const max = bounds[m][1];

                    r.ptr(i, m).* = min + n * (max - min) / @as(T, @floatFromInt(npoint));

                    const dki = 2 * std.math.pi / (max - min);

                    k.ptr(i, m).* = (if (j % npoint < npoint / 2) n else n - @as(T, @floatFromInt(npoint))) * dki;

                    j /= npoint;
                }
            }

            var shape = try gpa.alloc(i32, bounds.len);

            for (0..bounds.len) |i| {
                shape[i] = @as(i32, @intCast(npoint));
            }

            return .{ .shape = shape, .r = r, .k = k, .dr = dr, .dk = dk };
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            gpa.free(self.shape);

            self.r.deinit(gpa);
            self.k.deinit(gpa);
        }
    };
}

// WAVEFUNCTION ========================================================================================================

fn Wavefunction(comptime T: type) type {
    return struct {
        W: Matrix(Complex(T)),

        pub fn init(ic: InitialConditions, grid: Grid(T), nstate: usize, gpa: Allocator) !@This() {
            var W = try Matrix(Complex(T)).initZero(grid.r.nrow(), nstate, gpa);

            for (0..grid.r.nrow()) |i| {
                var exponent = Complex(T).init(0, 0);

                for (0..grid.r.ncol()) |j| {
                    const dx = grid.r.at(i, j) - ic.position[j];

                    exponent = exponent.add(Complex(T).init(-0.5 * ic.gamma[j] * dx * dx, ic.momentum[j] * dx));
                }

                W.ptr(i, ic.state).* = std.math.complex.exp(exponent);
            }

            var wfn: @This() = .{ .W = W };

            wfn.normalize(grid);

            return wfn;
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.W.deinit(gpa);
        }

        pub fn clone(self: @This(), gpa: Allocator) !@This() {
            return .{ .W = try self.W.clone(gpa) };
        }

        pub fn ekin(self: @This(), ham: Hamiltonian(T), grid: Grid(T)) T {
            var value: T = 0;

            for (0..self.W.ncol()) |j| for (0..self.W.nrow()) |i| {
                value += self.W.at(i, j).squaredMagnitude() * ham.K.at(i);
            };

            return value * grid.dk;
        }

        pub fn epot(self: @This(), ham: Hamiltonian(T), grid: Grid(T)) T {
            var value: T = 0;

            for (0..self.W.ncol()) |j| for (0..self.W.ncol()) |k| for (0..self.W.nrow()) |i| {
                const psi_j = self.W.at(i, j);
                const psi_k = self.W.at(i, k);

                const cross_re = psi_j.re * psi_k.re + psi_j.im * psi_k.im;

                value += cross_re * ham.V.at(j * self.W.ncol() + k, i);
            };

            return value * grid.dr;
        }

        pub fn fft(self: *@This(), grid: Grid(T), sign: i32) !void {
            for (0..self.W.ncol()) |j| {
                try fftn(Complex(T), self.W.colSlice(j), grid.shape, sign);
            }
        }

        pub fn mom(self: @This(), grid: Grid(T), gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(grid.k.ncol(), gpa);

            for (0..grid.r.ncol()) |j| for (0..self.W.ncol()) |k| for (0..self.W.nrow()) |i| {
                value.ptr(j).* += self.W.at(i, k).squaredMagnitude() * grid.k.at(i, j);
            };

            value.muls(grid.dk);

            return value;
        }

        pub fn norm(self: @This(), grid: Grid(T)) T {
            var value: T = 0;

            for (0..self.W.ncol()) |j| for (0..self.W.nrow()) |i| {
                value += self.W.at(i, j).squaredMagnitude();
            };

            return value * grid.dr;
        }

        pub fn normalize(self: *@This(), grid: Grid(T)) void {
            self.W.divs(Complex(T).init(std.math.sqrt(self.norm(grid)), 0));
        }

        pub fn overlap(self: @This(), other: @This(), grid: Grid(T)) Complex(T) {
            var value = Complex(T).init(0, 0);

            for (0..self.W.ncol()) |j| for (0..other.data.ncol()) |k| for (0..self.W.nrow()) |i| {
                value = value.add(self.W.at(i, j).conjugate().mul(other.data.at(i, k)));
            };

            return value.mul(Complex(T).init(grid.dr, 0));
        }

        pub fn pop(self: @This(), grid: Grid(T), gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(self.W.ncol(), gpa);

            for (0..self.W.ncol()) |j| for (0..self.W.nrow()) |i| {
                value.ptr(j).* += self.W.at(i, j).squaredMagnitude();
            };

            value.muls(grid.dr);

            return value;
        }

        pub fn pos(self: @This(), grid: Grid(T), gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(grid.r.ncol(), gpa);

            for (0..grid.r.ncol()) |j| for (0..self.W.ncol()) |k| for (0..self.W.nrow()) |i| {
                value.ptr(j).* += self.W.at(i, k).squaredMagnitude() * grid.r.at(i, j);
            };

            value.muls(grid.dr);

            return value;
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
            var V = try Matrix(T).init(pot.nstate() * pot.nstate(), grid.r.nrow(), gpa);
            var U = try Matrix(T).init(pot.nstate() * pot.nstate(), grid.r.nrow(), gpa);

            var W = try Matrix(T).init(pot.nstate(), grid.r.nrow(), gpa);

            var K = try Vector(T).initZero(grid.r.nrow(), gpa);

            pot.evalBatch(&V, grid.r);

            for (0..grid.r.ncol()) |j| for (0..grid.r.nrow()) |i| {
                K.ptr(i).* += 0.5 * grid.k.at(i, j) * grid.k.at(i, j) / m;
            };

            try eighMany(T, &W, &U, V);

            return .{ .V = V, .W = W, .U = U, .K = K };
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.V.deinit(gpa);
            self.W.deinit(gpa);
            self.U.deinit(gpa);
            self.K.deinit(gpa);
        }
    };
}

// QUANTUM SYSTEM ======================================================================================================

pub fn QuantumSystem(comptime T: type) type {
    return struct { grid: Grid(T), ham: Hamiltonian(T), wfn: Wavefunction(T) };
}

// PROPAGATOR ==========================================================================================================

fn Propagator(comptime T: type) type {
    return struct {
        R: Matrix(Complex(T)),
        K: Vector(Complex(T)),

        temp: []Complex(T),

        pub fn init(ham: Hamiltonian(T), dt: T, gpa: Allocator) !@This() {
            var R = try Matrix(Complex(T)).init(ham.V.nrow(), ham.V.ncol(), gpa);

            const nstate = std.math.sqrt(ham.V.nrow());

            for (0..R.ncol()) |i| for (0..nstate) |j| for (0..nstate) |k| {
                var sum = Complex(T).init(0, 0);

                for (0..nstate) |m| {
                    const U_jm = ham.U.at(j * nstate + m, i);
                    const U_km = ham.U.at(k * nstate + m, i);

                    const phase = std.math.complex.exp(Complex(T).init(0, -0.5 * ham.W.at(m, i) * dt));

                    sum = sum.add(Complex(T).init(phase.re * U_jm * U_km, phase.im * U_jm * U_km));
                }

                R.ptr(j * nstate + k, i).* = sum;
            };

            var K = try Vector(Complex(T)).initZero(ham.K.length(), gpa);

            for (0..ham.K.length()) |i| {
                K.ptr(i).* = std.math.complex.exp(Complex(T).init(0, -ham.K.at(i) * dt));
            }

            const temp = try gpa.alloc(Complex(T), nstate);

            return .{ .R = R, .K = K, .temp = temp };
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            gpa.free(self.temp);

            self.R.deinit(gpa);
            self.K.deinit(gpa);
        }

        pub fn step(self: @This(), wfn: *Wavefunction(T), grid: Grid(T)) !void {
            self.applyR(wfn);
            try self.applyK(wfn, grid);
            self.applyR(wfn);
        }

        fn applyK(self: @This(), wfn: *Wavefunction(T), grid: Grid(T)) !void {
            try wfn.fft(grid, -1);

            for (0..wfn.W.ncol()) |j| for (0..wfn.W.nrow()) |i| {
                wfn.W.ptr(i, j).* = self.K.at(i).mul(wfn.W.at(i, j));
            };

            try wfn.fft(grid, 1);
        }

        fn applyR(self: @This(), wfn: *Wavefunction(T)) void {
            const nstate = wfn.W.ncol();

            for (0..wfn.W.nrow()) |i| {
                for (0..nstate) |j| {
                    var sum = Complex(T).init(0, 0);

                    for (0..nstate) |k| {
                        sum = sum.add(self.R.at(j * nstate + k, i).mul(wfn.W.at(i, k)));
                    }

                    self.temp[j] = sum;
                }

                for (0..nstate) |j| {
                    wfn.W.ptr(i, j).* = self.temp[j];
                }
            }
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

        pub fn init(qsys: *QuantumSystem(T), write: Write, log: bool, gpa: Allocator) !@This() {
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
            if (calc.pop ) obs.pop  = try qsys.wfn.pop (          qsys.grid, gpa);
            // zig fmt: on

            const needs_fft = calc.mom or calc.ekin;

            if (needs_fft) {
                try qsys.wfn.fft(qsys.grid, -1);

                // zig fmt: off
                if (calc.mom ) obs.mom  = try qsys.wfn.mom (          qsys.grid, gpa);
                if (calc.ekin) obs.ekin =     qsys.wfn.epot(qsys.ham, qsys.grid     );
                // zig fmt: on

                try qsys.wfn.fft(qsys.grid, 1);
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

// LOGGERS =============================================================================================================

fn printHeader(io: std.Io, ndim: usize, nstate: usize) !void {
    const fmt = "\n{[0]s:8} {[1]s:12} {[2]s:12} {[3]s:12} {[4]s:[5]} {[6]s:[7]} {[8]s:[9]} {[10]s:10} {[11]s:4}\n";

    const tuple = .{
        "ITER",

        "EKIN (Eh)",
        "EPOT (Eh)",
        "ETOT (Eh)",

        // zig fmt: off
        "POS (a0)",    11 * ndim   + 1,
        "MOM (hb/a0)", 11 * ndim   + 1,
        "POP (-)",     10 * nstate + 2,
        // zig fmt: on

        "NORM (-)",
        "TIME",
    };

    try printf(io, fmt, tuple);
}

pub fn printIteration(comptime T: type, io: std.Io, obs: Observables(T), i: usize, timer: *std.Io.Timestamp) !void {
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

pub fn SimulationState(comptime T: type) type {
    return struct { qsys: QuantumSystem(T), prop: Propagator(T), pot: Potential(T) };
}

fn init(comptime T: type, opt: Options, arena: Allocator) !SimulationState(T) {
    const pot = Potential(T).init(opt.potential);

    const grid = try Grid(T).init(opt.grid.bounds, opt.grid.npoint, arena);
    const wfn = try Wavefunction(T).init(opt.initial_conditions, grid, pot.nstate(), arena);
    const ham = try Hamiltonian(T).init(grid, pot, opt.mass, arena);
    const prop = try Propagator(T).init(ham, opt.time_step, arena);

    return .{ .qsys = .{ .grid = grid, .ham = ham, .wfn = wfn }, .prop = prop, .pot = pot };
}

pub fn run(comptime T: type, io: std.Io, opt: Options, gpa: Allocator, arena: Allocator) !void {
    var sim = try init(T, opt, arena);

    try printHeader(io, sim.qsys.grid.r.ncol(), sim.qsys.wfn.W.ncol());

    var timer = std.Io.Timestamp.now(io, .real);

    for (0..opt.iterations + 1) |i| {
        // const time = (@as(T, @floatFromInt(i)) + 0.5) * opt.time_step;

        if (i > 0) try sim.prop.step(&sim.qsys.wfn, sim.qsys.grid);

        const is_log_step = (i % opt.log_interval == 0) or (i == opt.iterations);

        var obs = try Observables(T).init(&sim.qsys, opt.write, is_log_step, gpa);

        if (is_log_step) {
            try printIteration(T, io, obs, i, &timer);
        }

        obs.deinit(gpa);
    }
}
