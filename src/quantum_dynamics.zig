const std = @import("std");

const Complex = std.math.Complex;

const Matrix = @import("matrix.zig").Matrix;
const Potential = @import("potential.zig").Potential;
const PotentialOptions = @import("potential.zig").Options;
const Vector = @import("vector.zig").Vector;

const eighMany = @import("openblas.zig").eighMany;
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
    population: ?[]const u8 = null,
    position: ?[]const u8 = null,
    potential_energy: ?[]const u8 = null,
    wavefunction: ?[]const u8 = null,
    total_energy: ?[]const u8 = null,
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
        r: Matrix(T),
        k: Matrix(T),
        dr: T,
        dk: T,

        pub fn init(bounds: []const [2]T, npoint: u32, gpa: std.mem.Allocator) !@This() {
            const nrow = std.math.pow(usize, npoint, bounds.len);

            var r = try Matrix(T).init(nrow, bounds.len, gpa);
            var k = try Matrix(T).init(nrow, bounds.len, gpa);

            var dr: T = 1;
            var dk: T = 1;

            for (0..bounds.len) |i| {
                const min = bounds[i][0];
                const max = bounds[i][1];

                dr *= (max - min) / @as(T, @floatFromInt(npoint));
            }

            for (0..bounds.len) |i| {
                const min = bounds[i][0];
                const max = bounds[i][1];

                dk *= 2 * std.math.pi / (max - min);
            }

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

            return @This(){ .r = r, .k = k, .dr = dr, .dk = dk };
        }

        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            self.r.deinit(gpa);
            self.k.deinit(gpa);
        }
    };
}

// WAVEFUNCTION ========================================================================================================

fn Wavefunction(comptime T: type) type {
    return struct {
        data: Matrix(Complex(T)),
        decay: Vector(T),

        pub fn init(ic: InitialConditions, grid: Grid(T), nstate: usize, gpa: std.mem.Allocator) !@This() {
            var data = try Matrix(Complex(T)).initZero(grid.r.nrow(), nstate, gpa);

            for (0..grid.r.nrow()) |i| {
                var exponent = Complex(T).init(0, 0);

                for (0..grid.r.ncol()) |j| {
                    const dx = grid.r.at(i, j) - ic.position[j];

                    exponent = exponent.add(Complex(T).init(-0.5 * ic.gamma[j] * dx * dx, ic.momentum[j] * dx));
                }

                data.ptr(i, ic.state).* = std.math.complex.exp(exponent);
            }

            var wfn: @This() = .{ .data = data, .decay = try Vector(T).initZero(nstate, gpa) };

            wfn.normalize(grid);

            return wfn;
        }

        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            self.data.deinit(gpa);
            self.decay.deinit(gpa);
        }

        pub fn norm(self: @This(), grid: Grid(T)) T {
            var value: T = 0;

            for (0..self.data.ncol()) |j| for (0..self.data.nrow()) |i| {
                value += self.data.at(i, j).squaredMagnitude();
            };

            return value * grid.dr;
        }

        pub fn normalize(self: *@This(), grid: Grid(T)) void {
            self.data.divs(Complex(T).init(std.math.sqrt(self.norm(grid)), 0));
        }

        pub fn overlap(self: @This(), other: @This(), grid: Grid(T)) Complex(T) {
            var value = Complex(T).init(0, 0);

            for (0..self.data.ncol()) |j| for (0..other.data.ncol()) |k| for (0..self.data.nrow()) |i| {
                value = value.add(self.data.at(i, j).conjugate().mul(other.data.at(i, k)));
            };

            return value.mul(Complex(T).init(grid.dr, 0));
        }

        pub fn pe(self: @This(), ham: Hamiltonian(T), grid: Grid(T)) T {
            var value: T = 0;

            for (0..self.data.ncol()) |j| for (0..self.data.ncol()) |k| for (0..self.data.nrow()) |i| {
                const psi_j = self.data.at(i, j);
                const psi_k = self.data.at(i, k);

                const cross_re = psi_j.re * psi_k.re + psi_j.im * psi_k.im;
                
                value += cross_re * ham.V.at(j * self.data.ncol() + k, i);
            };

            return value * grid.dr;
        }

        pub fn pop(self: @This(), grid: Grid(T), gpa: std.mem.Allocator) !Vector(T) {
            var value = try Vector(T).initZero(self.data.ncol(), gpa);

            for (0..self.data.ncol()) |j| for (0..self.data.ncol()) |k| for (0..self.data.nrow()) |i| {
                const psi_j = self.data.at(i, j);
                const psi_k = self.data.at(i, k);

                value.ptr(j).* += psi_j.re * psi_k.re + psi_j.im * psi_k.im;
            };

            value.muls(grid.dr);

            return value;
        }

        pub fn pos(self: @This(), grid: Grid(T), gpa: std.mem.Allocator) !Vector(T) {
            var value = try Vector(T).initZero(grid.r.ncol(), gpa);

            for (0..grid.r.ncol()) |j| for (0..self.data.ncol()) |k| for (0..self.data.nrow()) |i| {
                value.ptr(j).* += self.data.at(i, k).squaredMagnitude() * grid.r.at(i, j);
            };

            value.muls(grid.dr);

            return value;
        }
    };
}

// HAMILTONIAN =========================================================================================================

pub fn Hamiltonian(comptime T: type) type {
    return struct {
        V: Matrix(T),
        W: Matrix(T),
        U: Matrix(T),
        K: Vector(T),

        pub fn init(grid: Grid(T), pot: Potential(T), m: T, gpa: std.mem.Allocator) !@This() {
            var V = try Matrix(T).init(pot.nstate() * pot.nstate(), grid.r.nrow(), gpa);
            var U = try Matrix(T).init(pot.nstate() * pot.nstate(), grid.r.nrow(), gpa);

            var W = try Matrix(T).init(pot.nstate(), grid.r.nrow(), gpa);

            var K = try Vector(T).initZero(grid.r.nrow(), gpa);

            pot.evalMany(&V, grid.r);

            for (0..grid.r.ncol()) |j| for (0..grid.r.nrow()) |i| {
                K.ptr(i).* += 0.5 * grid.k.at(i, j) * grid.k.at(i, j) / m;
            };

            try eighMany(T, &W, &U, V);

            return .{ .V = V, .W = W, .U = U, .K = K };
        }

        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            self.V.deinit(gpa);
            self.W.deinit(gpa);
            self.U.deinit(gpa);
            self.K.deinit(gpa);
        }
    };
}

// PROPAGATOR ==========================================================================================================

pub fn Propagator(comptime T: type) type {
    return struct {
        R: Matrix(Complex(T)),
        K: Vector(Complex(T)),

        pub fn init(ham: Hamiltonian(T), dt: T, gpa: std.mem.Allocator) !@This() {
            var R = try Matrix(Complex(T)).init(ham.V.nrow(), ham.V.ncol(), gpa);

            const nstate = std.math.sqrt(ham.V.nrow());

            for (0..ham.V.ncol()) |i| for (0..nstate) |j| for (0..nstate) |k| {
                var sum = Complex(T).init(0, 0);

                for (0..nstate) |m| {
                    const U_jm = ham.U.at(j * nstate + m, i);
                    const U_km = ham.U.at(k * nstate + m, i);

                    const phase = std.math.complex.exp(Complex(T).init(0, -ham.W.at(m, i) * dt));

                    sum = sum.add(Complex(T).init(phase.re * U_jm * U_km, phase.im * U_jm * U_km));
                }

                R.ptr(j * nstate + k, i).* = sum;
            };

            var K = try Vector(Complex(T)).initZero(ham.K.length(), gpa);

            for (0..ham.K.length()) |i| {
                K.ptr(i).* = std.math.complex.exp(Complex(T).init(0, -ham.K.at(i) * dt));
            }

            return .{ .R = R, .K = K };
        }

        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            self.R.deinit(gpa);
            self.K.deinit(gpa);
        }
    };
}

// LOGGERS =============================================================================================================

pub fn printHeader(io: std.Io, ndim: usize, nstate: usize) !void {
    try printf(io, "\n{s:8} ", .{"ITER"});
    try printf(io, "{s:12} {s:12} {s:12} ", .{ "KIN (Eh)", "POT (Eh)", "TOT (Eh)" });
    try printf(io, "{[value]s:[width]} ", .{ .value = "POS (a0)", .width = 11 * ndim + 1 });
    try printf(io, "{[value]s:[width]} ", .{ .value = "MOM (hb/a0)", .width = 11 * ndim + 1 });
    try printf(io, "{[value]s:[width]} ", .{ .value = "POPULATION", .width = 11 * nstate + 1 });
    try printf(io, "{s:10} {s:4}\n", .{ "NORM", "TIME" });
}

// RUN =================================================================================================================

pub fn run(comptime T: type, io: std.Io, opt: Options, gpa: std.mem.Allocator) !void {
    const pot = Potential(T).init(opt.potential);

    var grid = try Grid(T).init(opt.grid.bounds, opt.grid.npoint, gpa);
    defer grid.deinit(gpa);

    var wfn = try Wavefunction(T).init(opt.initial_conditions, grid, pot.nstate(), gpa);
    defer wfn.deinit(gpa);

    var ham = try Hamiltonian(T).init(grid, pot, opt.mass, gpa);
    defer ham.deinit(gpa);

    var prop = try Propagator(T).init(ham, opt.time_step, gpa);
    defer prop.deinit(gpa);
    
    try printHeader(io, grid.r.ncol(), wfn.data.ncol());

    const pe = wfn.pe(ham, grid);

    var pos = try wfn.pos(grid, gpa);
    defer pos.deinit(gpa);

    var pop = try wfn.pop(grid, gpa);
    defer pop.deinit(gpa);

    const norm = wfn.norm(grid);

    try printf(io, "{d:20.14}\n", .{ pe });
    try printf(io, "{any}\n", .{ pos.data });
    try printf(io, "{any}\n", .{ pop.data });
    try printf(io, "{d:20.14}\n", .{ norm });
}
