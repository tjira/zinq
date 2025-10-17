//! Main file for the vMCG program.

const std = @import("std");

const complex_gaussian = @import("complex_gaussian.zig");
const complex_matrix = @import("complex_matrix.zig");
const complex_vector = @import("complex_vector.zig");
const device_write = @import("device_write.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const electronic_potential = @import("electronic_potential.zig");
const global_variables = @import("global_variables.zig");
const grid_generator = @import("grid_generator.zig");
const grid_wavefunction = @import("grid_wavefunction.zig");
const harmonic_potential = @import("harmonic_potential.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");
const tully_potential = @import("tully_potential.zig");

const Complex = std.math.Complex;
const ComplexGaussian = complex_gaussian.ComplexGaussian;
const ComplexMatrix = complex_matrix.ComplexMatrix;
const ComplexVector = complex_vector.ComplexVector;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const GridWavefunction = grid_wavefunction.GridWavefunction;
const HarmonicPotential = harmonic_potential.HarmonicPotential;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;
const TullyPotential1 = tully_potential.TullyPotential1;

const exportRealMatrix = device_write.exportRealMatrix;
const exportRealMatrixWithLinspacedLeftColumn = device_write.exportRealMatrixWithLinspacedLeftColumn;
const fixGauge = eigenproblem_solver.fixGauge;
const mm = matrix_multiplication.mm;
const momentumGridAlloc = grid_generator.momentumGridAlloc;
const positionAtRow = grid_generator.positionAtRow;
const positionGridAlloc = grid_generator.positionGridAlloc;
const print = device_write.print;
const printJson = device_write.printJson;

/// The vMCG dynamics options struct.
pub fn Options(comptime T: type) type {
    return struct {
        pub const IntegrationGrid = struct {
            limits: []const []const T,
            points: u32
        };
        pub const LogIntervals = struct {
            iteration: u32 = 1
        };
        pub const Write = struct {
            population: ?[]const u8 = null,
            wavefunction: ?[]const u8 = null
        };
        pub const Method = union(enum) {
            vMCG: struct {
                position: []const T,
                gamma: []const T,
                momentum: []const T
            }
        };

        potential: ElectronicPotential(T),
        integration_grid: IntegrationGrid,
        method: Method,

        initial_state: u32,
        iterations: u32,
        time_step: T,
        mass: []const T,

        log_intervals: LogIntervals = .{},
        write: Write = .{},

        adiabatic: bool = false,
    };
}

/// The vMCG dynamics output struct.
pub fn Output(comptime T: type) type {
    return struct {
        kinetic_energy: T = undefined,
        population: RealMatrix(T),
        potential_energy: T = undefined,

        /// Allocate the output structure.
        pub fn init(nstate: usize, iterations: usize, allocator: std.mem.Allocator) !@This() {
            return @This(){
                .population = try RealMatrix(T).init(iterations + 1, nstate, allocator),
            };
        }

        /// Deallocate the output structure.
        pub fn deinit(self: @This()) void {
            self.population.deinit();
        }
    };
}

/// Container for custom structs related to the vMCG.
pub fn Custom(comptime T: type) type {
    return struct {

        /// Structure to hold information about each iteration used for logging.
        pub const IterationInfo = struct {
            iteration: usize,
            kinetic_energy: T,
            coefs: ComplexVector(T),
            momentum: RealVector(T),
            position: RealVector(T),
            potential_energy: T,
            time: T
        };
    };
}

/// Run quantum dynamics simulation.
pub fn run(comptime T: type, options: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(options);

    switch (options.method) {
        .vMCG => return try runSingleGaussian(T, options, enable_printing, allocator)
    }
}

/// Perform the simulation with a single gaussian.
pub fn runSingleGaussian(comptime T: type, options: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(options);

    var potential = options.potential;
    const ndim = potential.ndim();
    const nstate = potential.nstate();

    const file_potential = if (potential == .file) try potential.file.init(allocator) else null; defer if (file_potential) |U| U.deinit();

    var output = try Output(T).init(nstate, @intCast(options.iterations), allocator);

    const gaussian = try complex_gaussian.ComplexGaussian(T).init(options.method.vMCG.position, options.method.vMCG.gamma, options.method.vMCG.momentum, allocator); defer gaussian.deinit(allocator);

    const coefs = try ComplexVector(T).init(nstate, allocator); defer coefs.deinit(); coefs.ptr(options.initial_state).* = Complex(T).init(1, 0);

    var position = try RealVector(T).init(ndim, allocator); defer position.deinit();
    var momentum = try RealVector(T).init(ndim, allocator); defer momentum.deinit();

    try printIterationHeader(ndim, nstate);

    var timer = try std.time.Timer.start();

    for (0..options.iterations + 1) |i| {

        const time = @as(T, @floatFromInt(i)) * options.time_step;

        const V = try gaussian.potential(gaussian, potential, options.integration_grid.limits, options.integration_grid.points, 0, allocator); defer V.deinit();

        const potential_energy = try potentialEnergySingleGaussian(T, V, coefs);
        const kinetic_energy = (try gaussian.kinetic(gaussian, options.mass)).re;

        for (0..ndim) |j| {
            position.ptr(j).* = gaussian.center[j]; momentum.ptr(j).* = gaussian.momentum[j];
        }

        const info = Custom(T).IterationInfo{
            .iteration = i,
            .kinetic_energy = kinetic_energy,
            .coefs = coefs,
            .momentum = momentum,
            .position = position,
            .potential_energy = potential_energy,
            .time = time
        };

        try printIterationInfo(T, info, &timer);
    }

    output = output;

    return output;
}

/// Returns the potential energy of a single complex Gaussian given its matrix element of the potential operator and coefficients.
pub fn potentialEnergySingleGaussian(comptime T: type, V: ComplexMatrix(T), coefs: ComplexVector(T)) !T {
    var potential_energy: Complex(T) = Complex(T).init(0, 0);

    for (0..coefs.len) |i| for (0..coefs.len) |j| {
        potential_energy = potential_energy.add(coefs.at(i).conjugate().mul(V.at(i, j)).mul(coefs.at(j)));
    };

    return potential_energy.re;
}

/// Print header for iteration info.
pub fn printIterationHeader(ndim: usize, nstate: usize) !void {
    var buffer: [128]u8 = undefined;

    var writer = std.io.Writer.fixed(&buffer);

    try writer.print("\n{s:8} ", .{"ITER"});
    try writer.print("{s:12} {s:12} {s:12} ", .{"KINETIC", "POTENTIAL", "TOTAL"});
    try writer.print("{[value]s:[width]} ", .{.value = "POSITION", .width = 9 * ndim + 2 * (ndim - 1) + 2});
    try writer.print("{[value]s:[width]} ", .{.value = "MOMENTUM", .width = 9 * ndim + 2 * (ndim - 1) + 2});
    try writer.print("{[value]s:[width]} ", .{.value = "POPULATION", .width = 9 * nstate + 2 * (nstate - 1) + 2});
    try writer.print("{s:4}", .{"TIME"});

    try print("{s}\n", .{writer.buffered()});
}

/// Prints the iteration info to standard output.
pub fn printIterationInfo(comptime T: type, info: Custom(T).IterationInfo, timer: *std.time.Timer) !void {
    var buffer: [128]u8 = undefined;

    var writer = std.io.Writer.fixed(&buffer);

    try writer.print("{d:8} ", .{info.iteration});
    try writer.print("{d:12.6} {d:12.6} {d:12.6} ", .{info.kinetic_energy, info.potential_energy, info.kinetic_energy + info.potential_energy});

    try writer.print("[", .{});

    for (0..info.position.len) |i| {
        try writer.print("{d:9.4}{s}", .{info.position.at(i), if (i == info.position.len - 1) "" else ", "});
    }

    try writer.print("] [", .{});

    for (0..info.momentum.len) |i| {
        try writer.print("{d:9.4}{s}", .{info.momentum.at(i), if (i == info.momentum.len - 1) "" else ", "});
    }

    try writer.print("] [", .{});

    for (0..info.coefs.len) |i| {
        try writer.print("{d:9.4}{s}", .{info.coefs.at(i).magnitude(), if (i == info.coefs.len - 1) "" else ", "});
    }

    try writer.print("] {D}", .{timer.read()}); timer.reset();

    try print("{s}\n", .{writer.buffered()});
}
