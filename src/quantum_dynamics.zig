//! Target for performing quantum dynamics simulations.

const std = @import("std");

const complex_matrix = @import("complex_matrix.zig");
const complex_vector = @import("complex_vector.zig");
const device_write = @import("device_write.zig");
const electronic_potential = @import("electronic_potential.zig");
const global_variables = @import("global_variables.zig");
const grid_generator = @import("grid_generator.zig");
const grid_wavefunction = @import("grid_wavefunction.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");
const tully_potential = @import("tully_potential.zig");

const ComplexMatrix = complex_matrix.ComplexMatrix;
const ComplexVector = complex_vector.ComplexVector;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const GridWavefunction = grid_wavefunction.GridWavefunction;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;
const TullyPotential1 = tully_potential.TullyPotential1;

const momentumGridAlloc = grid_generator.momentumGridAlloc;
const positionGridAlloc = grid_generator.positionGridAlloc;
const print = device_write.print;

const TEST_TOLERANCE = global_variables.TEST_TOLERANCE;

/// The quantum dynamics options struct.
pub fn Options(comptime T: type) type {
    return struct {
        pub const Grid = struct {
            limits: []const []const T,
            points: u32
        };
        pub const InitialConditions = struct {
            adiabatic: bool = false,
            gamma: T = 2,
            mass: T,
            momentum: []const T,
            position: []const T,
            state: u32
        };
        pub const LogIntervals = struct {
            iteration: u32 = 1
        };

        potential: ElectronicPotential(T),
        grid: Grid,
        initial_conditions: InitialConditions,

        iterations: u32,
        time_step: T,

        log_intervals: LogIntervals = .{},

        adiabatic: bool = false,
        imaginary: bool = false,
    };
}

/// The quantum dynamics output struct.
pub fn Output(comptime T: type) type {
    return struct {
        population: RealMatrix(f64),

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

/// Container for custom structs related to the quantum dynamics.
pub fn Custom(comptime T: type) type {
    return struct {

        /// Structure to hold information about each iteration used for logging.
        pub const IterationInfo = struct {
            density_matrix: ComplexMatrix(T),
            iteration: usize,
            kinetic_energy: T,
            momentum: RealVector(T),
            position: RealVector(T),
            potential_energy: T,
            time: T
        };
    };
}

const printRealMatrix = device_write.printRealMatrix;

/// Run quantum dynamics simulation.
pub fn run(comptime T: type, options: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try print("\nRUNNING QUANTUM DYNAMICS SIMULATION ON A {d}-DIMENSIONAL GRID WITH {d} POINTS\n\n", .{options.grid.limits.len, options.grid.points});

    const ndim = options.potential.ndim();
    const nstate = options.potential.nstate();

    if (options.grid.limits.len != ndim) return error.IncompatibleDimension;
    if (options.initial_conditions.momentum.len != ndim) return error.IncompatibleDimension;
    if (options.initial_conditions.position.len != ndim) return error.IncompatibleDimension;

    var output = try Output(T).init(nstate, @intCast(options.iterations), allocator);

    const position_grid = try positionGridAlloc(T, options.grid.limits, @intCast(options.grid.points), allocator); defer position_grid.deinit();
    const momentum_grid = try momentumGridAlloc(T, options.grid.limits, @intCast(options.grid.points), allocator); defer momentum_grid.deinit();

    var wavefunction = try GridWavefunction(T).init(@intCast(options.grid.points), nstate, ndim, options.initial_conditions.mass, allocator); defer wavefunction.deinit();

    wavefunction.attachGrids(&position_grid, &momentum_grid);

    wavefunction.initialGaussian(options.initial_conditions.position, options.initial_conditions.momentum, options.initial_conditions.state, options.initial_conditions.gamma);

    if (options.adiabatic) try wavefunction.transformRepresentation(options.potential, 0, false);

    var temporary_wavefunction_column = try ComplexVector(T).init(wavefunction.data.rows, allocator); defer temporary_wavefunction_column.deinit();

    if (enable_printing) try printIterationHeader(ndim, nstate);

    for (0..options.iterations + 1) |i| {

        const time = @as(T, @floatFromInt(i)) * options.time_step;

        if (i > 0) try wavefunction.propagate(options.potential, time, options.time_step, options.imaginary, &temporary_wavefunction_column);

        const density_matrix = try wavefunction.density(options.potential, time, options.adiabatic); defer density_matrix.deinit();

        const potential_energy = try wavefunction.potentialEnergy(options.potential, time);
        const kinetic_energy = try wavefunction.kineticEnergy(&temporary_wavefunction_column);

        const position = try wavefunction.positionMean(); defer position.deinit();
        const momentum = try wavefunction.momentumMean(&temporary_wavefunction_column); defer momentum.deinit();

        for (0..nstate) |j| {
            output.population.ptr(i, j).* = density_matrix.at(j, j).re;
        }

        if (!enable_printing or (i > 0 and i % options.log_intervals.iteration != 0)) continue;

        const iteration_info = Custom(T).IterationInfo{
            .density_matrix = density_matrix,
            .iteration = i,
            .kinetic_energy = kinetic_energy,
            .momentum = momentum,
            .position = position,
            .potential_energy = potential_energy,
            .time = time
        };

        try printIterationInfo(T, iteration_info);
    }

    if (enable_printing) for (0..nstate) |i| {
        try print("{s}FINAL POPULATION OF STATE {d}: {d:.6}\n", .{if (i == 0) "\n" else "", i, output.population.at(options.iterations, i)});
    };

    return output;
}

/// Print header for iteration info.
pub fn printIterationHeader(ndim: usize, nstate: usize) !void {
    var buffer: [128]u8 = undefined;

    var writer = std.io.Writer.fixed(&buffer);

    try writer.print("{s:8} ", .{"ITER"});
    try writer.print("{s:12} {s:12} {s:12} ", .{"KINETIC", "POTENTIAL", "TOTAL"});
    try writer.print("{[value]s:[width]} ", .{.value = "POSITION", .width = 9 * ndim + 2 * (ndim - 1) + 2});
    try writer.print("{[value]s:[width]} ", .{.value = "MOMENTUM", .width = 9 * ndim + 2 * (ndim - 1) + 2});
    try writer.print("{[value]s:[width]}", .{.value = "POPULATION", .width = 9 * nstate + 2 * (nstate - 1) + 2});

    try print("{s}\n", .{writer.buffered()});
}

/// Prints the iteration info to standard output.
pub fn printIterationInfo(comptime T: type, info: Custom(T).IterationInfo) !void {
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

    for (0..info.density_matrix.rows) |i| {
        try writer.print("{d:9.4}{s}", .{info.density_matrix.at(i, i).re, if (i == info.density_matrix.rows - 1) "" else ", "});
    }

    try writer.print("]", .{});

    try print("{s}\n", .{writer.buffered()});
}

test "exact nonadiabatic dynamics on Tully's first potential" {
    const options = Options(f64){
        .potential = .{.tully_1 = TullyPotential1(f64){}},
        .grid = .{
            .limits = &.{&.{-24.0, 32.0}},
            .points = 2048
        },
        .initial_conditions = .{
            .mass = 2000,
            .momentum = &.{15.0},
            .position = &.{-10.0},
            .state = 1,
            .gamma = 2
        },
        .iterations = 350,
        .time_step = 10,
    };

    const output = try run(f64, options, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(@abs(output.population.at(options.iterations, 0) - 0.58949426578088) < TEST_TOLERANCE);
    try std.testing.expect(@abs(output.population.at(options.iterations, 1) - 0.41050573421579) < TEST_TOLERANCE);
}
