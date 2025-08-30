//! Target for performing quantum dynamics simulations.

const std = @import("std");

const device_write = @import("device_write.zig");
const electronic_potential = @import("electronic_potential.zig");
const grid_generator = @import("grid_generator.zig");
const real_matrix = @import("real_matrix.zig");
const complex_vector = @import("complex_vector.zig");
const grid_wavefunction = @import("grid_wavefunction.zig");
const complex_matrix = @import("complex_matrix.zig");
const real_vector = @import("real_vector.zig");
const grid_wavefunction_propagator = @import("grid_wavefunction_propagator.zig");

const ElectronicPotential = electronic_potential.ElectronicPotential;
const GridWavefunctionPropagator = grid_wavefunction_propagator.GridWavefunctionPropagator;
const GridWavefunction = grid_wavefunction.GridWavefunction;
const RealMatrix = real_matrix.RealMatrix;
const ComplexVector = complex_vector.ComplexVector;
const ComplexMatrix = complex_matrix.ComplexMatrix;
const RealVector = real_vector.RealVector;

const momentumGridAlloc = grid_generator.momentumGridAlloc;
const positionGridAlloc = grid_generator.positionGridAlloc;
const printRealMatrix = device_write.printRealMatrix;
const print = device_write.print;

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

        adiabatic: bool = true,
        imaginary: bool = false,
    };
}

/// The quantum dynamics output struct.
pub fn Output(comptime _: type) type {
    return struct {

        /// Allocate the output structure.
        pub fn init(allocator: std.mem.Allocator) !@This() {
            _ = allocator;

            return @This(){
            };
        }

        /// Deallocate the output structure.
        pub fn deinit(self: @This()) void {
            _ = self;
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

/// Run quantum dynamics simulation.
pub fn run(comptime T: type, options: Options(T), allocator: std.mem.Allocator) !Output(T) {
    try print("\nRUNNING QUANTUM DYNAMICS SIMULATION ON A {d}-DIMENSIONAL GRID WITH {d} POINTS\n\n", .{options.grid.limits.len, options.grid.points});

    const ndim = options.potential.ndim();
    const nstate = options.potential.nstate();

    var output = try Output(T).init(allocator); defer output.deinit();

    const position_grid = try positionGridAlloc(T, options.grid.limits, @intCast(options.grid.points), allocator); defer position_grid.deinit();
    const momentum_grid = try momentumGridAlloc(T, options.grid.limits, @intCast(options.grid.points), allocator); defer momentum_grid.deinit();

    var wavefunction = try GridWavefunction(T).init(@intCast(options.grid.points), nstate, ndim, options.initial_conditions.mass, allocator); defer wavefunction.deinit();

    wavefunction.attachGrids(&position_grid, &momentum_grid);

    wavefunction.initialGaussian(options.initial_conditions.position, options.initial_conditions.momentum, options.initial_conditions.state, options.initial_conditions.gamma);

    var propagator = try GridWavefunctionPropagator(T).init(@intCast(options.grid.points), ndim, nstate, allocator); defer propagator.deinit();

    try propagator.generate(wavefunction, options.potential, 0, options.time_step, options.imaginary);

    var temporary_wavefunction_column = try ComplexVector(T).init(wavefunction.data.rows, allocator); defer temporary_wavefunction_column.deinit();

    try printIterationHeader(ndim, nstate);

    for (0..options.iterations + 1) |i| {

        const time = @as(T, @floatFromInt(i)) * options.time_step;

        if (i > 0) try wavefunction.propagate(propagator, options.imaginary, &temporary_wavefunction_column);

        const density_matrix = try wavefunction.density(); defer density_matrix.deinit();

        const potential_energy = try wavefunction.potentialEnergy(options.potential, time);
        const kinetic_energy = try wavefunction.kineticEnergy(&temporary_wavefunction_column);

        const position = try wavefunction.positionMean(); defer position.deinit();
        const momentum = try wavefunction.momentumMean(&temporary_wavefunction_column); defer momentum.deinit();

        if (i > 0 and i % options.log_intervals.iteration != 0) continue;

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
