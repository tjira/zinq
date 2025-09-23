//! Target for performing quantum dynamics simulations.

const std = @import("std");

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
        pub const Write = struct {
            population: ?[]const u8 = null,
            wavefunction: ?[]const u8 = null
        };

        potential: ElectronicPotential(T),
        grid: Grid,
        initial_conditions: InitialConditions,

        iterations: u32,
        time_step: T,

        log_intervals: LogIntervals = .{},
        write: Write = .{},

        adiabatic: bool = false,
        imaginary: bool = false,
    };
}

/// The quantum dynamics output struct.
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
pub fn run(comptime T: type, options: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(options);

    var potential = options.potential;
    const ndim = potential.ndim();
    const nstate = potential.nstate();

    const file_potential = if (potential == .file) try potential.file.init(allocator) else null; defer if (file_potential) |U| U.deinit();

    var output = try Output(T).init(nstate, @intCast(options.iterations), allocator);

    var wavefunction = try GridWavefunction(T).init(@intCast(options.grid.points), nstate, ndim, options.grid.limits, options.initial_conditions.mass, allocator); defer wavefunction.deinit();

    try wavefunction.initialGaussian(options.initial_conditions.position, options.initial_conditions.momentum, options.initial_conditions.state, options.initial_conditions.gamma);

    if (options.adiabatic) try wavefunction.transformRepresentation(potential, 0, false);

    var wavefunction_dynamics: ?RealMatrix(T) = if (options.write.wavefunction) |_| try initializeWavefunctionDynamicsContainer(T, wavefunction, options.iterations, allocator) else null;

    var temporary_wavefunction_column = try ComplexVector(T).init(wavefunction.data.rows, allocator); defer temporary_wavefunction_column.deinit();

    if (enable_printing) try printIterationHeader(ndim, nstate);

    var timer = try std.time.Timer.start();

    for (0..options.iterations + 1) |i| {

        const time = @as(T, @floatFromInt(i)) * options.time_step;

        if (i > 0) try wavefunction.propagate(potential, time, options.time_step, options.imaginary, &temporary_wavefunction_column);

        const density_matrix = try wavefunction.density(potential, time, options.adiabatic); defer density_matrix.deinit();

        const potential_energy = try wavefunction.potentialEnergy(potential, time);
        const kinetic_energy = try wavefunction.kineticEnergy(&temporary_wavefunction_column);

        const position = try wavefunction.positionMean(); defer position.deinit();
        const momentum = try wavefunction.momentumMean(&temporary_wavefunction_column); defer momentum.deinit();

        for (0..nstate) |j| {
            output.population.ptr(i, j).* = density_matrix.at(j, j).re;
        }

        if (options.write.wavefunction != null) try assignWavefunctionStep(T, &wavefunction_dynamics.?, wavefunction, potential, i, time, options.adiabatic);

        if (i == options.iterations) {
            output.kinetic_energy = kinetic_energy;
            output.potential_energy = potential_energy;
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

        try printIterationInfo(T, iteration_info, &timer);
    }

    if (enable_printing) for (0..nstate) |i| {
        try print("{s}FINAL POPULATION OF STATE {d}: {d:.6}\n", .{if (i == 0) "\n" else "", i, output.population.at(options.iterations, i)});
    };

    const end_time = @as(T, @floatFromInt(options.iterations)) * options.time_step;

    if (options.write.population) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.population, 0, end_time, options.iterations + 1);
    if (options.write.wavefunction) |path| try exportRealMatrix(T, path, wavefunction_dynamics.?);

    if (wavefunction_dynamics != null) wavefunction_dynamics.?.deinit();

    return output;
}

/// Assign current wavefunction to the wavefunction dynamics matrix.
pub fn assignWavefunctionStep(comptime T: type, wavefunction_dynamics: *RealMatrix(T), wavefunction: GridWavefunction(T), potential: ElectronicPotential(T), iter: usize, time: T, adiabatic: bool) !void {
    var diabatic_potential = try RealMatrix(T).init(wavefunction.nstate, wavefunction.nstate, wavefunction.allocator); defer diabatic_potential.deinit();
    var adiabatic_potential = try RealMatrix(T).init(wavefunction.nstate, wavefunction.nstate, wavefunction.allocator); defer adiabatic_potential.deinit();
    var adiabatic_eigenvectors = try RealMatrix(T).init(wavefunction.nstate, wavefunction.nstate, wavefunction.allocator); defer adiabatic_eigenvectors.deinit();
    var previous_eigenvectors = try RealMatrix(T).init(wavefunction.nstate, wavefunction.nstate, wavefunction.allocator); defer previous_eigenvectors.deinit();

    var position_at_row = try RealVector(T).init(wavefunction.ndim, wavefunction.allocator); defer position_at_row.deinit();

    var wavefunction_row = try ComplexMatrix(T).init(wavefunction.nstate, 1, wavefunction.allocator); defer wavefunction_row.deinit();

    for (0..wavefunction.data.rows) |i| {

        for (0..wavefunction.nstate) |j| wavefunction_row.ptr(j, 0).* = wavefunction.data.at(i, j);

        if (adiabatic) {

            try adiabatic_eigenvectors.copyTo(&previous_eigenvectors);

            positionAtRow(T, &position_at_row, i, wavefunction.ndim, wavefunction.npoint, wavefunction.limits);

            try potential.evaluateEigensystem(&diabatic_potential, &adiabatic_potential, &adiabatic_eigenvectors, position_at_row, time);

            if (i > 0) try fixGauge(T, &adiabatic_eigenvectors, previous_eigenvectors);

            try mm(T, &wavefunction_row, adiabatic_eigenvectors, true, wavefunction.data.row(i).asMatrix(), false);
        }

        for (0..wavefunction.data.cols) |j| {
            wavefunction_dynamics.ptr(i, wavefunction.ndim + 2 * wavefunction.nstate * iter + 2 * j + 0).* = wavefunction_row.at(j, 0).re;
            wavefunction_dynamics.ptr(i, wavefunction.ndim + 2 * wavefunction.nstate * iter + 2 * j + 1).* = wavefunction_row.at(j, 0).im;
        }
    }
}

/// Initialize the container for the wavefunction dynamics.
pub fn initializeWavefunctionDynamicsContainer(comptime T: type, wavefunction: GridWavefunction(T), iterations: usize, allocator: std.mem.Allocator) !RealMatrix(T) {
    var wavefunction_dynamics: RealMatrix(T) = try RealMatrix(T).init(wavefunction.data.rows, wavefunction.ndim + 2 * wavefunction.nstate * (iterations + 1), allocator);

    var position_at_row = try RealVector(T).init(wavefunction.ndim, allocator); defer position_at_row.deinit();

    for (0..wavefunction.data.rows) |i| {

        positionAtRow(T, &position_at_row, i, wavefunction.ndim, wavefunction.npoint, wavefunction.limits);

        for (0..wavefunction.ndim) |j| {
            wavefunction_dynamics.ptr(i, j).* = position_at_row.at(j);
        }
    }

    return wavefunction_dynamics;
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

    for (0..info.density_matrix.rows) |i| {
        try writer.print("{d:9.4}{s}", .{info.density_matrix.at(i, i).re, if (i == info.density_matrix.rows - 1) "" else ", "});
    }

    try writer.print("] {D}", .{timer.read()}); timer.reset();

    try print("{s}\n", .{writer.buffered()});
}

test "Exact Dynamics on 1D Harmonic Potential" {
    const options = Options(f64){
        .grid = .{
            .limits = &.{&.{-8, 8}},
            .points = 64
        },
        .initial_conditions = .{
            .mass = 1,
            .momentum = &.{0},
            .position = &.{1},
            .state = 0,
            .gamma = 2
        },
        .potential = .{
            .harmonic = HarmonicPotential(f64){}
        },
        .iterations = 1000,
        .time_step = 0.1
    };

    const output = try run(f64, options, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(@abs(output.kinetic_energy - 0.52726330098766) < TEST_TOLERANCE);
    try std.testing.expect(@abs(output.potential_energy - 0.59766836993815) < TEST_TOLERANCE);
}

test "Exact Dynamics on 2D Harmonic Potential" {
    const options = Options(f64){
        .grid = .{
            .limits = &.{&.{-8, 8}, &.{-8, 8}},
            .points = 64
        },
        .initial_conditions = .{
            .mass = 1,
            .momentum = &.{0, 0},
            .position = &.{1, 1},
            .state = 0,
            .gamma = 2
        },
        .potential = .{
            .harmonic = HarmonicPotential(f64){
                .k = &.{1, 1}
            }
        },
        .iterations = 1000,
        .time_step = 0.1
    };

    const output = try run(f64, options, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(@abs(output.kinetic_energy - 1.05452660197613) < TEST_TOLERANCE);
    try std.testing.expect(@abs(output.potential_energy - 1.19533673987723) < TEST_TOLERANCE);
}

test "Exact Nonadiabatic Dynamics on Tully's First Potential" {
    const options = Options(f64){
        .grid = .{
            .limits = &.{&.{-24, 32}},
            .points = 2048
        },
        .initial_conditions = .{
            .mass = 2000,
            .momentum = &.{15},
            .position = &.{-10},
            .state = 1,
            .gamma = 2
        },
        .potential = .{
            .tully_1 = TullyPotential1(f64){}
        },
        .iterations = 350,
        .time_step = 10
    };

    const output = try run(f64, options, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(@abs(output.kinetic_energy - 0.06471011654226) < TEST_TOLERANCE);
    try std.testing.expect(@abs(output.population.at(options.iterations, 0) - 0.58949426578088) < TEST_TOLERANCE);
    try std.testing.expect(@abs(output.population.at(options.iterations, 1) - 0.41050573421579) < TEST_TOLERANCE);
    try std.testing.expect(@abs(output.potential_energy - 0.00178988577130) < TEST_TOLERANCE);
}

test "Imaginary Time Propagation on 1D Harmonic Potential" {
    const options = Options(f64){
        .grid = .{
            .limits = &.{&.{-8, 8}},
            .points = 64
        },
        .initial_conditions = .{
            .mass = 1,
            .momentum = &.{0},
            .position = &.{1},
            .state = 0,
            .gamma = 2
        },
        .potential = .{
            .harmonic = HarmonicPotential(f64){}
        },
        .iterations = 1000,
        .time_step = 0.1,
        .imaginary = true
    };

    const output = try run(f64, options, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(@abs(output.kinetic_energy - 0.25031230493126) < TEST_TOLERANCE);
    try std.testing.expect(@abs(output.potential_energy - 0.24968808471946) < TEST_TOLERANCE);
}

test "Imaginary Time Propagation on 2D Harmonic Potential" {
    const options = Options(f64){
        .grid = .{
            .limits = &.{&.{-8, 8}, &.{-8, 8}},
            .points = 64
        },
        .initial_conditions = .{
            .mass = 1,
            .momentum = &.{0, 0},
            .position = &.{1, 0},
            .state = 0,
            .gamma = 2
        },
        .potential = .{
            .harmonic = HarmonicPotential(f64){
                .k = &.{1, 1}
            }
        },
        .iterations = 1000,
        .time_step = 0.1,
        .imaginary = true
    };

    const output = try run(f64, options, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(@abs(output.kinetic_energy - 0.50062460986252) < TEST_TOLERANCE);
    try std.testing.expect(@abs(output.potential_energy - 0.49937616943892) < TEST_TOLERANCE);
}
