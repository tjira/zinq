//! Code to propagate classical particles.

const std = @import("std");

const classical_particle = @import("classical_particle.zig");
const complex_runge_kutta = @import("complex_runge_kutta.zig");
const complex_vector = @import("complex_vector.zig");
const derivative_coupling = @import("derivative_coupling.zig");
const device_write = @import("device_write.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const electronic_potential = @import("electronic_potential.zig");
const error_handling = @import("error_handling.zig");
const fewest_switches = @import("fewest_switches.zig");
const harmonic_potential = @import("harmonic_potential.zig");
const landau_zener = @import("landau_zener.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const norm_preserving_interpolation = @import("norm_preserving_interpolation.zig");
const object_array = @import("object_array.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");
const ring_buffer = @import("ring_buffer.zig");
const surface_hopping_algorithm = @import("surface_hopping_algorithm.zig");
const tully_potential = @import("tully_potential.zig");

const ClassicalParticle = classical_particle.ClassicalParticle;
const Complex = std.math.complex.Complex;
const ComplexRungeKutta = complex_runge_kutta.ComplexRungeKutta;
const ComplexVector = complex_vector.ComplexVector;
const DerivativeCoupling = derivative_coupling.DerivativeCoupling;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const FewestSwitches = fewest_switches.FewestSwitches;
const HarmonicPotential = harmonic_potential.HarmonicPotential;
const LandauZener = landau_zener.LandauZener;
const NormPreservingInterpolation = norm_preserving_interpolation.NormPreservingInterpolation;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;
const RingBufferArray = object_array.RingBufferArray;
const SurfaceHoppingAlgorithm = surface_hopping_algorithm.SurfaceHoppingAlgorithm;
const TullyPotential1 = tully_potential.TullyPotential1;

const exportRealMatrix = device_write.exportRealMatrix;
const exportRealMatrixWithLinspacedLeftColumn = device_write.exportRealMatrixWithLinspacedLeftColumn;
const fixGauge = eigenproblem_solver.fixGauge;
const mm = matrix_multiplication.mm;
const print = device_write.print;
const printJson = device_write.printJson;
const throw = error_handling.throw;

/// Classical dynamics option struct.
pub fn Options(comptime T: type) type {
    return struct {
        pub const InitialConditions = struct {
            mass: []const T,
            momentum_mean: []const T,
            momentum_std: []const T,
            position_mean: []const T,
            position_std: []const T,
            state: u32
        };
        pub const LogIntervals = struct {
            trajectory: u32 = 1,
            iteration: u32 = 1
        };
        pub const Write = struct {
            population_mean: ?[]const u8 = null
        };

        potential: ElectronicPotential(T),
        initial_conditions: InitialConditions,

        trajectories: u32,
        iterations: u32,
        time_step: T,

        derivative_coupling: ?DerivativeCoupling(T) = null,
        log_intervals: LogIntervals = .{},
        surface_hopping: ?SurfaceHoppingAlgorithm(T) = null,
        write: Write = .{},

        seed: u32 = 0,
    };
}

/// Structure that hold the output of the simulation.
pub fn Output(comptime T: type) type {
    return struct {
        population_mean: RealMatrix(T),

        /// Allocate the output structure.
        pub fn init(nstate: usize, iterations: usize, allocator: std.mem.Allocator) !@This() {
            return @This(){
                .population_mean = try RealMatrix(T).initZero(iterations + 1, nstate, allocator)
            };
        }

        /// Free the memory allocated for the output structure.
        pub fn deinit(self: @This()) void {
            self.population_mean.deinit();
        }
    };
}

/// Container for custom structs related to the classical dynamics.
pub fn Custom(comptime T: type) type {
    return struct {

        /// Structure to hold information about each iteration used for logging.
        pub const IterationInfo = struct {
            iteration: usize,
            kinetic_energy: T,
            potential_energy: T,
            state: usize,
            system: ClassicalParticle(T),
            time: T,
            trajectory: usize,
            amplitudes: ComplexVector(T)
        };

        /// Structure to hold a single trajectory output.
        pub const TrajectoryOutput = struct {
            population: RealMatrix(T),

            /// Allocate the trajectory output structure.
            pub fn init(nstate: usize, iterations: usize, allocator: std.mem.Allocator) !@This() {
                return @This(){
                    .population = try RealMatrix(T).initZero(iterations + 1, nstate, allocator)
                };
            }

            /// Free the memory allocated for the trajectory output structure.
            pub fn deinit(self: @This()) void {
                self.population.deinit();
            }
        };
    };
}

/// Run classical dynamics simulation.
pub fn run(comptime T: type, options: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(options);

    var potential = options.potential;
    const ndim = potential.ndim();
    const nstate = potential.nstate();

    const file_potential = if (potential == .file) try potential.file.init(allocator) else null; defer if (file_potential) |U| U.deinit();

    var output = try Output(T).init(nstate, options.iterations, allocator);

    var split_mix = std.Random.SplitMix64.init(options.seed); var rng = std.Random.DefaultPrng.init(split_mix.next()); var random = rng.random();

    if (enable_printing) try printIterationHeader(T, ndim, nstate, options.surface_hopping);

    for (0..options.trajectories) |i| {

        var system = try ClassicalParticle(T).initZero(ndim, options.initial_conditions.mass, allocator); defer system.deinit();

        try system.setPositionRandn(options.initial_conditions.position_mean, options.initial_conditions.position_std, &random);
        try system.setMomentumRandn(options.initial_conditions.momentum_mean, options.initial_conditions.momentum_std, &random);

        var options_copy = options; options_copy.potential = potential;

        const trajectory_output = try runTrajectory(T, options_copy, &system, i, enable_printing, allocator); defer trajectory_output.deinit();

        try output.population_mean.add(trajectory_output.population);
    }

    output.population_mean.divs(@as(T, @floatFromInt(options.trajectories)));

    const population_error = 1.96 * std.math.sqrt(output.population_mean.at(options.iterations, 0) * (1 - output.population_mean.at(options.iterations, 0)) / @as(T, @floatFromInt(options.trajectories)));

    if (enable_printing) for (0..nstate) |i| {
        try print("{s}FINAL POPULATION OF STATE {d:2}: {d:.6} ± {:.6}\n", .{if (i == 0) "\n" else "", i, output.population_mean.at(options.iterations, i), population_error});
    };

    const end_time = @as(T, @floatFromInt(options.iterations)) * options.time_step;

    if (options.write.population_mean) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.population_mean, 0, end_time, options.iterations + 1);

    return output;
}

/// Run a single trajectory.
pub fn runTrajectory(comptime T: type, options: Options(T), system: *ClassicalParticle(T), index: usize, enable_printing: bool, allocator: std.mem.Allocator) !Custom(T).TrajectoryOutput {
    const nstate = options.potential.nstate();

    var output = try Custom(T).TrajectoryOutput.init(nstate, options.iterations, allocator);

    var current_state: usize = @intCast(options.initial_conditions.state);

    if (current_state >= nstate) return throw(Custom(T).TrajectoryOutput, "ACTIVE STATE MUST NOT BE HIGHER THAN THE TOTAL NUMBER OF STATES", .{});

    var diabatic_potential = try RealMatrix(T).init(nstate, nstate, allocator); defer diabatic_potential.deinit();
    var adiabatic_potential = try RealMatrix(T).init(nstate, nstate, allocator); defer adiabatic_potential.deinit();
    var adiabatic_eigenvectors = try RealMatrix(T).init(nstate, nstate, allocator); defer adiabatic_eigenvectors.deinit();
    var previous_eigenvectors = try RealMatrix(T).init(nstate, nstate, allocator); defer previous_eigenvectors.deinit();
    var eigenvector_overlap = try RealMatrix(T).init(nstate, nstate, allocator); defer eigenvector_overlap.deinit();
    var time_derivative_coupling = try RealMatrix(T).initZero(nstate, nstate, allocator); defer time_derivative_coupling.deinit();

    var energy_gaps = try RingBufferArray(T).init((2 * nstate - 2) / 2, .{.max_len = 5}, allocator); defer energy_gaps.deinit();
    var jump_probabilities = try RealVector(T).init(nstate, allocator); defer jump_probabilities.deinit();
    var runge_kutta_solver = try ComplexRungeKutta(T).init(nstate, allocator); defer runge_kutta_solver.deinit();
    var amplitudes = try ComplexVector(T).initZero(nstate, allocator); defer amplitudes.deinit();

    amplitudes.ptr(current_state).* = Complex(T).init(1, 0);

    const fs_parameters: fewest_switches.Parameters(T) = .{
        .adiabatic_potential = adiabatic_potential,
        .amplitudes = &amplitudes,
        .derivative_coupling = time_derivative_coupling,
        .runge_kutta = runge_kutta_solver,
        .time_step = options.time_step
    };

    const lz_parameters: landau_zener.Parameters(T) = .{
        .energy_gaps = energy_gaps,
        .time_step = options.time_step
    };

    const surface_hopping_parameters: surface_hopping_algorithm.Parameters(T) = .{
        .fs_parameters = fs_parameters,
        .lz_parameters = lz_parameters
    };

    var split_mix = std.Random.SplitMix64.init(options.seed + index); var rng = std.Random.DefaultPrng.init(split_mix.next()); var random = rng.random();

    var timer = try std.time.Timer.start();

    for (0..options.iterations + 1) |i| {

        const time = @as(T, @floatFromInt(i)) * options.time_step;

        try adiabatic_eigenvectors.copyTo(&previous_eigenvectors);

        if (i > 0) {
            try system.propagateVelocityVerlet(options.potential, &adiabatic_potential, time, current_state, options.time_step);
        }

        try options.potential.evaluateEigensystem(&diabatic_potential, &adiabatic_potential, &adiabatic_eigenvectors, system.position, time);

        if (i > 0) {

            try fixGauge(T, &adiabatic_eigenvectors, previous_eigenvectors);
            try mm(T, &eigenvector_overlap, previous_eigenvectors, true, adiabatic_eigenvectors, false);

            if (options.derivative_coupling) |tdc_algorithm| {
                try tdc_algorithm.evaluate(&time_derivative_coupling, eigenvector_overlap, options.time_step);
            }
        }

        for (0..nstate) |j| for (j + 1..nstate) |k| {
            energy_gaps.ptr(j + k - 1).append(adiabatic_potential.at(k, k) - adiabatic_potential.at(j, j));
        };

        if (options.surface_hopping) |algorithm| if (i > 1) {
            current_state = algorithm.jump(system, &jump_probabilities, surface_hopping_parameters, adiabatic_potential, current_state, &random);
        };

        const kinetic_energy = system.kineticEnergy();
        const potential_energy = adiabatic_potential.at(current_state, current_state);

        output.population.ptr(i, current_state).* += 1;

        if (!enable_printing or (index > 0 and (index + 1) % options.log_intervals.trajectory != 0) or (i > 0 and i % options.log_intervals.iteration != 0)) continue;

        const iteration_info = Custom(T).IterationInfo{
            .iteration = i,
            .kinetic_energy = kinetic_energy,
            .potential_energy = potential_energy,
            .state = current_state,
            .system = system.*,
            .time = time,
            .trajectory = index,
            .amplitudes = amplitudes
        };

        try printIterationInfo(T, iteration_info, options.surface_hopping, &timer);
    }

    return output;
}

/// Print header for iteration info.
pub fn printIterationHeader(comptime T: type, ndim: usize, nstate: usize, surface_hopping: ?SurfaceHoppingAlgorithm(T)) !void {
    var buffer: [128]u8 = undefined;

    var writer = std.io.Writer.fixed(&buffer);

    try writer.print("\n{s:8} {s:8} ", .{"TRAJ", "ITER"});
    try writer.print("{s:12} {s:12} {s:12} ", .{"KINETIC", "POTENTIAL", "TOTAL"});
    try writer.print("{s:5} ", .{"STATE"});
    try writer.print("{[value]s:[width]} ", .{.value = "POSITION", .width = 9 * ndim + 2 * (ndim - 1) + 2});
    try writer.print("{[value]s:[width]} ", .{.value = "MOMENTUM", .width = 9 * ndim + 2 * (ndim - 1) + 2});

    if (surface_hopping) |algorithm| switch (algorithm) {
        .fewest_switches => try writer.print("{[value]s:[width]} ", .{.value = "|COEFS|^2", .width = 7 * nstate + 2 * (nstate - 1) + 2}),
        else => {}
    };

    try writer.print("{s:4}", .{"TIME"});

    try print("{s}\n", .{writer.buffered()});
}

/// Prints the iteration info to standard output.
pub fn printIterationInfo(comptime T: type, info: Custom(T).IterationInfo, surface_hopping: ?SurfaceHoppingAlgorithm(T), timer: *std.time.Timer) !void {
    var buffer: [128]u8 = undefined;

    var writer = std.io.Writer.fixed(&buffer);

    try writer.print("{d:8} {d:8} ", .{info.trajectory + 1, info.iteration});
    try writer.print("{d:12.6} {d:12.6} {d:12.6} ", .{info.kinetic_energy, info.potential_energy, info.kinetic_energy + info.potential_energy});
    try writer.print("{d:5} ", .{info.state});

    try writer.print("[", .{});

    for (0..info.system.ndim) |i| {
        try writer.print("{d:9.4}{s}", .{info.system.position.at(i), if (i == info.system.ndim - 1) "" else ", "});
    }
    
    try writer.print("] [", .{});

    for (0..info.system.ndim) |i| {
        try writer.print("{d:9.4}{s}", .{info.system.velocity.at(i) * info.system.masses.at(i), if (i == info.system.ndim - 1) "" else ", "});
    }

    if (surface_hopping) |algorithm| {

        if (algorithm == .fewest_switches) try writer.print("] [", .{});

        if (algorithm == .fewest_switches) for (0..info.amplitudes.len) |i| {
            try writer.print("{d:7.4}{s}", .{std.math.pow(T, info.amplitudes.at(i).magnitude(), 2), if (i == info.amplitudes.len - 1) "" else ", "});
        };
    }

    try writer.print("] {D}", .{timer.read()}); timer.reset();

    try print("{s}\n", .{writer.buffered()});
}

test "Fewest Switches Surface Hopping on Tully's First Potential" {
    const options = Options(f64){
        .derivative_coupling = .{
            .npi = NormPreservingInterpolation(f64){}
        },
        .initial_conditions = .{
            .mass = &.{2000},
            .momentum_mean = &.{15},
            .momentum_std = &.{1},
            .position_mean = &.{-10},
            .position_std = &.{0.5},
            .state = 1
        },
        .potential = .{
            .tully_1 = TullyPotential1(f64){}
        },
        .surface_hopping = .{
            .fewest_switches = FewestSwitches(f64){}
        },
        .iterations = 3500,
        .time_step = 1,
        .trajectories = 1000
    };

    const output = try run(f64, options, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(output.population_mean.at(options.iterations, 0) == 0.379);
    try std.testing.expect(output.population_mean.at(options.iterations, 1) == 0.621);
}

test "Landau-Lener Surface Hopping on Tully's First Potential" {
    const options = Options(f64){
        .initial_conditions = .{
            .mass = &.{2000},
            .momentum_mean = &.{15},
            .momentum_std = &.{1},
            .position_mean = &.{-10},
            .position_std = &.{0.5},
            .state = 1
        },
        .potential = .{
            .tully_1 = TullyPotential1(f64){}
        },
        .surface_hopping = .{
            .landau_zener = LandauZener(f64){}
        },
        .iterations = 3500,
        .time_step = 1,
        .trajectories = 1000
    };

    const output = try run(f64, options, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(output.population_mean.at(options.iterations, 0) == 0.498);
    try std.testing.expect(output.population_mean.at(options.iterations, 1) == 0.502);
}
