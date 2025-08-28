//! Code to propagate classical particles.

const std = @import("std");

const classical_particle = @import("classical_particle.zig");
const device_write = @import("device_write.zig");
const electronic_potential = @import("electronic_potential.zig");
const landau_zener = @import("landau_zener.zig");
const linear_algebra = @import("linear_algebra.zig");
const object_array = @import("object_array.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");
const ring_buffer = @import("ring_buffer.zig");
const surface_hopping_algorithm = @import("surface_hopping_algorithm.zig");

const ClassicalParticle = classical_particle.ClassicalParticle;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;
const RingBufferArray = object_array.RingBufferArray;
const SurfaceHoppingAlgorithm = surface_hopping_algorithm.SurfaceHoppingAlgorithm;

const print = device_write.print;

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

        potential: ElectronicPotential(T),
        initial_conditions: InitialConditions,

        trajectories: u32,
        iterations: u32,
        time_step: T,

        log_intervals: LogIntervals = .{},
        surface_hopping: ?SurfaceHoppingAlgorithm(T) = null,

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
                .population_mean = try RealMatrix(T).init(iterations + 1, nstate, allocator)
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
            state: u32,
            system: ClassicalParticle(T),
            time: T,
            trajectory: usize,
        };

        /// Structure to hold a single trajectory output.
        pub const TrajectoryOutput = struct {
            population: RealMatrix(T),

            /// Allocate the trajectory output structure.
            pub fn init(nstate: usize, iterations: usize, allocator: std.mem.Allocator) !@This() {
                return @This(){
                    .population = try RealMatrix(T).init(iterations + 1, nstate, allocator)
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
pub fn run(comptime T: type, options: Options(T), allocator: std.mem.Allocator) !*anyopaque {
    try print("\nRUNNING CLASSICAL DYNAMICS WITH {d} TRAJECTORIES OF {d} ITERATIONS EACH\n\n", .{options.trajectories, options.iterations});

    const ndim = options.potential.ndim();
    const nstate = options.potential.nstate();

    var output = try Output(T).init(nstate, options.iterations, allocator);

    var split_mix = std.Random.SplitMix64.init(options.seed); var rng = std.Random.DefaultPrng.init(split_mix.next()); var random = rng.random();

    try printIterationHeader(ndim);

    for (0..options.trajectories) |i| {

        var system = try ClassicalParticle(T).initZero(ndim, options.initial_conditions.mass, allocator); defer system.deinit();

        try system.setPositionRandn(options.initial_conditions.position_mean, options.initial_conditions.position_std, &random);
        try system.setMomentumRandn(options.initial_conditions.momentum_mean, options.initial_conditions.momentum_std, &random);

        const trajectory_output = try runTrajectory(T, options, &system, i, allocator); defer trajectory_output.deinit();

        output.population_mean.add(trajectory_output.population);
    }

    output.population_mean.divs(@as(T, @floatFromInt(options.trajectories)));

    for (0..nstate) |i| {
        try print("{s}FINAL POPULATION OF STATE {d}: {d:.6}\n", .{if (i == 0) "\n" else "", i, output.population_mean.at(options.iterations, i)});
    }

    return @ptrCast(&output);
}

/// Run a single trajectory.
pub fn runTrajectory(comptime T: type, options: Options(T), system: *ClassicalParticle(T), index: usize, allocator: std.mem.Allocator) !Custom(T).TrajectoryOutput {
    const nstate = options.potential.nstate();

    var output = try Custom(T).TrajectoryOutput.init(nstate, options.iterations, allocator);

    var current_state = options.initial_conditions.state;

    var diabatic_potential = try RealMatrix(T).init(nstate, nstate, allocator); defer diabatic_potential.deinit();
    var adiabatic_potential = try RealMatrix(T).init(nstate, nstate, allocator); defer adiabatic_potential.deinit();
    var adiabatic_eigenvectors = try RealMatrix(T).init(nstate, nstate, allocator); defer adiabatic_eigenvectors.deinit();

    var potential_eigensystem = .{
        .diabatic_potential = diabatic_potential,
        .adiabatic_potential = adiabatic_potential,
        .adiabatic_eigenvectors = adiabatic_eigenvectors
    };

    var energy_gaps = try RingBufferArray(T).init((2 * nstate - 2) / 2, .{.max_len = 5}, allocator); defer energy_gaps.deinit();
    var jump_probabilities = try RealVector(T).init(nstate, allocator); defer jump_probabilities.deinit();

    const lz_parameters: landau_zener.Parameters(T) = .{
        .energy_gaps = energy_gaps,
        .time_step = options.time_step,
    };

    const surface_hopping_parameters: surface_hopping_algorithm.Parameters(T) = .{
        .adiabatic_potential = adiabatic_potential,
        .lz_parameters = lz_parameters,
    };

    var split_mix = std.Random.SplitMix64.init(options.seed + index); var rng = std.Random.DefaultPrng.init(split_mix.next()); var random = rng.random();

    for (0..options.iterations + 1) |i| {

        const time = @as(T, @floatFromInt(i)) * options.time_step;

        if (i > 0) {
            try system.propagateVelocityVerlet(options.potential, &adiabatic_potential, time, current_state, options.time_step);
        }

        try options.potential.evaluateEigensystem(&potential_eigensystem, system.position, time);

        for (0..nstate) |j| for (j + 1..nstate) |k| {
            energy_gaps.ptr(j + k - 1).append(adiabatic_potential.at(k, k) - adiabatic_potential.at(j, j));
        };

        if (options.surface_hopping) |algorithm| if (i > 2) {
            current_state = algorithm.jump(system, &jump_probabilities, surface_hopping_parameters, current_state, &random);
        };

        const kinetic_energy = system.kineticEnergy();
        const potential_energy = adiabatic_potential.at(current_state, current_state);

        output.population.ptr(i, current_state).* += 1;

        if ((index > 0 and (index + 1) % options.log_intervals.trajectory != 0) or (i > 0 and i % options.log_intervals.iteration != 0)) continue;

        const iteration_info = Custom(T).IterationInfo{
            .iteration = i,
            .kinetic_energy = kinetic_energy,
            .potential_energy = potential_energy,
            .state = current_state,
            .system = system.*,
            .time = time,
            .trajectory = index,
        };

        try printIterationInfo(T, iteration_info);
    }

    return output;
}

/// Print header for iteration info.
pub fn printIterationHeader(ndim: usize) !void {
    var buffer: [128]u8 = undefined;

    var writer = std.io.Writer.fixed(&buffer);

    try writer.print("{s:8} {s:8} ", .{"TRAJ", "ITER"});
    try writer.print("{s:12} {s:12} {s:12} ", .{"KINETIC", "POTENTIAL", "TOTAL"});
    try writer.print("{s:5} ", .{"STATE"});
    try writer.print("{[value]s:[width]} ", .{.value = "POSITION", .width = 9 * ndim + 2 * (ndim - 1) + 2});
    try writer.print("{[value]s:[width]}", .{.value = "MOMENTUM", .width = 9 * ndim + 2 * (ndim - 1) + 2});

    try print("{s}\n", .{writer.buffered()});
}

/// Prints the iteration info to standard output.
pub fn printIterationInfo(comptime T: type, info: Custom(T).IterationInfo) !void {
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

    try writer.print("]", .{});

    try print("{s}\n", .{writer.buffered()});
}
