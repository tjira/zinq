//! Code to propagate classical particles.

const std = @import("std");

const classical_particle = @import("classical_particle.zig");
const device_write = @import("device_write.zig");
const electronic_potential = @import("electronic_potential.zig");
const linear_algebra = @import("linear_algebra.zig");
const math_functions = @import("math_functions.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const ClassicalParticle = classical_particle.ClassicalParticle;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const eigensystemSymmetric = linear_algebra.eigensystemSymmetric;
const print = device_write.print;
const printRealMatrix = device_write.printRealMatrix;
const printRealVector = device_write.printRealVector;
const sgn = math_functions.sgn;

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

        seed: u32 = 0,
    };
}

/// Structure that hold info about the current iteration.
pub fn IterationInfo(comptime T: type) type {
    return struct {
        iteration: usize,
        kinetic_energy: T,
        potential_energy: T,
        state: u32,
        system: ClassicalParticle(T),
        time: T,
        trajectory: usize,
    };
}

/// Run classical dynamics simulation.
pub fn run(comptime T: type, options: Options(T), allocator: std.mem.Allocator) !void {
    const ndim = options.potential.ndim();

    var rng = std.Random.DefaultPrng.init(options.seed); var random = rng.random();

    try print("\nRUNNING CLASSICAL DYNAMICS WITH {d} TRAJECTORIES OF {d} ITERATIONS EACH\n\n", .{options.trajectories, options.iterations});

    try printIterationHeader(ndim);

    for (0..options.trajectories) |i| {

        var system = try ClassicalParticle(T).initZero(ndim, options.initial_conditions.mass, allocator); defer system.deinit();

        try system.setPositionRandn(options.initial_conditions.position_mean, options.initial_conditions.position_std, &random);
        try system.setMomentumRandn(options.initial_conditions.momentum_mean, options.initial_conditions.momentum_std, &random);

        try runTrajectory(T, options, &system, i, allocator);
    }
}

/// Run a single trajectory.
pub fn runTrajectory(comptime T: type, options: Options(T), system: *ClassicalParticle(T), index: usize, allocator: std.mem.Allocator) !void {
    const nstate = options.potential.nstate();

    const current_state = options.initial_conditions.state;

    var diabatic_potential = try RealMatrix(T).init(nstate, nstate, allocator); defer diabatic_potential.deinit();
    var adiabatic_potential = try RealMatrix(T).init(nstate, nstate, allocator); defer adiabatic_potential.deinit();
    var adiabatic_eigenvectors = try RealMatrix(T).init(nstate, nstate, allocator); defer adiabatic_eigenvectors.deinit();

    var potential_eigensystem = .{
        .diabatic_potential = diabatic_potential,
        .adiabatic_potential = adiabatic_potential,
        .adiabatic_eigenvectors = adiabatic_eigenvectors
    };

    for (0..options.iterations + 1) |i| {

        const time = @as(T, @floatFromInt(i)) * options.time_step;

        if (i > 0) {
            try system.propagateVelocityVerlet(options.potential, &adiabatic_potential, time, current_state, options.time_step);
        }

        try options.potential.evaluateEigensystem(&potential_eigensystem, system.position, time);

        const kinetic_energy = system.kineticEnergy();
        const potential_energy = adiabatic_potential.at(current_state, current_state);

        if (index % options.log_intervals.trajectory != 0 or (i > 0 and i % options.log_intervals.iteration != 0)) continue;

        const iteration_info = IterationInfo(T){
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
pub fn printIterationInfo(comptime T: type, info: IterationInfo(T)) !void {
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
