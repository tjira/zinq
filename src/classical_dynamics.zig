//! Code to propagate classical particles.

const std = @import("std");

const classical_particle = @import("classical_particle.zig");
const device_write = @import("device_write.zig");
const electronic_potential = @import("electronic_potential.zig");
const math_functions = @import("math_functions.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const ClassicalParticle = classical_particle.ClassicalParticle;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

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

        electronic_potential: ElectronicPotential(T),
        initial_conditions: InitialConditions,

        trajectories: u32,
        iterations: u32,
        time_step: T,

        log_intervals: LogIntervals = .{},

        seed: u32 = 0,
    };
}

/// Run classical dynamics simulation.
pub fn run(comptime T: type, options: Options(T), allocator: std.mem.Allocator) !void {
    const ndim = options.electronic_potential.ndim();

    var rng = std.Random.DefaultPrng.init(options.seed); var random = rng.random();

    for (0..options.trajectories) |i| {

        var system = try ClassicalParticle(T).initZero(ndim, options.initial_conditions.mass, allocator); defer system.deinit();

        try system.setPositionRandn(options.initial_conditions.position_mean, options.initial_conditions.position_std, &random);
        try system.setMomentumRandn(options.initial_conditions.momentum_mean, options.initial_conditions.momentum_std, &random);

        try runTrajectory(T, options, &system, i, allocator);
    }
}

/// Run a single trajectory.
pub fn runTrajectory(comptime T: type, options: Options(T), system: *ClassicalParticle(T), index: usize, allocator: std.mem.Allocator) !void {
    // const ndim = options.electronic_potential.ndim();
    const nstate = options.electronic_potential.nstate();

    const current_state = options.initial_conditions.state;

    // var diabatic_potential = try RealMatrix(T).init(nstate, nstate, allocator); defer diabatic_potential.deinit();
    var adiabatic_potential = try RealMatrix(T).init(nstate, nstate, allocator); defer adiabatic_potential.deinit();
    // var potential_eigenvectors = try RealMatrix(T).init(nstate, nstate, allocator); defer potential_eigenvectors.deinit();

    for (0..options.iterations) |i| {

        const time = @as(T, @floatFromInt(i)) * options.time_step;

        if (i > 0) {
            try system.propagateVelocityVerlet(options.electronic_potential, current_state, options.time_step);
        }

        options.electronic_potential.evaluateAdiabatic(&adiabatic_potential, system.*, time);

        std.debug.print("{d:8} {d:8} {d} {d} {d}\n", .{index + 1, i,adiabatic_potential.at(0, 0), system.position.at(0), system.velocity.at(0)});
    }
}
