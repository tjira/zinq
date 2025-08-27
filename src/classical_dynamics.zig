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
            trajectory: u32 = 1, iteration: u32 = 1
        };

        electronic_potential: ElectronicPotential(T),
        initial_conditions: InitialConditions,

        trajectories: u32,
        iterations: u32,
        time_step: T,

        log_intervals: LogIntervals = .{},

        finite_differences_step: T = 1e-8,
        seed: u32 = 1
    };
}

/// Run classical dynamics simulation.
pub fn run(comptime T: type, options: Options(T), allocator: std.mem.Allocator) !void {
    const ndim = options.electronic_potential.ndim();
    const nstate = options.electronic_potential.nstate();

    var system = try ClassicalParticle(T).initZero(ndim, options.initial_conditions.mass, allocator); defer system.deinit();

    try system.setPositionRandn(options.initial_conditions.position_mean, options.initial_conditions.position_std, options.seed + 0);
    try system.setMomentumRandn(options.initial_conditions.momentum_mean, options.initial_conditions.momentum_std, options.seed + 1);

    var U = try RealMatrix(T).init(nstate, nstate, allocator); defer U.deinit();

    try options.electronic_potential.eval(&U, system, 0);

    try printRealMatrix(T, U);
}
