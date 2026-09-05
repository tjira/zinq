//! Implements thermostat algorithms for canonical ensemble molecular dynamics.

const std = @import("std");

const AU2K = @import("constant.zig").AU2K;
const Matrix = @import("tensor.zig").Matrix;

/// Tagged union for thermostat configuration options.
pub const Options = union(enum) {
    langevin: LangevinOptions,
};

/// Configuration options for the Langevin thermostat.
pub const LangevinOptions = struct {
    temperature: f64,
    gamma: f64 = 1,
    seed: u32 = 1,
};

/// Implements Langevin dynamics thermalization via Ornstein-Uhlenbeck stochastic integration.
pub fn Langevin(comptime T: type) type {
    return struct {
        c1: T,
        c2: T,

        rng: std.Random.DefaultPrng,

        /// Initializes Langevin thermal parameters and pseudo-random number generator.
        pub fn init(opt: LangevinOptions, dt: T) @This() {
            var split_mix = std.Random.SplitMix64.init(opt.seed);
            const rng = std.Random.DefaultPrng.init(split_mix.next());

            const c1 = std.math.exp(-opt.gamma * dt);

            return .{ .c1 = c1, .c2 = std.math.sqrt(opt.temperature * (1 - c1 * c1) / AU2K), .rng = rng };
        }

        /// Applies the stochastic Langevin friction and random thermal kick to trajectory momenta.
        pub fn apply(self: *@This(), p: *Matrix(T), m: []const T) void {
            const random = self.rng.random();

            for (0..p.nrow()) |i| for (0..p.ncol()) |j| {
                p.ptr(i, j).* = self.c1 * p.at(i, j) + self.c2 * @sqrt(m[j]) * random.floatNorm(T);
            };
        }
    };
}

/// Generic thermostat wrapper providing unified interface for canonical ensemble sampling.
pub fn Thermostat(comptime T: type) type {
    return union(enum) {
        langevin: Langevin(T),

        /// Initializes the configured thermostat method.
        pub fn init(opt: Options, dt: T) @This() {
            return switch (opt) {
                .langevin => |lopt| .{ .langevin = Langevin(T).init(lopt, dt) },
            };
        }

        /// Applies the thermalization update to the trajectory ensemble momenta.
        pub fn apply(self: *@This(), p: *Matrix(T), m: []const T) void {
            switch (self.*) {
                inline else => |*t| t.apply(p, m),
            }
        }
    };
}
