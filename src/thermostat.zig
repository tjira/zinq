//! Implements thermostat algorithms for canonical ensemble molecular dynamics.

const std = @import("std");

const Matrix = @import("tensor.zig").Matrix;

const AU2K = @import("constant.zig").AU2K;

/// Tagged union for thermostat configuration options.
pub const Options = union(enum) {
    berendsen: BerendsenOptions,
    langevin: LangevinOptions,
};

/// Configuration options for the Berendsen velocity rescaling thermostat.
pub const BerendsenOptions = struct {
    tau: f64 = 100,
    temperature: f64,
};

/// Configuration options for the Langevin thermostat.
pub const LangevinOptions = struct {
    gamma: f64 = 1,
    seed: u32 = 1,
    temperature: f64,
};

/// Implements Berendsen velocity rescaling thermostat for canonical temperature coupling.
pub fn Berendsen(comptime T: type) type {
    return struct {
        dt_over_tau: T,
        target_temp: T,

        /// Initializes Berendsen thermostat parameters from configuration options.
        pub fn init(opt: BerendsenOptions, dt: T) @This() {
            return .{ .dt_over_tau = dt / opt.tau, .target_temp = opt.temperature / AU2K };
        }

        /// Applies velocity rescaling to trajectory momenta to couple with thermal bath.
        pub fn apply(self: *@This(), p: *Matrix(T), m: []const T) void {
            const ndim = p.ncol();

            for (0..p.nrow()) |i| {
                var ekin: T = 0;

                for (0..ndim) |j| {
                    ekin += (p.at(i, j) * p.at(i, j)) / (2 * m[j]);
                }

                const t_inst = 2 * ekin / @as(T, @floatFromInt(ndim));

                if (t_inst > 1e-12) {
                    const factor = 1 + self.dt_over_tau * (self.target_temp / t_inst - 1);

                    for (0..ndim) |j| {
                        p.ptr(i, j).* *= @sqrt(@max(0, factor));
                    }
                }
            }
        }
    };
}

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
        berendsen: Berendsen(T),
        langevin: Langevin(T),

        /// Initializes the configured thermostat method.
        pub fn init(opt: Options, dt: T) @This() {
            return switch (opt) {
                .berendsen => |bopt| .{ .berendsen = Berendsen(T).init(bopt, dt) },
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
