const std = @import("std");

const Allocator = std.mem.Allocator;

// zig fmt: off
const Ensemble       = @import("classical_dynamics.zig").      Ensemble;
const GradientBuffer = @import("classical_dynamics.zig").GradientBuffer;
const Matrix         = @import("tensor.zig"            ).        Matrix;
const Potential      = @import("potential.zig"         ).     Potential;
const Vector         = @import("tensor.zig"            ).        Vector;
// zig fmt: on

// OPTIONS =============================================================================================================

pub const Options = union(enum) {
    // zig fmt: off
    landau_zener: LandauZenerOptions,
    // zig fmt: on
};

pub const LandauZenerOptions = struct {
    seed: u32 = 1,
};

// GENERIC SURFACE HOPPING =============================================================================================

pub fn SurfaceHopping(comptime T: type) type {
    return struct {
        // zig fmt: off
        rng: std.Random.DefaultPrng, probs: Matrix(T), method: Method, adia: bool,
        // zig fmt: on

        pub const Method = union(enum) {
            landau_zener: LandauZener(T),
        };

        pub fn init(options: Options, nstate: usize, ntraj: usize, adia: bool, gpa: Allocator) !@This() {
            const seed = switch (options) {
                .landau_zener => |field| field.seed,
            };

            var split_mix = std.Random.SplitMix64.init(seed);

            const method: Method = switch (options) {
                .landau_zener => .{ .landau_zener = try LandauZener(T).init(nstate, ntraj, adia, gpa) },
            };

            const probs = try Matrix(T).init(ntraj, nstate, gpa);

            return .{ .rng = std.Random.DefaultPrng.init(split_mix.next()), .probs = probs, .method = method, .adia = adia };
        }

        pub fn hop(self: *@This(), ensemble: *Ensemble(T), V: Matrix(T), W: Matrix(T), _: Matrix(T), dt: T) !void {
            self.update(V, W);
            self.probs.zero();

            switch (self.method) {
                inline else => |*field| try field.calcProbs(&self.probs, ensemble, dt),
            }

            for (0..ensemble.s.length()) |i| {
                // zig fmt: off
                var sum: T = 0; var accum: T = 0;
                // zig fmt: on

                for (0..self.probs.ncol()) |j| {
                    sum += self.probs.at(i, j);
                }

                if (sum > 1) for (0..self.probs.ncol()) |j| {
                    self.probs.ptr(i, j).* /= sum;
                };

                const rv = self.rng.random().float(T);

                // zig fmt: off
                const c = ensemble.s.at(i); var ns = c;
                // zig fmt: on

                for (0..self.probs.ncol()) |j| {
                    accum += self.probs.at(i, j);

                    if (rv < accum) {
                        // zig fmt: off
                        ns = @intCast(j); break;
                        // zig fmt: on
                    }
                }

                if (self.adia and ns != c) {
                    // zig fmt: off
                    const E_new = W.at(i, ns);
                    const E_old = W.at(i, c );
                    // zig fmt: on

                    if (rescaleMomentumIsotropic(T, ensemble, i, E_new - E_old)) {
                        ensemble.s.ptr(i).* = ns;
                    }
                }

                if (!self.adia and ns != c) {
                    ensemble.s.ptr(i).* = ns;
                }
            }
        }

        pub fn update(self: *@This(), V: Matrix(T), W: Matrix(T)) void {
            switch (self.method) {
                .landau_zener => |*field| field.update(V, W),
            }
        }
    };
}

// SPECIFIC METHODS ====================================================================================================

pub fn LandauZener(comptime T: type) type {
    return struct {
        // zig fmt: off
        history: [3]Matrix(T), adia: bool,
        // zig fmt: on

        pub fn init(nstate: usize, ntraj: usize, adia: bool, gpa: Allocator) !@This() {
            var history: [3]Matrix(T) = undefined;

            const cols = if (adia) nstate else nstate * nstate;

            history[0] = try Matrix(T).init(ntraj, cols, gpa);
            history[1] = try Matrix(T).init(ntraj, cols, gpa);
            history[2] = try Matrix(T).init(ntraj, cols, gpa);

            history[0].fill(std.math.nan(T));
            history[1].fill(std.math.nan(T));
            history[2].fill(std.math.nan(T));

            return .{ .history = history, .adia = adia };
        }

        pub fn calcProbs(self: *@This(), probs: *Matrix(T), ensemble: *Ensemble(T), dt: T) !void {
            if (self.adia) try self.calcProbsAdia(probs, ensemble, dt) else try self.calcProbsDia(probs, ensemble, dt);
        }

        pub fn calcProbsAdia(self: *@This(), probs: *Matrix(T), ensemble: *Ensemble(T), dt: T) !void {
            if (std.math.isNan(self.history[0].at(0, 0))) return;

            for (0..self.history[0].nrow()) |i| {
                const c = ensemble.s.at(i);

                for (0..self.history[0].ncol()) |j| {
                    if (c == j) continue;

                    const z0 = @abs(self.history[0].at(i, j) - self.history[0].at(i, c));
                    const z1 = @abs(self.history[1].at(i, j) - self.history[1].at(i, c));
                    const z2 = @abs(self.history[2].at(i, j) - self.history[2].at(i, c));

                    const dz1 = (z1 - z0) / dt;
                    const dz2 = (z2 - z1) / dt;

                    if (dz1 < 0 and dz2 > 0) {
                        const d2z = (dz2 - dz1) / dt;

                        if (d2z > 1e-14) {
                            // zig fmt: off
                            const b = (z2 - z0) / (2 * dt); const a = d2z / 2; var z_min = z1 - (b * b) / (4 * a);
                            // zig fmt: on

                            if (z_min < 0) z_min = 0;

                            const p = std.math.exp(-(std.math.pi / 2.0) * std.math.sqrt((z_min * z_min * z_min) / d2z));

                            probs.ptr(i, j).* = p;
                        }
                    }
                }
            }
        }

        pub fn calcProbsDia(self: *@This(), probs: *Matrix(T), ensemble: *Ensemble(T), dt: T) !void {
            if (std.math.isNan(self.history[1].at(0, 0))) return;

            for (0..self.history[0].nrow()) |i| {
                const c = ensemble.s.at(i);

                for (0..ensemble.nstate) |j| {
                    if (c == j) continue;

                    const v_cc_old = self.history[1].at(i, c * ensemble.nstate + c);
                    const v_cc_new = self.history[2].at(i, c * ensemble.nstate + c);
                    const v_jj_old = self.history[1].at(i, j * ensemble.nstate + j);
                    const v_jj_new = self.history[2].at(i, j * ensemble.nstate + j);

                    const gap_old  = v_cc_old - v_jj_old;
                    const gap_new  = v_cc_new - v_jj_new;

                    if ((gap_old > 0 and gap_new <= 0) or (gap_old < 0 and gap_new >= 0)) {
                        const d_gap = @abs(gap_new - gap_old) / dt;

                        if (d_gap > 1e-14) {
                            const v_cj = @abs(self.history[2].at(i, c * ensemble.nstate + j));

                            probs.ptr(i, j).* = 1 - std.math.exp(-2 * std.math.pi * v_cj * v_cj / d_gap);
                        }
                    }
                }
            }
        }

        pub fn update(self: *@This(), V: Matrix(T), W: Matrix(T)) void {
            const data = if (self.adia) W else V;

            for (0..data.nrow()) |i| for (0..data.ncol()) |j| {
                self.history[0].ptr(i, j).* = self.history[1].at(i, j);
                self.history[1].ptr(i, j).* = self.history[2].at(i, j);

                self.history[2].ptr(i, j).* = data.at(i, j);
            };
        }
    };
}

// HELPER FUNCTIONS ====================================================================================================

fn rescaleMomentumIsotropic(comptime T: type, ensemble: *Ensemble(T), i: usize, dE: T) bool {
    var ekin_old: T = 0;

    for (0..ensemble.p.ncol()) |j| {
        ekin_old += ensemble.p.at(i, j) * ensemble.p.at(i, j) / (2 * ensemble.mass);
    }

    const ekin_new = ekin_old - dE;

    if (ekin_new >= 0) {
        for (0..ensemble.p.ncol()) |j| {
            ensemble.p.ptr(i, j).* *= std.math.sqrt(ekin_new / ekin_old);
        }

        return true;
    }

    for (0..ensemble.p.ncol()) |j| {
        ensemble.p.ptr(i, j).* *= -1;
    }

    return false;
}
