const std = @import("std");

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;

const Ensemble = @import("classical_dynamics.zig").Ensemble;
const GradientBuffer = @import("classical_dynamics.zig").GradientBuffer;
const Integrator = @import("integrator.zig").Integrator;
const Matrix = @import("tensor.zig").Matrix;
const Potential = @import("potential.zig").Potential;
const Vector = @import("tensor.zig").Vector;

// OPTIONS =============================================================================================================

pub const Options = union(enum) {
    fewest_switches: FewestSwitchesOptions,
    landau_zener: LandauZenerOptions,
};

pub const FewestSwitchesOptions = struct {
    integrator: std.meta.Tag(Integrator(f64).Method) = .rk4,

    seed: u32 = 1,
    nstep: u32 = 10,
};

pub const LandauZenerOptions = struct {
    seed: u32 = 1,
};

// GENERIC SURFACE HOPPING =============================================================================================

pub fn SurfaceHopping(comptime T: type) type {
    return struct {
        rng: std.Random.DefaultPrng,
        probs: Matrix(T),
        method: Method,

        nstep: usize,
        targets: []usize,
        adia: bool,

        pub const Method = union(enum) {
            fewest_switches: FewestSwitches(T),
            landau_zener: LandauZener(T),
        };

        pub fn init(options: Options, nstate: usize, ntraj: usize, istate: usize, adia: bool, gpa: Allocator) !@This() {
            const seed = switch (options) {
                inline else => |field| field.seed,
            };

            const targets = try gpa.alloc(usize, ntraj);
            errdefer gpa.free(targets);

            const nstep = switch (options) {
                .fewest_switches => |field| field.nstep,
                inline else => 1,
            };

            var split_mix = std.Random.SplitMix64.init(seed);
            const rng = std.Random.DefaultPrng.init(split_mix.next());

            var method: Method = switch (options) {
                .fewest_switches => |opt| .{
                    .fewest_switches = try FewestSwitches(T).init(opt, nstate, ntraj, istate, adia, gpa),
                },
                .landau_zener => |opt| .{
                    .landau_zener = try LandauZener(T).init(opt, nstate, ntraj, istate, adia, gpa),
                },
            };

            errdefer switch (method) {
                inline else => |*field| field.deinit(gpa),
            };

            const probs = try Matrix(T).init(ntraj, nstate, gpa);

            return .{ .rng = rng, .probs = probs, .method = method, .nstep = nstep, .targets = targets, .adia = adia };
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.probs.deinit(gpa);
            gpa.free(self.targets);

            switch (self.method) {
                inline else => |*field| field.deinit(gpa),
            }
        }

        pub fn hop(self: *@This(), ensemble: *Ensemble(T), V: Matrix(T), W: Matrix(T), U: Matrix(T), dt: T) !void {
            self.update(if (self.adia) W else V, U);

            const subdt = dt / @as(T, @floatFromInt(self.nstep));

            for (0..self.targets.len) |i| {
                self.targets[i] = ensemble.s.at(i);
            }

            for (0..self.nstep) |_| {
                try self.calcProbs(&self.probs, ensemble, self.nstep, subdt);

                self.calcTargetStates(ensemble);
            }

            self.applyTargets(ensemble, W);
        }

        pub fn update(self: *@This(), H: Matrix(T), U: Matrix(T)) void {
            switch (self.method) {
                inline else => |*field| field.update(H, U),
            }
        }

        fn applyTargets(self: *@This(), ensemble: *Ensemble(T), W: Matrix(T)) void {
            for (0..ensemble.s.length()) |i| {
                const c = ensemble.s.at(i);

                if (self.adia and self.targets[i] != c) {
                    const E_new = W.at(i, self.targets[i]);

                    if (rescaleMomentumIsotropic(T, ensemble, i, E_new - W.at(i, c))) {
                        ensemble.s.ptr(i).* = self.targets[i];
                    }
                }

                if (!self.adia and self.targets[i] != c) {
                    ensemble.s.ptr(i).* = self.targets[i];
                }
            }
        }

        fn calcProbs(self: *@This(), probs: *Matrix(T), ensemble: *Ensemble(T), nstep: usize, dt: T) !void {
            self.probs.zero();

            if (self.adia) try self.calcProbsAdia(probs, ensemble, nstep, dt);
            if (!self.adia) try self.calcProbsDia(probs, ensemble, nstep, dt);
        }

        fn calcProbsAdia(self: *@This(), probs: *Matrix(T), ensemble: *Ensemble(T), nstep: usize, dt: T) !void {
            switch (self.method) {
                inline else => |*field| try field.calcProbsAdia(probs, ensemble, nstep, dt),
            }
        }

        fn calcProbsDia(self: *@This(), probs: *Matrix(T), ensemble: *Ensemble(T), nstep: usize, dt: T) !void {
            switch (self.method) {
                inline else => |*field| try field.calcProbsDia(probs, ensemble, nstep, dt),
            }
        }

        fn calcTargetStates(self: *@This(), ensemble: *Ensemble(T)) void {
            for (0..ensemble.s.length()) |i| {
                if (self.targets[i] != ensemble.s.at(i)) continue;

                var sum: T = 0;
                var accum: T = 0;

                for (0..self.probs.ncol()) |j| {
                    sum += self.probs.at(i, j);
                }

                if (sum > 1) for (0..self.probs.ncol()) |j| {
                    self.probs.ptr(i, j).* /= sum;
                };

                const rv = self.rng.random().float(T);

                for (0..self.probs.ncol()) |j| {
                    accum += self.probs.at(i, j);

                    if (rv < accum) {
                        self.targets[i] = j;
                        break;
                    }
                }
            }
        }
    };
}

// SPECIFIC METHODS ====================================================================================================

pub fn FewestSwitches(comptime T: type) type {
    return struct {
        coef: Matrix(Complex(T)),
        itg: Integrator(Complex(T)),
        ham: Matrix(T),
        sigma: Matrix(T),
        uhist: [2]Matrix(T),

        pub fn init(opt: anytype, nstate: usize, ntraj: usize, istate: usize, adia: bool, gpa: Allocator) !@This() {
            const cols = if (adia) nstate else nstate * nstate;

            var ham = try Matrix(T).init(ntraj, cols, gpa);
            errdefer ham.deinit(gpa);

            var uhist: [2]Matrix(T) = undefined;

            uhist[0] = try Matrix(T).init(ntraj, nstate * nstate, gpa);
            errdefer uhist[0].deinit(gpa);

            uhist[1] = try Matrix(T).init(ntraj, nstate * nstate, gpa);
            errdefer uhist[1].deinit(gpa);

            uhist[0].fill(std.math.nan(T));
            uhist[1].fill(std.math.nan(T));

            var coef = try Matrix(Complex(T)).init(ntraj, nstate, gpa);
            errdefer coef.deinit(gpa);

            var sigma = try Matrix(T).init(ntraj, nstate * nstate, gpa);
            errdefer sigma.deinit(gpa);

            coef.zero();

            for (0..ntraj) |i| {
                coef.ptr(i, istate).* = Complex(T).init(1, 0);
            }

            const itg_tag = switch (opt.integrator) {
                inline else => |tag| @field(Integrator(Complex(T)).Method, @tagName(tag)),
            };

            const itg = try Integrator(Complex(T)).init(itg_tag, nstate, gpa);

            return .{ .coef = coef, .ham = ham, .uhist = uhist, .sigma = sigma, .itg = itg };
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.coef.deinit(gpa);
            self.ham.deinit(gpa);
            self.sigma.deinit(gpa);

            self.uhist[0].deinit(gpa);
            self.uhist[1].deinit(gpa);

            self.itg.deinit(gpa);
        }

        pub fn calcProbsAdia(self: *@This(), probs: *Matrix(T), ensemble: *Ensemble(T), nstep: usize, dt: T) !void {
            if (std.math.isNan(self.uhist[0].at(0, 0))) return;

            std.debug.assert(probs.nrow() == ensemble.s.length());
            std.debug.assert(probs.ncol() == ensemble.nstate);
            std.debug.assert(self.sigma.nrow() == ensemble.s.length());
            std.debug.assert(self.sigma.ncol() == ensemble.nstate * ensemble.nstate);

            for (0..ensemble.s.length()) |i| {
                const c = ensemble.s.at(i);

                const row_sigma = self.sigma.rowSlice(i);

                const row_u_new = self.uhist[1].rowSlice(i);
                const row_u_old = self.uhist[0].rowSlice(i);

                hammesSchifferTully(T, row_sigma, row_u_new, row_u_old, dt * @as(T, @floatFromInt(nstep)));

                const int_ctx = .{ .ham = self.ham.rowSlice(i), .sigma = row_sigma, .nstate = ensemble.nstate };

                self.itg.step(self.coef.rowSlice(i), dt, int_ctx, coefDerAdia);

                const rho_cc = self.coef.at(i, c).squaredMagnitude();

                if (rho_cc < 1e-14) continue;

                for (0..ensemble.nstate) |j| {
                    if (c == j) continue;

                    const coef_c = self.coef.at(i, c);
                    const coef_j = self.coef.at(i, j);

                    const flux = self.sigma.at(i, c * ensemble.nstate + j) * coef_c.conjugate().mul(coef_j).re;

                    probs.ptr(i, j).* += @max(0, 2 * dt * flux / rho_cc);
                }
            }
        }

        pub fn calcProbsDia(self: *@This(), probs: *Matrix(T), ensemble: *Ensemble(T), _: usize, dt: T) !void {
            std.debug.assert(probs.nrow() == ensemble.s.length());
            std.debug.assert(probs.ncol() == ensemble.nstate);
            std.debug.assert(self.ham.nrow() == ensemble.s.length());
            std.debug.assert(self.ham.ncol() == ensemble.nstate * ensemble.nstate);

            for (0..ensemble.s.length()) |i| {
                const c = ensemble.s.at(i);

                const int_ctx = .{ .ham = self.ham.rowSlice(i), .nstate = ensemble.nstate };

                self.itg.step(self.coef.rowSlice(i), dt, int_ctx, coefDerDia);

                const rho_cc = self.coef.at(i, c).squaredMagnitude();

                if (rho_cc < 1e-14) continue;

                for (0..ensemble.nstate) |j| {
                    if (c == j) continue;

                    const coef_c = self.coef.at(i, c);
                    const coef_j = self.coef.at(i, j);

                    const flux = self.ham.at(i, c * ensemble.nstate + j) * coef_c.mul(coef_j.conjugate()).im;

                    probs.ptr(i, j).* += @max(0, 2 * dt * flux / rho_cc);
                }
            }
        }

        pub fn update(self: *@This(), H: Matrix(T), U: Matrix(T)) void {
            for (0..H.nrow()) |i| for (0..H.ncol()) |j| {
                self.ham.ptr(i, j).* = H.at(i, j);
            };

            for (0..self.uhist[0].nrow()) |i| for (0..self.uhist[0].ncol()) |j| {
                self.uhist[0].ptr(i, j).* = self.uhist[1].at(i, j);

                self.uhist[1].ptr(i, j).* = U.at(i, j);
            };

            if (!std.math.isNan(self.uhist[0].at(0, 0))) {
                const nstate = self.coef.ncol();

                for (0..self.uhist[0].nrow()) |i| for (0..nstate) |j| {
                    var overlap: T = 0;

                    for (0..nstate) |k| {
                        overlap += self.uhist[1].at(i, k * nstate + j) * self.uhist[0].at(i, k * nstate + j);
                    }

                    if (overlap < 0) for (0..nstate) |k| {
                        self.uhist[1].ptr(i, k * nstate + j).* = -self.uhist[1].at(i, k * nstate + j);
                    };
                };
            }
        }

        fn coefDerAdia(ctx: anytype, y: []const Complex(T), dy: []Complex(T)) void {
            const nstate = ctx.nstate;

            for (0..nstate) |i| {
                var sum_sigma = Complex(T).init(0, 0);

                for (0..nstate) |j| {
                    const sig = Complex(T).init(ctx.sigma[i * nstate + j], 0);

                    sum_sigma = sum_sigma.add(sig.mul(y[j]));
                }

                dy[i] = Complex(T).init(0, -ctx.ham[i]).mul(y[i]).sub(sum_sigma);
            }
        }

        fn coefDerDia(ctx: anytype, y: []const Complex(T), dy: []Complex(T)) void {
            const nstate = ctx.nstate;

            for (0..nstate) |i| {
                var sum = Complex(T).init(0, 0);

                for (0..nstate) |j| {
                    const term = y[j].mul(Complex(T).init(0, -ctx.ham[i * nstate + j]));

                    sum = sum.add(term);
                }

                dy[i] = sum;
            }
        }
    };
}

pub fn LandauZener(comptime T: type) type {
    return struct {
        history: [3]Matrix(T),

        pub fn init(_: anytype, nstate: usize, ntraj: usize, _: usize, adia: bool, gpa: Allocator) !@This() {
            var history: [3]Matrix(T) = undefined;

            const cols = if (adia) nstate else nstate * nstate;

            history[0] = try Matrix(T).init(ntraj, cols, gpa);
            errdefer history[0].deinit(gpa);

            history[1] = try Matrix(T).init(ntraj, cols, gpa);
            errdefer history[1].deinit(gpa);

            history[2] = try Matrix(T).init(ntraj, cols, gpa);

            history[0].fill(std.math.nan(T));
            history[1].fill(std.math.nan(T));
            history[2].fill(std.math.nan(T));

            return .{ .history = history };
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.history[0].deinit(gpa);
            self.history[1].deinit(gpa);
            self.history[2].deinit(gpa);
        }

        pub fn calcProbsAdia(self: *@This(), probs: *Matrix(T), ensemble: *Ensemble(T), _: usize, dt: T) !void {
            if (std.math.isNan(self.history[0].at(0, 0))) return;

            std.debug.assert(probs.nrow() == ensemble.s.length());
            std.debug.assert(probs.ncol() == ensemble.nstate);
            std.debug.assert(self.history[0].nrow() == ensemble.s.length());
            std.debug.assert(self.history[0].ncol() == ensemble.nstate);

            for (0..self.history[0].nrow()) |i| {
                const c = ensemble.s.at(i);

                for (0..self.history[0].ncol()) |j| {
                    if (c == j) continue;

                    const z0 = @abs(self.history[0].at(i, j) - self.history[0].at(i, c));
                    const z1 = @abs(self.history[1].at(i, j) - self.history[1].at(i, c));
                    const z2 = @abs(self.history[2].at(i, j) - self.history[2].at(i, c));

                    const dz1 = (z1 - z0) / dt;
                    const dz2 = (z2 - z1) / dt;

                    if (dz1 < 0 and dz2 >= 0) {
                        const d2z = (dz2 - dz1) / dt;

                        if (d2z > 1e-14) {
                            const b = (z2 - z0) / (2 * dt);
                            const a = d2z / 2;

                            var z_min = z1 - (b * b) / (4 * a);

                            if (z_min < 0) z_min = 0;

                            const p = std.math.exp(-(std.math.pi / 2.0) * std.math.sqrt((z_min * z_min * z_min) / d2z));

                            probs.ptr(i, j).* = p;
                        }
                    }
                }
            }
        }

        pub fn calcProbsDia(self: *@This(), probs: *Matrix(T), ensemble: *Ensemble(T), _: usize, dt: T) !void {
            if (std.math.isNan(self.history[1].at(0, 0))) return;

            std.debug.assert(probs.nrow() == ensemble.s.length());
            std.debug.assert(probs.ncol() == ensemble.nstate);
            std.debug.assert(self.history[0].nrow() == ensemble.s.length());
            std.debug.assert(self.history[0].ncol() == ensemble.nstate * ensemble.nstate);

            for (0..self.history[0].nrow()) |i| {
                const c = ensemble.s.at(i);

                for (0..ensemble.nstate) |j| {
                    if (c == j) continue;

                    const v_cc_old = self.history[1].at(i, c * ensemble.nstate + c);
                    const v_cc_new = self.history[2].at(i, c * ensemble.nstate + c);
                    const v_jj_old = self.history[1].at(i, j * ensemble.nstate + j);
                    const v_jj_new = self.history[2].at(i, j * ensemble.nstate + j);

                    const gap_old = v_cc_old - v_jj_old;
                    const gap_new = v_cc_new - v_jj_new;

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

        pub fn update(self: *@This(), H: Matrix(T), _: Matrix(T)) void {
            for (0..H.nrow()) |i| for (0..H.ncol()) |j| {
                self.history[0].ptr(i, j).* = self.history[1].at(i, j);
                self.history[1].ptr(i, j).* = self.history[2].at(i, j);

                self.history[2].ptr(i, j).* = H.at(i, j);
            };
        }
    };
}

// TIME DERIVATIVE COUPLINGS ===========================================================================================

pub fn hammesSchifferTully(comptime T: type, sigma: []T, U_new: []const T, U_old: []const T, dt: T) void {
    const nstate = std.math.sqrt(sigma.len);

    for (0..nstate) |i| for (0..nstate) |j| {
        var S_ij: T = 0;
        var S_ji: T = 0;

        for (0..nstate) |k| {
            const U_t_ki = U_new[k * nstate + i];
            const U_t_kj = U_new[k * nstate + j];

            const U_old_ki = U_old[k * nstate + i];
            const U_old_kj = U_old[k * nstate + j];

            S_ij += U_t_ki * U_old_kj;
            S_ji += U_t_kj * U_old_ki;
        }

        sigma[i * nstate + j] = (S_ji - S_ij) / (2 * dt);
    };
}

// HELPER FUNCTIONS ====================================================================================================

fn rescaleMomentumIsotropic(comptime T: type, ensemble: *Ensemble(T), i: usize, dE: T) bool {
    var ekin_old: T = 0;

    for (0..ensemble.p.ncol()) |j| {
        ekin_old += ensemble.p.at(i, j) * ensemble.p.at(i, j) / (2 * ensemble.mass);
    }

    if (ekin_old < 1e-14) {
        return false;
    }

    const ekin_new = ekin_old - dE;

    if (ekin_new >= 0) {
        for (0..ensemble.p.ncol()) |j| {
            ensemble.p.ptr(i, j).* *= std.math.sqrt(ekin_new / ekin_old);
        }

        return true;
    }

    return false;
}
