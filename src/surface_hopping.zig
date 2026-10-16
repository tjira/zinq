//! Trajectory surface hopping algorithms for nonadiabatic molecular dynamics transitions.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;
const DefaultPrng = std.Random.DefaultPrng;

const Ensemble = @import("classical_dynamics.zig").Ensemble;
const Integrator = @import("integrator.zig").Integrator;
const Matrix = @import("tensor.zig").Matrix;

/// Options for selecting Fewest Switches Tully hopping or Landau-Zener hopping models.
pub const Options = union(enum) {
    fewest_switches: FewestSwitchesOptions,
    landau_zener: LandauZenerOptions,
};

/// Tagged union selecting the time-derivative nonadiabatic coupling scheme.
pub const CouplingOptions = union(enum) {
    hst: HstOptions,
    nacv: NacvOptions,
};

/// Options for Hammes-Schiffer-Tully antisymmetric wavefunction overlap coupling.
pub const HstOptions = struct {};

/// Options for nonadiabatic coupling vectors time-derivative coupling.
pub const NacvOptions = struct {};

/// Options configuration for Tully's fewest switches surface hopping.
const FewestSwitchesOptions = struct {
    alpha: ?f64 = null,
    coupling: ?CouplingOptions = null,
    integrator: std.meta.Tag(Integrator(f64).Method) = .rk4,
    norm_tolerance: ?f64 = 0.05,
    nstep: u32 = 10,
    reflect_frustrated: bool = false,
    seed: u32 = 1,
};

/// Options configuration for Landau-Zener surface hopping.
const LandauZenerOptions = struct {
    seed: u32 = 1,
};

/// Implements nonadiabatic state transitions via surface hopping algorithms.
pub fn SurfaceHopping(comptime T: type) type {
    return struct {
        adia_alg: bool,
        method: Method,
        nosteps: usize,

        reflect_frustrated: bool,

        probs: Matrix(T),
        targets: []usize,
        rng: DefaultPrng,

        /// Union type representing the specific surface hopping method used.
        pub const Method = union(enum) {
            fewest_switches: FewestSwitches(T),
            landau_zener: LandauZener(T),
        };

        /// Initializes the surface hopping solver and the chosen transition probability method.
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

            const reflect_frustrated = switch (options) {
                .fewest_switches => |field| field.reflect_frustrated,
                inline else => false,
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
            errdefer probs.deinit(gpa);

            var sh: @This() = undefined;

            sh.adia_alg = adia;
            sh.method = method;
            sh.nosteps = nstep;

            sh.reflect_frustrated = reflect_frustrated;

            sh.rng, sh.probs, sh.targets = .{ rng, probs, targets };

            return sh;
        }

        /// Deallocates surface hopping state variables and method buffers.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.probs.deinit(gpa);
            gpa.free(self.targets);

            switch (self.method) {
                inline else => |*field| field.deinit(gpa),
            }
        }

        /// Executes a surface hopping step, calculating transition probabilities and updating states.
        pub fn hop(self: *@This(), ensemble: *Ensemble(T), V: Matrix(T), W: Matrix(T), U: Matrix(T), grad_V: Matrix(T), is_td: bool, dt: T) !bool {
            self.update(if (self.adia_alg) W else V, U);

            const subdt = dt / @as(T, @floatFromInt(self.nosteps));

            for (0..self.targets.len) |i| {
                self.targets[i] = ensemble.s.at(i);
            }

            for (0..self.nosteps) |step_idx| {
                try self.calcProbs(ensemble, grad_V, is_td, step_idx, subdt);

                self.calcTargetStates(ensemble);
            }

            return self.applyTargets(ensemble, W, grad_V);
        }

        /// Updates Hamiltonian and unitary transformation histories for hopping evaluation.
        pub fn update(self: *@This(), H: Matrix(T), U: Matrix(T)) void {
            switch (self.method) {
                inline else => |*field| field.update(H, U),
            }
        }

        /// Applies state transitions with momentum rescaling along coupling vectors or isotropically.
        fn applyTargets(self: *@This(), ensemble: *Ensemble(T), W: Matrix(T), grad_V: Matrix(T)) bool {
            var any_hop, const nstate = .{ false, W.ncol() };

            for (0..ensemble.s.length()) |i| {
                const c = ensemble.s.at(i);

                if (self.adia_alg and self.targets[i] != c) {
                    const target = self.targets[i];

                    const dE = W.at(i, target) - W.at(i, c);

                    const enough_energy = rescaleMomentum(T, ensemble, i, c, target, dE, grad_V, nstate, self.reflect_frustrated);

                    if (enough_energy) {
                        ensemble.s.ptr(i).* = target;

                        self.targets[i], any_hop = .{ c, true };
                    }

                    if (!enough_energy) {
                        self.targets[i] = c;
                    }
                }

                if (!self.adia_alg and self.targets[i] != c) {
                    ensemble.s.ptr(i).* = self.targets[i];

                    self.targets[i], any_hop = .{ c, true };
                }
            }

            return any_hop;
        }

        /// Computes state transition probabilities at the current dynamics step.
        fn calcProbs(self: *@This(), ensemble: *Ensemble(T), grad_V: Matrix(T), is_td: bool, step_idx: usize, dt: T) !void {
            self.probs.zero();

            if (self.adia_alg) {
                try self.calcProbsAdia(ensemble, grad_V, is_td, step_idx, dt);
            }

            if (!self.adia_alg) {
                try self.calcProbsDia(ensemble, dt);
            }
        }

        /// Computes adiabatic representation state transition probabilities.
        fn calcProbsAdia(self: *@This(), ensemble: *Ensemble(T), grad_V: Matrix(T), is_td: bool, step_idx: usize, dt: T) !void {
            switch (self.method) {
                inline else => |*field| {
                    try field.calcProbsAdia(&self.probs, ensemble, grad_V, is_td, step_idx, self.nosteps, dt);
                },
            }
        }

        /// Computes diabatic representation state transition probabilities.
        fn calcProbsDia(self: *@This(), ensemble: *Ensemble(T), dt: T) !void {
            switch (self.method) {
                inline else => |*field| try field.calcProbsDia(&self.probs, ensemble, self.nosteps, dt),
            }
        }

        /// Performs stochastic sampling to select potential target electronic states.
        fn calcTargetStates(self: *@This(), ensemble: *Ensemble(T)) void {
            for (0..ensemble.s.length()) |i| {
                if (self.targets[i] != ensemble.s.at(i)) continue;

                var sum: T = 0;

                for (0..self.probs.ncol()) |j| {
                    sum += self.probs.at(i, j);
                }

                var accum: T = 0;

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

/// Implements Tully's fewest switches surface hopping algorithm.
fn FewestSwitches(comptime T: type) type {
    return struct {
        coefics: Matrix(Complex(T)),
        itg: Integrator(Complex(T)),

        alpha: ?T,
        normt: ?T,

        hamilton: Matrix(T),
        sigmatdc: Matrix(T),
        uhist: [2]Matrix(T),

        coupling: ?std.meta.Tag(CouplingOptions),

        /// Allocates and initializes wavefunction coefficients and Tully integrator.
        pub fn init(opt: anytype, nstate: usize, ntraj: usize, istate: usize, adia: bool, gpa: Allocator) !@This() {
            if (opt.alpha) |a| {
                if (!adia or a < 0) return error.InvalidInput;
            }

            if (adia and opt.coupling == null) {
                return error.InvalidInput;
            }

            if (!adia and opt.coupling != null) {
                return error.InvalidInput;
            }

            const coupling: ?std.meta.Tag(CouplingOptions) = if (opt.coupling) |c| switch (c) {
                inline else => |_, tag| tag,
            } else null;

            if (opt.reflect_frustrated and (coupling == null or coupling.? != .nacv)) {
                return error.InvalidInput;
            }

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
            errdefer itg.deinit(gpa);

            return .{
                .alpha = opt.alpha,
                .coefics = coef,
                .coupling = coupling,
                .hamilton = ham,
                .itg = itg,
                .normt = if (opt.norm_tolerance) |tol| @floatCast(tol) else null,
                .sigmatdc = sigma,
                .uhist = uhist,
            };
        }

        /// Deallocates Tully fewest switches coefficients and integrator resources.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.coefics.deinit(gpa);

            self.hamilton.deinit(gpa);
            self.sigmatdc.deinit(gpa);

            self.uhist[0].deinit(gpa);
            self.uhist[1].deinit(gpa);

            self.itg.deinit(gpa);
        }

        /// Computes adiabatic transition probabilities using time-derivative couplings.
        pub fn calcProbsAdia(self: *@This(), probs: *Matrix(T), ensemble: *Ensemble(T), grad_V: Matrix(T), is_td: bool, step_idx: usize, nstep: usize, dt: T) !void {
            if (self.coupling.? == .hst and std.math.isNan(self.uhist[0].at(0, 0))) return;

            for (0..ensemble.s.length()) |i| {
                const c = ensemble.s.at(i);

                const row_sigma = self.sigmatdc.rowSlice(i);

                const row_u_new = self.uhist[1].rowSlice(i);
                const row_u_old = self.uhist[0].rowSlice(i);

                const eff_dt = dt * @as(T, @floatFromInt(nstep));

                const ham = self.hamilton.rowSlice(i);

                if (step_idx == 0) {
                    switch (self.coupling.?) {
                        .hst => hammesSchifferTully(T, row_sigma, row_u_new, row_u_old, eff_dt),
                        .nacv => nacvCoupling(T, row_sigma, grad_V.rowSlice(i), ensemble.p.rowSlice(i), ensemble.m),
                    }
                }

                const int_ctx = .{ .ham = ham, .is_td = is_td, .sigma = row_sigma };

                self.itg.step(self.coefics.rowSlice(i), dt, int_ctx, coefDerAdia);

                try checkNorm(T, self.coefics.rowSlice(i), self.normt);

                const rho_cc = self.coefics.at(i, c).squaredMagnitude();

                if (rho_cc >= 1e-14) {
                    for (0..probs.ncol()) |j| {
                        if (c == j) continue;

                        const coef_c = self.coefics.at(i, c);
                        const coef_j = self.coefics.at(i, j);

                        const flux = self.sigmatdc.at(i, c * probs.ncol() + j) * coef_c.conjugate().mul(coef_j).re;

                        probs.ptr(i, j).* += @max(0, 2 * dt * flux / rho_cc);
                    }
                }

                self.applyDecoherence(ensemble, i, dt);
            }
        }

        /// Computes diabatic transition probabilities using diabatic coupling elements.
        pub fn calcProbsDia(self: *@This(), probs: *Matrix(T), ensemble: *Ensemble(T), _: usize, dt: T) !void {
            for (0..ensemble.s.length()) |i| {
                const c = ensemble.s.at(i);

                const int_ctx = .{ .ham = self.hamilton.rowSlice(i) };

                self.itg.step(self.coefics.rowSlice(i), dt, int_ctx, coefDerDia);

                try checkNorm(T, self.coefics.rowSlice(i), self.normt);

                const rho_cc = self.coefics.at(i, c).squaredMagnitude();

                if (rho_cc < 1e-14) continue;

                for (0..probs.ncol()) |j| {
                    if (c == j) continue;

                    const coef_c = self.coefics.at(i, c);
                    const coef_j = self.coefics.at(i, j);

                    const flux = self.hamilton.at(i, c * probs.ncol() + j) * coef_c.mul(coef_j.conjugate()).im;

                    probs.ptr(i, j).* += @max(0, 2 * dt * flux / rho_cc);
                }
            }
        }

        /// Updates Hamiltonian history and eigenvector history for transition probability evaluation.
        pub fn update(self: *@This(), H: Matrix(T), U: Matrix(T)) void {
            for (0..H.nrow()) |i| for (0..H.ncol()) |j| {
                self.hamilton.ptr(i, j).* = H.at(i, j);
            };

            if (self.coupling == null) return;

            for (0..self.uhist[0].nrow()) |i| for (0..self.uhist[0].ncol()) |j| {
                self.uhist[0].ptr(i, j).* = self.uhist[1].at(i, j);

                self.uhist[1].ptr(i, j).* = U.at(i, j);
            };
        }

        /// Damps nonadiabatic electronic state amplitudes according to the Granucci-Persico decoherence time.
        fn applyDecoherence(self: *@This(), ensemble: *Ensemble(T), traj: usize, dt: T) void {
            const alpha = self.alpha orelse return;

            var ekin: T = 0;

            for (0..ensemble.p.ncol()) |j| {
                ekin += ensemble.p.at(traj, j) * ensemble.p.at(traj, j) / (2 * ensemble.m[j]);
            }

            if (ekin <= 1e-4) return;

            const c = ensemble.s.at(traj);

            const e_curr = self.hamilton.at(traj, c);

            for (0..self.coefics.ncol()) |j| {
                if (j == c) continue;

                const delta_e = @abs(self.hamilton.at(traj, j) - e_curr);

                if (delta_e < 1e-14) continue;

                const tau = (1 + alpha / ekin) / delta_e;
                const d = std.math.exp(-0.5 * dt / tau);

                self.coefics.ptr(traj, j).re *= d;
                self.coefics.ptr(traj, j).im *= d;
            }

            var sum_norm: T = 1;

            for (0..self.coefics.ncol()) |j| if (j != c) {
                sum_norm -= self.coefics.at(traj, j).squaredMagnitude();
            };

            const renorm = std.math.sqrt(@max(0, sum_norm) / (self.coefics.at(traj, c).squaredMagnitude() + 1e-14));

            self.coefics.ptr(traj, c).re *= renorm;
            self.coefics.ptr(traj, c).im *= renorm;
        }
    };
}

/// Implements the Landau-Zener model for nonadiabatic transition probabilities.
fn LandauZener(comptime T: type) type {
    return struct {
        history: [3]Matrix(T),

        /// Allocates history matrices for tracking energy levels.
        pub fn init(_: anytype, nstate: usize, ntraj: usize, _: usize, adia: bool, gpa: Allocator) !@This() {
            var history: [3]Matrix(T) = undefined;

            const cols = if (adia) nstate else nstate * nstate;

            history[0] = try Matrix(T).init(ntraj, cols, gpa);
            errdefer history[0].deinit(gpa);

            history[1] = try Matrix(T).init(ntraj, cols, gpa);
            errdefer history[1].deinit(gpa);

            history[2] = try Matrix(T).init(ntraj, cols, gpa);
            errdefer history[2].deinit(gpa);

            history[0].fill(std.math.nan(T));
            history[1].fill(std.math.nan(T));
            history[2].fill(std.math.nan(T));

            return .{ .history = history };
        }

        /// Deallocates history matrices used in Landau-Zener calculations.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.history[0].deinit(gpa);
            self.history[1].deinit(gpa);
            self.history[2].deinit(gpa);
        }

        /// Computes transition probabilities from adiabatic energy level curvatures.
        pub fn calcProbsAdia(self: *@This(), probs: *Matrix(T), ensemble: *Ensemble(T), _: Matrix(T), _: bool, _: usize, _: usize, dt: T) !void {
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

        /// Computes transition probabilities from diabatic energy level crossings.
        pub fn calcProbsDia(self: *@This(), probs: *Matrix(T), ensemble: *Ensemble(T), _: usize, dt: T) !void {
            if (std.math.isNan(self.history[1].at(0, 0))) return;

            for (0..self.history[0].nrow()) |i| {
                const c = ensemble.s.at(i);

                for (0..probs.ncol()) |j| {
                    if (c == j) continue;

                    const v_cc_old = self.history[1].at(i, c * probs.ncol() + c);
                    const v_cc_new = self.history[2].at(i, c * probs.ncol() + c);
                    const v_jj_old = self.history[1].at(i, j * probs.ncol() + j);
                    const v_jj_new = self.history[2].at(i, j * probs.ncol() + j);

                    const gap_old = v_cc_old - v_jj_old;
                    const gap_new = v_cc_new - v_jj_new;

                    if ((gap_old > 0 and gap_new <= 0) or (gap_old < 0 and gap_new >= 0)) {
                        const d_gap = @abs(gap_new - gap_old) / dt;

                        if (d_gap > 1e-14) {
                            const v_cj = @abs(self.history[2].at(i, c * probs.ncol() + j));

                            probs.ptr(i, j).* = 1 - std.math.exp(-2 * std.math.pi * v_cj * v_cj / d_gap);
                        }
                    }
                }
            }
        }

        /// Updates Hamiltonian history matrices.
        pub fn update(self: *@This(), H: Matrix(T), _: Matrix(T)) void {
            for (0..H.nrow()) |i| for (0..H.ncol()) |j| {
                self.history[0].ptr(i, j).* = self.history[1].at(i, j);
                self.history[1].ptr(i, j).* = self.history[2].at(i, j);

                self.history[2].ptr(i, j).* = H.at(i, j);
            };
        }
    };
}

/// Validates electronic state coefficient normalization against a tolerance threshold.
pub fn checkNorm(comptime T: type, coefs: []const Complex(T), tol: ?T) !void {
    const t = tol orelse return;

    var norm_val: T = 0;

    for (coefs) |c| {
        norm_val += c.squaredMagnitude();
    }

    if (std.math.isNan(norm_val) or @abs(norm_val - 1) > t) {
        return error.NormDivergence;
    }
}

/// Computes time derivative of electronic coefficients in the adiabatic basis.
pub fn coefDerAdia(ctx: anytype, y: anytype, dy: anytype) void {
    const e_ref = if (ctx.is_td) 0 else ctx.ham[0];

    const T = @TypeOf(y[0].re);

    for (0..ctx.ham.len) |i| {
        var sum_sigma = Complex(T).init(0, 0);

        for (0..ctx.ham.len) |j| {
            const sig = Complex(T).init(ctx.sigma[i * ctx.ham.len + j], 0);

            sum_sigma = sum_sigma.add(sig.mul(y[j]));
        }

        dy[i] = Complex(T).init(0, -(ctx.ham[i] - e_ref)).mul(y[i]).sub(sum_sigma);
    }
}

/// Computes time derivative of electronic coefficients in the diabatic basis.
pub fn coefDerDia(ctx: anytype, y: anytype, dy: anytype) void {
    const nstate = std.math.sqrt(ctx.ham.len);

    const T = @TypeOf(y[0].re);

    for (0..nstate) |i| {
        var sum = Complex(T).init(0, 0);

        for (0..nstate) |j| {
            const term = y[j].mul(Complex(T).init(0, -ctx.ham[i * nstate + j]));

            sum = sum.add(term);
        }

        dy[i] = sum;
    }
}

/// Calculates nonadiabatic coupling terms using Hammes-Schiffer-Tully approximation.
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

/// Evaluates nonadiabatic time-derivative coupling matrix elements from nonadiabatic coupling vectors.
pub fn nacvCoupling(comptime T: type, sigma: []T, grad_V: []const T, p: []const T, m: []const T) void {
    const nstate = std.math.sqrt(sigma.len);

    for (0..nstate) |i| for (0..nstate) |j| {
        if (i == j) {
            sigma[i * nstate + j] = 0;

            continue;
        }

        var v_dot_grad: T = 0;

        for (0..p.len) |alpha| {
            const v_alpha = p[alpha] / m[alpha];

            v_dot_grad += v_alpha * grad_V[alpha * nstate * nstate + i * nstate + j];
        }

        sigma[i * nstate + j] = v_dot_grad;
    };
}

/// Rescales nuclear momentum along nonadiabatic coupling vectors or isotropically to conserve total energy.
fn rescaleMomentum(comptime T: type, ensemble: *Ensemble(T), traj: usize, c: usize, k: usize, dE: T, grad_V: Matrix(T), nstate: usize, reflect: bool) bool {
    var a: T = 0;
    var b: T = 0;

    for (0..ensemble.p.ncol()) |j| {
        const d_j, const m_j = .{ grad_V.at(traj, j * nstate * nstate + c * nstate + k), ensemble.m[j] };

        a += 0.5 * d_j * d_j / m_j;

        b += (ensemble.p.at(traj, j) / m_j) * d_j;
    }

    if (a < 1e-14) {
        return rescaleMomentumIsotropic(T, ensemble, traj, dE);
    }

    const discriminant = b * b - 4 * a * dE;

    if (discriminant >= 0) {
        const sgn_b: T = (if (b < 0) -1 else 1);

        const gamma = (-b + sgn_b * std.math.sqrt(discriminant)) / (2 * a);

        for (0..ensemble.p.ncol()) |j| {
            const d_j = grad_V.at(traj, j * nstate * nstate + c * nstate + k);

            ensemble.p.ptr(traj, j).* += gamma * d_j;
        }

        return true;
    }

    if (reflect) {
        const factor = b / a;

        for (0..ensemble.p.ncol()) |j| {
            const d_j = grad_V.at(traj, j * nstate * nstate + c * nstate + k);

            ensemble.p.ptr(traj, j).* -= factor * d_j;
        }
    }

    return false;
}

/// Rescales nuclear momentum isotropically to conserve total energy after hopping.
fn rescaleMomentumIsotropic(comptime T: type, ensemble: *Ensemble(T), i: usize, dE: T) bool {
    var ekin_old: T = 0;

    for (0..ensemble.p.ncol()) |j| {
        ekin_old += ensemble.p.at(i, j) * ensemble.p.at(i, j) / (2 * ensemble.m[j]);
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
