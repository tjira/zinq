//! Trajectory surface hopping algorithms for nonadiabatic molecular dynamics transitions.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;
const DefaultPrng = std.Random.DefaultPrng;

const electronic_propagation = @import("electronic_propagation.zig");

const Ensemble = @import("ensemble.zig").Ensemble;
const Integrator = @import("integrator.zig").Integrator;
const Matrix = @import("tensor.zig").Matrix;

pub const CouplingOptions = electronic_propagation.CouplingOptions;
pub const HstOptions = electronic_propagation.HstOptions;
pub const NacvOptions = electronic_propagation.NacvOptions;
pub const checkNorm = electronic_propagation.checkNorm;
pub const coefDerAdia = electronic_propagation.coefDerAdia;
pub const coefDerDia = electronic_propagation.coefDerDia;
pub const hammesSchifferTully = electronic_propagation.hammesSchifferTully;
pub const nacvCoupling = electronic_propagation.nacvCoupling;

/// Options for selecting Fewest Switches Tully hopping or Landau-Zener hopping models.
pub const Options = union(enum) {
    fewest_switches: FewestSwitchesOptions,
    landau_zener: LandauZenerOptions,
};

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

        /// Returns electronic wavefunction coefficients for nonadiabatic dynamics if active.
        pub fn coefs(self: @This()) ?Matrix(Complex(T)) {
            switch (self.method) {
                .fewest_switches => |fs| return fs.coefics,
                inline else => return null,
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

            const any_hop = self.applyTargets(ensemble, W, grad_V);

            if (self.adia_alg) {
                self.applyDecoherence(ensemble, dt);
            }

            return any_hop;
        }

        /// Updates Hamiltonian history and eigenvector history for transition probability evaluation.
        pub fn update(self: *@This(), H: Matrix(T), U: Matrix(T)) void {
            switch (self.method) {
                inline else => |*field| field.update(H, U),
            }
        }

        /// Applies energy-based decoherence damping to non-active electronic state amplitudes.
        fn applyDecoherence(self: *@This(), ensemble: *Ensemble(T), dt: T) void {
            switch (self.method) {
                .fewest_switches => |*fs| fs.applyDecoherence(ensemble, dt),
                inline else => {},
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

                self.targets[i] = sampleTargetState(T, self.probs.rowSlice(i), ensemble.s.at(i), self.rng.random().float(T));
            }
        }
    };
}

/// Implements Tully's fewest switches surface hopping algorithm for nonadiabatic dynamics.
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
                if (!adia or std.math.isNan(a) or a <= 0) return error.InvalidInput;
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
                inline else => |tag| @field(std.meta.Tag(Integrator(Complex(T)).Method), @tagName(tag)),
            };

            const itg = try Integrator(Complex(T)).init(itg_tag, nstate, gpa);
            errdefer itg.deinit(gpa);

            return .{
                .alpha = if (opt.alpha) |a| @floatCast(a) else null,
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

        /// Damps nonadiabatic electronic state amplitudes according to the Granucci-Persico decoherence time.
        pub fn applyDecoherence(self: *@This(), ensemble: *Ensemble(T), dt: T) void {
            const alpha = self.alpha orelse return;

            for (0..ensemble.s.length()) |traj| {
                var ekin: T = 0;

                for (0..ensemble.p.ncol()) |j| {
                    ekin += ensemble.p.at(traj, j) * ensemble.p.at(traj, j) / (2 * ensemble.m[j]);
                }

                if (ekin <= 1e-4) continue;

                dampStateAmplitudes(T, self.coefics.rowSlice(traj), self.hamilton.rowSlice(traj), ensemble.s.at(traj), ekin, alpha, dt);
            }
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

/// Evaluates Granucci-Persico decoherence time from active-target energy gap and kinetic energy.
fn computeDecoherenceTime(comptime T: type, delta_e: T, ekin: T, alpha: T) T {
    return (1 + alpha / ekin) / delta_e;
}

/// Damps non-active electronic state amplitudes and renormalizes active state population.
fn dampStateAmplitudes(comptime T: type, coefs: []Complex(T), ham: []const T, c: usize, ekin: T, alpha: T, dt: T) void {
    if (ekin <= 1e-4) return;

    const rho_cc = coefs[c].squaredMagnitude();

    if (rho_cc < 1e-14) return;

    const e_curr = ham[c];

    var delta_rho: T = 0;

    for (0..coefs.len) |j| {
        if (j == c) continue;

        const delta_e = @abs(ham[j] - e_curr);

        if (delta_e < 1e-14) continue;

        const tau = computeDecoherenceTime(T, delta_e, ekin, alpha);

        const d = std.math.exp(-0.5 * dt / tau);

        const rho_j = coefs[j].squaredMagnitude();

        delta_rho += rho_j * (1 - d * d);

        coefs[j].re *= d;
        coefs[j].im *= d;
    }

    const renorm = std.math.sqrt((rho_cc + delta_rho) / rho_cc);

    coefs[c].re *= renorm;
    coefs[c].im *= renorm;
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

    return rescaleMomentumDirectional(T, ensemble, traj, c, k, a, b, dE, grad_V, nstate, reflect);
}

/// Applies directional momentum rescaling along nonadiabatic coupling vector or reflects on frustrated hops.
fn rescaleMomentumDirectional(comptime T: type, ensemble: *Ensemble(T), traj: usize, c: usize, k: usize, a: T, b: T, dE: T, grad_V: Matrix(T), nstate: usize, reflect: bool) bool {
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

    if (std.math.isNan(ekin_old) or ekin_old < 1e-14 or std.math.isNan(dE)) {
        return false;
    }

    const ekin_new = ekin_old - dE;

    if (std.math.isNan(ekin_new) or ekin_new < 1e-14) {
        return false;
    }

    const factor = std.math.sqrt(ekin_new / ekin_old);

    if (std.math.isNan(factor)) {
        return false;
    }

    for (0..ensemble.p.ncol()) |j| {
        ensemble.p.ptr(i, j).* *= factor;
    }

    return true;
}

/// Stochastically samples the target electronic state from accumulated transition probabilities.
fn sampleTargetState(comptime T: type, probs_row: []T, current: usize, rv: T) usize {
    var sum: T = 0;

    for (probs_row) |p| sum += p;

    if (sum > 1.0) {
        for (probs_row) |*p| p.* /= sum;
    }

    var accum: T = 0;

    for (probs_row, 0..) |p, j| {
        accum += p;

        if (rv < accum) return j;
    }

    return current;
}
