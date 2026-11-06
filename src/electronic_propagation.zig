//! Solves the time-dependent Schrödinger equation for electronic coefficients along classical trajectories.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;

const Ensemble = @import("ensemble.zig").Ensemble;
const Integrator = @import("integrator.zig").Integrator;
const Matrix = @import("tensor.zig").Matrix;

/// Tagged union selecting the time-derivative nonadiabatic coupling scheme.
pub const CouplingOptions = union(enum) {
    hst: HstOptions,
    nacv: NacvOptions,
};

/// Options for Hammes-Schiffer-Tully antisymmetric wavefunction overlap coupling.
pub const HstOptions = struct {};

/// Options for nonadiabatic coupling vectors time-derivative coupling.
pub const NacvOptions = struct {};

/// Manages electronic coefficients, Hamiltonians, and Runge-Kutta propagation along classical trajectories.
pub fn ElectronicPropagator(comptime T: type) type {
    return struct {
        coefics: Matrix(Complex(T)),
        itg: Integrator(Complex(T)),

        ham: Matrix(T),
        norm_toler: ?T,
        nosteps: usize,

        sigmatdc: Matrix(T),
        uhist: [2]Matrix(T),

        coupling: ?std.meta.Tag(CouplingOptions),

        /// Allocates electronic coefficients, Hamiltonian histories, and Runge-Kutta workspaces.
        pub fn init(coupling_opt: ?CouplingOptions, itg_tag: std.meta.Tag(Integrator(f64).Method), norm_tol: ?f64, nstep: u32, nstate: usize, ntraj: usize, istate: usize, adia: bool, gpa: Allocator) !@This() {
            if (adia and coupling_opt == null) {
                return error.InvalidInput;
            }

            if (!adia and coupling_opt != null) {
                return error.InvalidInput;
            }

            var coef = try Matrix(Complex(T)).init(ntraj, nstate, gpa);
            errdefer coef.deinit(gpa);

            coef.zero();

            for (0..ntraj) |i| {
                coef.ptr(i, istate).* = Complex(T).init(1, 0);
            }

            const cols = if (adia) nstate else nstate * nstate;

            var ham = try Matrix(T).init(ntraj, cols, gpa);
            errdefer ham.deinit(gpa);

            const coupling: ?std.meta.Tag(CouplingOptions) = if (coupling_opt) |c| switch (c) {
                inline else => |_, tag| tag,
            } else null;

            var uhist: [2]Matrix(T) = undefined;

            uhist[0] = try Matrix(T).init(ntraj, nstate * nstate, gpa);
            errdefer uhist[0].deinit(gpa);

            uhist[1] = try Matrix(T).init(ntraj, nstate * nstate, gpa);
            errdefer uhist[1].deinit(gpa);

            uhist[0].fill(std.math.nan(T));
            uhist[1].fill(std.math.nan(T));

            var sigmatdc = try Matrix(T).init(ntraj, nstate * nstate, gpa);
            errdefer sigmatdc.deinit(gpa);

            const tag = switch (itg_tag) {
                inline else => |t| @field(std.meta.Tag(Integrator(Complex(T)).Method), @tagName(t)),
            };

            const itg = try Integrator(Complex(T)).init(tag, nstate, gpa);
            errdefer itg.deinit(gpa);

            return .{
                .coefics = coef,
                .coupling = coupling,
                .ham = ham,
                .itg = itg,
                .norm_toler = if (norm_tol) |tol| @floatCast(tol) else null,
                .nosteps = nstep,
                .sigmatdc = sigmatdc,
                .uhist = uhist,
            };
        }

        /// Deallocates coefficients, Hamiltonians, coupling histories, and integrator workspaces.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.coefics.deinit(gpa);

            self.ham.deinit(gpa);
            self.itg.deinit(gpa);

            self.sigmatdc.deinit(gpa);
            self.uhist[0].deinit(gpa);
            self.uhist[1].deinit(gpa);
        }

        /// Returns electronic wavefunction coefficients for nonadiabatic dynamics.
        pub fn coefs(self: @This()) Matrix(Complex(T)) {
            return self.coefics;
        }

        /// Integrates electronic coefficients in the adiabatic basis using time-derivative couplings.
        pub fn stepAdia(self: *@This(), ensemble: *Ensemble(T), grad_V: Matrix(T), is_td: bool, dt: T) !void {
            if (self.coupling.? == .hst and std.math.isNan(self.uhist[0].at(0, 0))) return;

            const subdt = dt / @as(T, @floatFromInt(self.nosteps));

            for (0..ensemble.s.length()) |i| {
                const row_sigma = self.sigmatdc.rowSlice(i);
                const row_u_new = self.uhist[1].rowSlice(i);
                const row_u_old = self.uhist[0].rowSlice(i);

                const ham = self.ham.rowSlice(i);

                const p, const m = .{ ensemble.p.rowSlice(i), ensemble.m };

                switch (self.coupling.?) {
                    .hst => hammesSchifferTully(T, row_sigma, row_u_new, row_u_old, dt),
                    .nacv => nacvCoupling(T, row_sigma, grad_V.rowSlice(i), p, m),
                }

                const int_ctx = .{ .ham = ham, .is_td = is_td, .sigma = row_sigma };

                for (0..self.nosteps) |_| {
                    self.itg.step(self.coefics.rowSlice(i), subdt, int_ctx, coefDerAdia);
                }

                try checkNorm(T, self.coefics.rowSlice(i), self.norm_toler);
            }
        }

        /// Integrates electronic coefficients in the diabatic basis using the diabatic potential matrix.
        pub fn stepDia(self: *@This(), V: Matrix(T), dt: T) !void {
            for (0..V.nrow()) |i| for (0..V.ncol()) |j| {
                self.ham.ptr(i, j).* = V.at(i, j);
            };

            const subdt = dt / @as(T, @floatFromInt(self.nosteps));

            for (0..self.nosteps) |_| for (0..self.coefics.nrow()) |i| {
                self.itg.step(self.coefics.rowSlice(i), subdt, .{ .ham = self.ham.rowSlice(i) }, coefDerDia);
            };

            for (0..self.coefics.nrow()) |i| {
                try checkNorm(T, self.coefics.rowSlice(i), self.norm_toler);
            }
        }

        /// Updates Hamiltonian history and eigenvector history for coupling evaluation.
        pub fn update(self: *@This(), H: Matrix(T), U: Matrix(T)) void {
            for (0..H.nrow()) |i| for (0..H.ncol()) |j| {
                self.ham.ptr(i, j).* = H.at(i, j);
            };

            if (self.coupling == null) return;

            for (0..self.uhist[0].nrow()) |i| for (0..self.uhist[0].ncol()) |j| {
                self.uhist[0].ptr(i, j).* = self.uhist[1].at(i, j);

                self.uhist[1].ptr(i, j).* = U.at(i, j);
            };
        }
    };
}

/// Verifies electronic wavefunction norm conservation within specified tolerance.
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

/// Evaluates adiabatic electronic coefficient time-derivatives including nonadiabatic couplings.
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

/// Evaluates diabatic electronic coefficient time-derivatives from diabatic Hamiltonian matrix elements.
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

/// Computes time-derivative couplings via Hammes-Schiffer-Tully wave function overlap approximation.
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

/// Computes time-derivative nonadiabatic couplings from nonadiabatic coupling vectors and nuclear velocities.
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
