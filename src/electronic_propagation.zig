//! Solves the time-dependent Schrödinger equation for electronic coefficients along classical trajectories.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;

const Ensemble = @import("ensemble.zig").Ensemble;
const Integrator = @import("integrator.zig").Integrator;
const Matrix = @import("tensor.zig").Matrix;

const logmOrthogonalWithWorkspace = @import("matrix_logarithm.zig").logmOrthogonalWithWorkspace;
const logmWorkspaceSize = @import("matrix_logarithm.zig").logmWorkspaceSize;

/// Tagged union selecting the time-derivative nonadiabatic coupling scheme.
pub const CouplingOptions = union(enum) {
    hst: HstOptions,
    nacv: NacvOptions,
    npi: NpiOptions,
};

/// Options for Hammes-Schiffer-Tully antisymmetric wavefunction overlap coupling.
pub const HstOptions = struct {};

/// Options for nonadiabatic coupling vectors time-derivative coupling.
pub const NacvOptions = struct {};

/// Options for norm-preserving interpolation overlap matrix logarithm coupling.
pub const NpiOptions = struct {};

/// Manages electronic coefficients, Hamiltonians, and Runge-Kutta propagation along classical trajectories.
pub fn ElectronicPropagator(comptime T: type) type {
    return struct {
        coefics: Matrix(Complex(T)),
        itg: Integrator(Complex(T)),

        ham: Matrix(T),
        norm_toler: ?T,
        nosteps: usize,

        npi_work: ?[]T,

        sigmatdc: Matrix(T),
        uhist: [2]Matrix(T),

        coupling: ?std.meta.Tag(CouplingOptions),

        /// Allocates electronic coefficients, Hamiltonian histories, and Runge-Kutta workspaces.
        pub fn init(cpl: ?CouplingOptions, itg_m: anytype, tol: ?f64, nstep: u32, nstate: usize, ntraj: usize, istate: usize, adia: bool, gpa: Allocator) !@This() {
            if (adia and cpl == null) {
                return error.InvalidInput;
            }

            if (!adia and cpl != null) {
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

            const coupling: ?std.meta.Tag(CouplingOptions) = if (cpl) |c| switch (c) {
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

            const tag = switch (itg_m) {
                inline else => |t| @field(std.meta.Tag(Integrator(Complex(T)).Method), @tagName(t)),
            };

            var itg = try Integrator(Complex(T)).init(tag, nstate, gpa);
            errdefer itg.deinit(gpa);

            var npi_work: ?[]T = null;

            if (coupling != null and coupling.? == .npi) {
                npi_work = try gpa.alloc(T, npiWorkspaceSize(nstate));
            }

            errdefer if (npi_work) |w| gpa.free(w);

            return .{
                .norm_toler = if (tol) |t| @floatCast(t) else null,

                .coefics = coef,

                .ham = ham,
                .itg = itg,

                .nosteps = nstep,

                .coupling = coupling,
                .npi_work = npi_work,
                .sigmatdc = sigmatdc,

                .uhist = uhist,
            };
        }

        /// Deallocates coefficients, Hamiltonians, coupling histories, and integrator workspaces.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.coefics.deinit(gpa);

            self.ham.deinit(gpa);
            self.itg.deinit(gpa);

            if (self.npi_work) |w| gpa.free(w);
            self.sigmatdc.deinit(gpa);
            self.uhist[0].deinit(gpa);
            self.uhist[1].deinit(gpa);
        }

        /// Returns electronic wavefunction coefficients for nonadiabatic dynamics.
        pub fn coefs(self: @This()) Matrix(Complex(T)) {
            return self.coefics;
        }

        /// Evaluates time-derivative nonadiabatic couplings for a single trajectory.
        pub fn calcCoupling(self: *@This(), traj: usize, grad_V_row: []const T, p: []const T, m: []const T, dt: T) !void {
            const row_sigma = self.sigmatdc.rowSlice(traj);
            const row_u_new = self.uhist[1].rowSlice(traj);
            const row_u_old = self.uhist[0].rowSlice(traj);

            switch (self.coupling.?) {
                .hst => hammesSchifferTully(T, row_sigma, row_u_new, row_u_old, dt),
                .nacv => nacvCoupling(T, row_sigma, grad_V_row, p, m),
                .npi => try normPreservingInterpolation(T, row_sigma, row_u_new, row_u_old, dt, self.npi_work.?),
            }
        }

        /// Checks whether eigenvector history is available for overlap-based coupling methods.
        pub fn hasHistory(self: @This()) bool {
            if (self.coupling == null) return true;

            if ((self.coupling.? == .hst or self.coupling.? == .npi) and std.math.isNan(self.uhist[0].at(0, 0))) return false;

            return true;
        }

        /// Evaluates time-derivative nonadiabatic couplings for all trajectories.
        pub fn calcCouplings(self: *@This(), ensemble: *Ensemble(T), grad_V: Matrix(T), dt: T) !void {
            if (!self.hasHistory()) return;

            for (0..ensemble.s.length()) |i| {
                try self.calcCoupling(i, grad_V.rowSlice(i), ensemble.p.rowSlice(i), ensemble.m, dt);
            }
        }

        /// Integrates electronic coefficients for a single trajectory in the adiabatic basis.
        pub fn stepTrajectoryAdia(self: *@This(), traj: usize, is_td: bool, dt: T) !void {
            const int_ctx = .{ .ham = self.ham.rowSlice(traj), .is_td = is_td, .sigma = self.sigmatdc.rowSlice(traj) };

            self.itg.step(self.coefics.rowSlice(traj), dt, int_ctx, coefDerAdia);

            try checkNorm(T, self.coefics.rowSlice(traj), self.norm_toler);
        }

        /// Integrates electronic coefficients for a single trajectory in the diabatic basis.
        pub fn stepTrajectoryDia(self: *@This(), traj: usize, dt: T) !void {
            const int_ctx = .{ .ham = self.ham.rowSlice(traj) };

            self.itg.step(self.coefics.rowSlice(traj), dt, int_ctx, coefDerDia);

            try checkNorm(T, self.coefics.rowSlice(traj), self.norm_toler);
        }

        /// Integrates electronic coefficients in the adiabatic basis using time-derivative couplings.
        pub fn stepAdia(self: *@This(), ensemble: *Ensemble(T), grad_V: Matrix(T), is_td: bool, dt: T) !void {
            if (!self.hasHistory()) return;

            try self.calcCouplings(ensemble, grad_V, dt);

            const subdt = dt / @as(T, @floatFromInt(self.nosteps));

            for (0..ensemble.s.length()) |i| for (0..self.nosteps) |_| {
                try self.stepTrajectoryAdia(i, is_td, subdt);
            };
        }

        /// Integrates electronic coefficients in the diabatic basis using the diabatic potential matrix.
        pub fn stepDia(self: *@This(), V: Matrix(T), dt: T) !void {
            for (0..V.nrow()) |i| for (0..V.ncol()) |j| {
                self.ham.ptr(i, j).* = V.at(i, j);
            };

            const subdt = dt / @as(T, @floatFromInt(self.nosteps));

            for (0..self.coefics.nrow()) |i| for (0..self.nosteps) |_| {
                try self.stepTrajectoryDia(i, subdt);
            };
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

/// Computes the required workspace slice length for norm-preserving interpolation.
pub fn npiWorkspaceSize(nstate: usize) usize {
    return logmWorkspaceSize(nstate) + nstate * nstate;
}

/// Computes time-derivative couplings via norm-preserving interpolation using the overlap matrix logarithm.
pub fn normPreservingInterpolation(comptime T: type, sigma: []T, U_new: []const T, U_old: []const T, dt: T, work: []T) !void {
    const nstate = std.math.sqrt(sigma.len);

    if (nstate == 1) {
        sigma[0] = 0;

        return;
    }

    if (nstate == 2) {
        var p01 = U_old[0] * U_new[1] + U_old[2] * U_new[3];
        var p11 = U_old[1] * U_new[1] + U_old[3] * U_new[3];

        const p00 = U_old[0] * U_new[0] + U_old[2] * U_new[2];
        const p10 = U_old[1] * U_new[0] + U_old[3] * U_new[2];

        if (@abs(p01) > 0.5 and @abs(p10) > 0.5 and p01 * p10 > 0) {
            p01 = -p01;
            p11 = -p11;
        }

        const theta = std.math.atan2(p10 - p01, p00 + p11);

        sigma[1] = (0 - theta) / dt;
        sigma[2] = (0 + theta) / dt;

        sigma[0] = 0;
        sigma[3] = 0;

        return;
    }

    std.debug.assert(work.len >= npiWorkspaceSize(nstate));

    const P, const logm_work = .{ work[0 .. nstate * nstate], work[nstate * nstate ..] };

    for (0..nstate) |i| for (0..nstate) |j| {
        var sum: T = 0;

        for (0..nstate) |k| {
            sum += U_old[k * nstate + i] * U_new[k * nstate + j];
        }

        P[i * nstate + j] = sum;
    };

    for (0..nstate) |i| for (i + 1..nstate) |j| {
        const pij = P[i * nstate + j];
        const pji = P[j * nstate + i];

        if (@abs(pij) > 0.5 and @abs(pji) > 0.5 and pij * pji > 0) {
            for (0..nstate) |row| {
                P[row * nstate + j] = -P[row * nstate + j];
            }
        }
    };

    var L_mat = Matrix(T){ .data = sigma, .shape = .{ nstate, nstate } };

    const P_matrix = Matrix(T){ .data = P, .shape = .{ nstate, nstate } };

    try logmOrthogonalWithWorkspace(T, &L_mat, P_matrix, logm_work);

    const inv_dt = 1 / dt;

    for (sigma) |*elem| {
        elem.* *= inv_dt;
    }
}
