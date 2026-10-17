//! Ehrenfest dynamics propagation of electronic coefficients along classical trajectories.

const std = @import("std");

const Complex = std.math.Complex;

const CouplingOptions = @import("surface_hopping.zig").CouplingOptions;
const Ensemble = @import("classical_dynamics.zig").Ensemble;
const Integrator = @import("integrator.zig").Integrator;
const Matrix = @import("tensor.zig").Matrix;

const checkNorm = @import("surface_hopping.zig").checkNorm;
const coefDerAdia = @import("surface_hopping.zig").coefDerAdia;
const coefDerDia = @import("surface_hopping.zig").coefDerDia;
const hammesSchifferTully = @import("surface_hopping.zig").hammesSchifferTully;
const nacvCoupling = @import("surface_hopping.zig").nacvCoupling;

/// Configuration options for Ehrenfest molecular dynamics.
pub const Options = struct {
    coupling: ?CouplingOptions = null,
    integrator: std.meta.Tag(Integrator(f64).Method) = .rk4,
    norm_tolerance: ?f64 = 0.05,
    nstep: u32 = 10,
};

/// Manages the electronic wavefunction coefficients and TDSE propagation for Ehrenfest dynamics.
pub fn Ehrenfest(comptime T: type) type {
    return struct {
        coefics: Matrix(Complex(T)),
        itg: Integrator(Complex(T)),

        ham: Matrix(T),
        norm_toler: ?T,
        nosteps: usize,

        sigmatdc: Matrix(T),
        uhist: [2]Matrix(T),

        coupling: ?std.meta.Tag(CouplingOptions),

        /// Allocates and initializes electronic coefficients, Hamiltonian histories, and Runge-Kutta integrator.
        pub fn init(opt: Options, nstate: usize, ntraj: usize, adia: bool, gpa: std.mem.Allocator) !@This() {
            if (adia and opt.coupling == null) {
                return error.InvalidInput;
            }

            if (!adia and opt.coupling != null) {
                return error.InvalidInput;
            }

            var coef = try Matrix(Complex(T)).init(ntraj, nstate, gpa);
            errdefer coef.deinit(gpa);

            coef.zero();

            const cols = if (adia) nstate else nstate * nstate;

            var ham = try Matrix(T).init(ntraj, cols, gpa);
            errdefer ham.deinit(gpa);

            const coupling: ?std.meta.Tag(CouplingOptions) = if (opt.coupling) |c| switch (c) {
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

            const tag = switch (opt.integrator) {
                inline else => |t| @field(std.meta.Tag(Integrator(Complex(T)).Method), @tagName(t)),
            };

            const itg = try Integrator(Complex(T)).init(tag, nstate, gpa);
            errdefer itg.deinit(gpa);

            return .{
                .coefics = coef,
                .coupling = coupling,
                .ham = ham,
                .itg = itg,
                .norm_toler = if (opt.norm_tolerance) |tol| @floatCast(tol) else null,
                .nosteps = opt.nstep,
                .sigmatdc = sigmatdc,
                .uhist = uhist,
            };
        }

        /// Deallocates coefficients, Hamiltonians, coupling histories, and integrator workspaces.
        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            self.coefics.deinit(gpa);

            self.ham.deinit(gpa);
            self.itg.deinit(gpa);

            self.sigmatdc.deinit(gpa);
            self.uhist[0].deinit(gpa);
            self.uhist[1].deinit(gpa);
        }

        /// Returns electronic wavefunction coefficients for nonadiabatic Ehrenfest dynamics.
        pub fn coefs(self: @This()) Matrix(Complex(T)) {
            return self.coefics;
        }

        /// Initializes the active electronic state amplitude with unit norm.
        pub fn setInitialState(self: *@This(), istate: usize) void {
            self.coefics.zero();

            for (0..self.coefics.nrow()) |i| {
                self.coefics.ptr(i, istate).* = Complex(T).init(1, 0);
            }
        }

        /// Propagates electronic coefficients over nuclear time step dt in adiabatic or diabatic representation.
        pub fn step(self: *@This(), ensemble: *Ensemble(T), V: Matrix(T), W: Matrix(T), U: Matrix(T), grad_V: Matrix(T), is_td: bool, dt: T) !void {
            if (self.coupling != null) {
                self.update(W, U);

                try self.stepAdia(ensemble, grad_V, is_td, dt);
            }

            if (self.coupling == null) {
                try self.stepDia(V, dt);
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

        /// Integrates electronic coefficients in the adiabatic basis using time-derivative couplings.
        fn stepAdia(self: *@This(), ensemble: *Ensemble(T), grad_V: Matrix(T), is_td: bool, dt: T) !void {
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
        fn stepDia(self: *@This(), V: Matrix(T), dt: T) !void {
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
    };
}
