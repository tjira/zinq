//! Computes potential energy gradients, nonadiabatic couplings, and nuclear forces for trajectory ensembles.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;

const Ensemble = @import("ensemble.zig").Ensemble;
const Matrix = @import("tensor.zig").Matrix;
const Potential = @import("potential.zig").Potential;
const ScalarDual = @import("dual.zig").ScalarDual;
const Vector = @import("tensor.zig").Vector;

const syevdBatch = @import("linear_algebra.zig").syevdBatch;

/// Helper struct managing memory for potential energy gradients and wavefunctions.
pub fn GradientBuffer(comptime T: type) type {
    return struct {
        adia: bool,

        r_dual: Matrix(ScalarDual(T)),
        V_dual: Matrix(ScalarDual(T)),

        V: Matrix(T),
        W: Matrix(T),
        U: Matrix(T),

        U_prev: Matrix(T),
        grad_V: Matrix(T),

        /// Allocates memory for dual numbers and gradients needed in force evaluations.
        pub fn init(ndim: usize, nstate: usize, ntraj: usize, adia: bool, gpa: Allocator) !@This() {
            var r_dual = try Matrix(ScalarDual(T)).init(ntraj, ndim, gpa);
            errdefer r_dual.deinit(gpa);

            var V_dual = try Matrix(ScalarDual(T)).init(ntraj, nstate * nstate, gpa);
            errdefer V_dual.deinit(gpa);

            var V = try Matrix(T).init(ntraj, nstate * nstate, gpa);
            errdefer V.deinit(gpa);

            var W = try Matrix(T).init(ntraj, nstate, gpa);
            errdefer W.deinit(gpa);

            var U = try Matrix(T).init(ntraj, nstate * nstate, gpa);
            errdefer U.deinit(gpa);

            var U_prev = try Matrix(T).init(ntraj, nstate * nstate, gpa);
            errdefer U_prev.deinit(gpa);

            U_prev.fill(std.math.nan(T));

            const grad_V = try Matrix(T).init(ntraj, ndim * nstate * nstate, gpa);
            errdefer grad_V.deinit(gpa);

            return .{
                .r_dual = r_dual,
                .V_dual = V_dual,
                .U_prev = U_prev,
                .grad_V = grad_V,

                .V = V,
                .W = W,
                .U = U,

                .adia = adia,
            };
        }

        /// Deallocates gradient matrices, eigenvalues, and eigenvectors.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.r_dual.deinit(gpa);
            self.V_dual.deinit(gpa);

            self.V.deinit(gpa);
            self.W.deinit(gpa);
            self.U.deinit(gpa);

            self.U_prev.deinit(gpa);
            self.grad_V.deinit(gpa);
        }

        /// Computes and applies nuclear forces to trajectories based on potential gradients.
        pub fn apply(self: *@This(), ensemble: *Ensemble(T), pot: Potential(T), coefs: ?*const Matrix(Complex(T))) void {
            const nstate = pot.nstate();

            for (0..ensemble.r.nrow()) |i| for (0..ensemble.r.ncol()) |j| {
                if (coefs) |coef| {
                    if (self.adia) {
                        applyAdiabaticEhrenfestForces(T, ensemble, self.grad_V, self.W, coef, i, j, nstate);
                    }

                    if (!self.adia) {
                        applyDiabaticEhrenfestForces(T, ensemble, self.grad_V, coef, i, j, nstate);
                    }
                }

                if (coefs == null and self.adia) {
                    applyAdiabaticSingleForces(T, ensemble, self.grad_V, i, j, nstate);
                }

                if (coefs == null and !self.adia) {
                    applyDiabaticSingleForces(T, ensemble, self.grad_V, i, j, nstate);
                }
            };
        }

        /// Calculates average potential energy from cached gradient buffer matrices.
        pub fn epot(self: @This(), ensemble: Ensemble(T), nstate: usize, coefs: ?Matrix(Complex(T))) T {
            var sum: T = 0;

            if (coefs) |coef| for (0..ensemble.r.nrow()) |i| {
                var traj_epot: T = 0;

                if (self.adia) {
                    for (0..nstate) |k| {
                        const rho_kk = coef.at(i, k).squaredMagnitude();

                        traj_epot += rho_kk * self.W.at(i, k);
                    }
                }

                if (!self.adia) {
                    for (0..nstate) |k| for (0..nstate) |l| {
                        const rho_kl = coef.at(i, k).conjugate().mul(coef.at(i, l)).re;

                        traj_epot += rho_kl * self.V.at(i, k * nstate + l);
                    };
                }

                sum += traj_epot;
            };

            if (coefs == null and self.adia) for (0..ensemble.r.nrow()) |i| {
                sum += self.W.at(i, ensemble.s.at(i));
            };

            if (coefs == null and !self.adia) for (0..ensemble.r.nrow()) |i| {
                sum += self.V.at(i, ensemble.s.at(i) * nstate + ensemble.s.at(i));
            };

            return sum / @as(T, @floatFromInt(ensemble.r.nrow()));
        }

        /// Computes the ensemble-averaged potential energy for each electronic state.
        pub fn epotState(self: @This(), nstate: usize, gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(nstate, gpa);

            for (0..self.V.nrow()) |i| for (0..nstate) |k| {
                const v_k = if (self.adia) self.W.at(i, k) else self.V.at(i, k * nstate + k);

                value.ptr(k).* += v_k;
            };

            value.divs(@floatFromInt(self.V.nrow()));

            return value;
        }

        /// Evaluates the potential and its gradients using automatic differentiation or direct ab initio evaluation.
        pub fn update(self: *@This(), r: Matrix(T), pot: Potential(T), time: T, states: ?[]const usize) !void {
            switch (pot) {
                .ab_initio => |ab| {
                    try ab.updateGradients(&self.V, &self.grad_V, r, states);

                    if (self.adia) for (0..self.V.nrow()) |i| for (0..pot.nstate()) |k| {
                        self.W.ptr(i, k).* = self.V.at(i, k * pot.nstate() + k);

                        for (0..pot.nstate()) |l| {
                            self.U.ptr(i, k * pot.nstate() + l).* = if (k == l) 1 else 0;
                        }
                    };
                },

                inline else => {
                    try pot.evalBatch(T, &self.V, r, time);

                    if (self.adia) {
                        try syevdBatch(T, &self.W, &self.U, self.V);

                        alignEigenvectorPhases(T, &self.U, &self.U_prev, pot.nstate());
                    }

                    try self.evalDualGradients(r, pot, time);
                },
            }
        }

        /// Evaluates potential gradients across all coordinates using automatic differentiation.
        fn evalDualGradients(self: *@This(), r: Matrix(T), pot: Potential(T), time: T) !void {
            const nstate = pot.nstate();

            for (0..r.nrow()) |i| for (0..r.ncol()) |j| {
                for (0..r.ncol()) |k| {
                    self.r_dual.ptr(i, k).* = ScalarDual(T).init(r.at(i, k), if (k == j) 1 else 0);
                }

                const v_row, const r_row = .{ self.V_dual.rowSlice(i), self.r_dual.rowSlice(i) };

                try pot.eval(ScalarDual(T), v_row, r_row, ScalarDual(T).init(time, 0));

                if (self.adia) {
                    transformAdiabaticGradients(T, &self.grad_V, self.V_dual, self.W, self.U, i, j, nstate);
                }

                if (!self.adia) {
                    for (0..nstate * nstate) |k| {
                        self.grad_V.ptr(i, j * nstate * nstate + k).* = self.V_dual.at(i, k).der;
                    }
                }
            };
        }
    };
}

/// Aligns eigenvector column signs between consecutive steps to ensure phase continuity.
fn alignEigenvectorPhases(comptime T: type, U: *Matrix(T), U_prev: *Matrix(T), nstate: usize) void {
    if (std.math.isNan(U_prev.at(0, 0))) {
        @memcpy(U_prev.data, U.data);

        return;
    }

    for (0..U.nrow()) |i| for (0..nstate) |k| {
        var overlap: T = 0;

        for (0..nstate) |a| {
            const u_k = U.at(i, a * nstate + k);

            overlap += u_k * U_prev.at(i, a * nstate + k);
        }

        if (overlap < 0) for (0..nstate) |a| {
            U.ptr(i, a * nstate + k).* = -U.at(i, a * nstate + k);
        };
    };

    @memcpy(U_prev.data, U.data);
}

/// Applies adiabatic mean-field Ehrenfest forces to trajectories.
fn applyAdiabaticEhrenfestForces(comptime T: type, ensemble: *Ensemble(T), grad_V: Matrix(T), W: Matrix(T), coef: *const Matrix(Complex(T)), i: usize, j: usize, nstate: usize) void {
    var der: T = 0;

    for (0..nstate) |m| for (0..nstate) |n| {
        const rho_mn = coef.at(i, m).conjugate().mul(coef.at(i, n)).re;

        if (m == n) {
            der += rho_mn * grad_V.at(i, j * nstate * nstate + m * nstate + n);
        }

        if (m != n) {
            const dE = W.at(i, n) - W.at(i, m);

            der += rho_mn * dE * grad_V.at(i, j * nstate * nstate + m * nstate + n);
        }
    };

    ensemble.a.ptr(i, j).* = -der / ensemble.m[j];
}

/// Applies single-state adiabatic forces to trajectories.
fn applyAdiabaticSingleForces(comptime T: type, ensemble: *Ensemble(T), grad_V: Matrix(T), i: usize, j: usize, nstate: usize) void {
    const s = ensemble.s.at(i);

    const der = grad_V.at(i, j * nstate * nstate + s * nstate + s);

    ensemble.a.ptr(i, j).* = -der / ensemble.m[j];
}

/// Applies diabatic mean-field Ehrenfest forces to trajectories.
fn applyDiabaticEhrenfestForces(comptime T: type, ensemble: *Ensemble(T), grad_V: Matrix(T), coef: *const Matrix(Complex(T)), i: usize, j: usize, nstate: usize) void {
    var der: T = 0;

    for (0..nstate) |k| for (0..nstate) |l| {
        const rho_kl = coef.at(i, k).conjugate().mul(coef.at(i, l)).re;

        der += rho_kl * grad_V.at(i, j * nstate * nstate + k * nstate + l);
    };

    ensemble.a.ptr(i, j).* = -der / ensemble.m[j];
}

/// Applies single-state diabatic forces to trajectories.
fn applyDiabaticSingleForces(comptime T: type, ensemble: *Ensemble(T), grad_V: Matrix(T), i: usize, j: usize, nstate: usize) void {
    const k = ensemble.s.at(i) * nstate + ensemble.s.at(i);

    ensemble.a.ptr(i, j).* = -grad_V.at(i, j * nstate * nstate + k) / ensemble.m[j];
}

/// Transforms diabatic dual gradients into adiabatic diagonal gradients and nonadiabatic derivative couplings.
fn transformAdiabaticGradients(comptime T: type, grad_V: *Matrix(T), V_dual: Matrix(ScalarDual(T)), W: Matrix(T), U: Matrix(T), i: usize, j: usize, nstate: usize) void {
    for (0..nstate) |k| for (0..nstate) |l| {
        var sum: T = 0;

        for (0..nstate) |a| for (0..nstate) |b| {
            const u_ak = U.at(i, a * nstate + k);
            const u_bl = U.at(i, b * nstate + l);

            const d_ab = V_dual.at(i, a * nstate + b).der;

            sum += u_ak * u_bl * d_ab;
        };

        if (k == l) {
            grad_V.ptr(i, j * nstate * nstate + k * nstate + l).* = sum;
        }

        if (k != l) {
            const dE = W.at(i, l) - W.at(i, k);

            const val = if (@abs(dE) < 1e-14) 0 else sum / dE;

            grad_V.ptr(i, j * nstate * nstate + k * nstate + l).* = val;
        }
    };
}
