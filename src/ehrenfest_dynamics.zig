//! Ehrenfest dynamics propagation of electronic coefficients along classical trajectories.

const std = @import("std");

const Complex = std.math.Complex;

const CouplingOptions = @import("electronic_propagation.zig").CouplingOptions;
const ElectronicPropagator = @import("electronic_propagation.zig").ElectronicPropagator;
const Ensemble = @import("ensemble.zig").Ensemble;
const Integrator = @import("integrator.zig").Integrator;
const Matrix = @import("tensor.zig").Matrix;

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
        prop: ElectronicPropagator(T),

        /// Allocates and initializes electronic coefficients, Hamiltonian histories, and Runge-Kutta integrator.
        pub fn init(opt: Options, nstate: usize, ntraj: usize, istate: usize, adia: bool, gpa: std.mem.Allocator) !@This() {
            const prop = try ElectronicPropagator(T).init(opt.coupling, opt.integrator, opt.norm_tolerance, opt.nstep, nstate, ntraj, istate, adia, gpa);

            return .{ .coefics = prop.coefics, .prop = prop };
        }

        /// Deallocates coefficients, Hamiltonians, coupling histories, and integrator workspaces.
        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            self.prop.deinit(gpa);
        }

        /// Returns electronic wavefunction coefficients for nonadiabatic Ehrenfest dynamics.
        pub fn coefs(self: @This()) Matrix(Complex(T)) {
            return self.prop.coefs();
        }

        /// Propagates electronic coefficients over nuclear time step dt in adiabatic or diabatic representation.
        pub fn step(self: *@This(), ensemble: *Ensemble(T), V: Matrix(T), W: Matrix(T), U: Matrix(T), grad_V: Matrix(T), is_td: bool, dt: T) !void {
            if (self.prop.coupling != null) {
                self.prop.update(W, U);

                try self.prop.stepAdia(ensemble, grad_V, is_td, dt);
            }

            if (self.prop.coupling == null) {
                try self.prop.stepDia(V, dt);
            }
        }

        /// Updates Hamiltonian history and eigenvector history for coupling evaluation.
        pub fn update(self: *@This(), H: Matrix(T), U: Matrix(T)) void {
            self.prop.update(H, U);
        }
    };
}
