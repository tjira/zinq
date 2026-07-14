//! Ehrenfest dynamics propagation of electronic coefficients along classical trajectories.

const std = @import("std");

const Complex = std.math.Complex;

const Integrator = @import("integrator.zig").Integrator;
const Matrix = @import("tensor.zig").Matrix;

/// Configuration options for Ehrenfest molecular dynamics.
pub const Options = struct {
    integrator: std.meta.Tag(Integrator(f64).Method) = .rk4,
    nstep: u32 = 10,
};

/// Manages the electronic wavefunction coefficients and TDSE propagation for Ehrenfest dynamics.
pub fn Ehrenfest(comptime T: type) type {
    return struct {
        coefics: Matrix(Complex(T)),
        itg: Integrator(Complex(T)),

        ham: Matrix(T),
        nosteps: usize,

        /// Initializes the electronic coefficient and Hamiltonian matrices.
        pub fn init(opt: Options, nstate: usize, ntraj: usize, gpa: std.mem.Allocator) !@This() {
            var coef = try Matrix(Complex(T)).init(ntraj, nstate, gpa);
            errdefer coef.deinit(gpa);

            var ham = try Matrix(T).init(ntraj, nstate * nstate, gpa);
            errdefer ham.deinit(gpa);

            coef.zero();

            const tag = switch (opt.integrator) {
                inline else => |t| @field(std.meta.Tag(Integrator(Complex(T)).Method), @tagName(t)),
            };

            const itg = try Integrator(Complex(T)).init(tag, nstate, gpa);
            errdefer itg.deinit(gpa);

            return .{ .coefics = coef, .itg = itg, .ham = ham, .nosteps = opt.nstep };
        }

        /// Deallocates coefficients, Hamiltonians, and internal integrator workspaces.
        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            self.coefics.deinit(gpa);

            self.ham.deinit(gpa);
            self.itg.deinit(gpa);
        }

        /// Sets the initial electronic state coefficients using eigenvectors or Kronecker deltas.
        pub fn setInitialState(self: *@This(), istate: usize, adia: bool, U: Matrix(T)) void {
            const nstate = self.coefics.ncol();

            for (0..self.coefics.nrow()) |i| {
                if (adia) for (0..nstate) |k| {
                    self.coefics.ptr(i, k).* = Complex(T).init(U.at(i, k * nstate + istate), 0);
                };

                if (!adia) {
                    self.coefics.ptr(i, istate).* = Complex(T).init(1, 0);
                }
            }
        }

        /// Propagates the electronic coefficients over a time step dt using a series of sub-steps.
        pub fn step(self: *@This(), H: Matrix(T), dt: T) !void {
            for (0..H.nrow()) |i| for (0..H.ncol()) |j| {
                self.ham.ptr(i, j).* = H.at(i, j);
            };

            const subdt = dt / @as(T, @floatFromInt(self.nosteps));

            for (0..self.nosteps) |_| for (0..self.coefics.nrow()) |i| {
                self.itg.step(self.coefics.rowSlice(i), subdt, .{ .ham = self.ham.rowSlice(i) }, coefDerDia);
            };
        }

        /// Computes the time derivative of electronic coefficients in the diabatic basis.
        fn coefDerDia(ctx: anytype, y: []const Complex(T), dy: []Complex(T)) void {
            const nstate = std.math.sqrt(ctx.ham.len);

            for (0..nstate) |i| {
                var sum = Complex(T).init(0, 0);

                for (0..nstate) |j| {
                    sum = sum.add(y[j].mul(Complex(T).init(0, -ctx.ham[i * nstate + j])));
                }

                dy[i] = sum;
            }
        }
    };
}
