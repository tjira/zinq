//! Computes quantum dynamics scattering probabilities and cross sections using time-to-energy Fourier transform flux analysis.

const std = @import("std");

const fftw = @import("cimport.zig").fftw;

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;

const Grid = @import("wavepacket.zig").Grid;
const Matrix = @import("tensor.zig").Matrix;
const FftPlan = @import("fourier_transform.zig").FftPlan;
const Potential = @import("potential.zig").Potential;

const eighSlice = @import("linear_algebra.zig").eighSlice;
const writeMatrixLspace = @import("read_write.zig").writeMatrixLspace;

/// Holds all parameter configuration, input matrices, and grid references required for computing transition probabilities.
pub fn FluxAnalysisContext(comptime T: type) type {
    return struct {
        mass: T,

        fsurf: T,
        e_min: T,
        e_max: T,

        e_step: T,
        Vreact: T,

        istate: usize,
        grid: Grid(T),

        dt: T,

        /// Extracts necessary parameters and computes reactant potential energy directly during initialization.
        pub fn init(opt: anytype, grid: Grid(T), pot: Potential(T), gpa: Allocator) !@This() {
            const flux_opt = opt.flux_analysis.?;

            const V_arr = try gpa.alloc(T, pot.nstate() * pot.nstate());
            defer gpa.free(V_arr);

            pot.eval(T, V_arr, opt.initial_conditions.position, 0);

            var Vreact: T = 0;

            if (opt.initial_conditions.adiabatic) {
                const U_arr = try gpa.alloc(T, pot.nstate() * pot.nstate());
                defer gpa.free(U_arr);

                const W_arr = try gpa.alloc(T, pot.nstate());
                defer gpa.free(W_arr);

                try eighSlice(T, W_arr, U_arr, V_arr);

                Vreact = W_arr[opt.initial_conditions.state];
            }

            if (!opt.initial_conditions.adiabatic) {
                Vreact = V_arr[opt.initial_conditions.state * pot.nstate() + opt.initial_conditions.state];
            }

            return .{
                .fsurf = flux_opt.flux_surface,
                .e_min = flux_opt.e_min,
                .e_max = flux_opt.e_max,
                .e_step = flux_opt.e_step,
                .mass = opt.mass[0],
                .dt = opt.time_step,
                .istate = opt.initial_conditions.state,
                .Vreact = Vreact,
                .grid = grid,
            };
        }

        /// Computes transition probabilities by integrating flux of energy-resolved wavefunctions at a dividing surface.
        pub fn analyze(self: @This(), wfn_init: Matrix(Complex(T)), flux_acc: Matrix(Complex(T)), gpa: Allocator) !Matrix(T) {
            const dx, const dt, const mass = .{ self.grid.dr, self.dt, self.mass };

            const j0: usize = @intFromFloat(@round((self.fsurf - self.grid.r.at(0, 0)) / dx));

            if (j0 == 0 or j0 >= self.grid.r.nrow() - 1) {
                std.log.err("FLUX SURFACE IS OUTSIDE THE GRID INTERIOR", .{});

                return error.InvalidInput;
            }

            const ne = @as(usize, @intFromFloat(@round((self.e_max - self.e_min) / self.e_step))) + 1;

            var prob_matrix = try Matrix(T).initZero(ne, wfn_init.nrow(), gpa);
            errdefer prob_matrix.deinit(gpa);

            var temp_phi = try gpa.alloc(Complex(T), self.grid.r.nrow());
            defer gpa.free(temp_phi);

            const shape = try gpa.alloc(i32, 1);
            defer gpa.free(shape);

            shape[0] = @as(i32, @intCast(self.grid.r.nrow()));

            const ffft_plan = try FftPlan(Complex(T)).init(temp_phi, shape, -1, fftw.FFTW_ESTIMATE);
            defer ffft_plan.deinit();

            const ifft_plan = try FftPlan(Complex(T)).init(temp_phi, shape, 1, fftw.FFTW_ESTIMATE);
            defer ifft_plan.deinit();

            for (0..ne) |ei| {
                const E = self.e_min + @as(T, @floatFromInt(ei)) * self.e_step;

                if (E <= self.Vreact) {
                    continue;
                }

                const k_inc, var a_k = .{ std.math.sqrt(2 * self.mass * (E - self.Vreact)), Complex(T).init(0, 0) };

                for (0..self.grid.r.nrow()) |j| {
                    const phase = std.math.complex.exp(Complex(T).init(0, -k_inc * self.grid.r.at(j, 0)));

                    a_k = a_k.add(phase.mul(wfn_init.at(self.istate, j)));
                }

                if (a_k.squaredMagnitude() < 1e-6) {
                    continue;
                }

                for (0..wfn_init.nrow()) |f| {
                    const row = ei * wfn_init.nrow() + f;
                    const phi_mid = flux_acc.at(row, j0);

                    for (0..self.grid.r.nrow()) |j| {
                        temp_phi[j] = flux_acc.at(row, j);
                    }

                    ffft_plan.execute(temp_phi);

                    for (0..self.grid.r.nrow()) |m| {
                        temp_phi[m] = temp_phi[m].mul(Complex(T).init(0, self.grid.k.at(m, 0)));
                    }

                    ifft_plan.execute(temp_phi);

                    const factor = (k_inc * dt * dt) / (mass * mass * dx * dx * a_k.squaredMagnitude());

                    prob_matrix.ptr(ei, f).* = factor * phi_mid.conjugate().mul(temp_phi[j0]).im;
                }
            }

            return prob_matrix;
        }
    };
}
