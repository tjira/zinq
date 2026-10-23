//! Computes quantum dynamics scattering probabilities and cross sections using time-to-energy Fourier transform flux analysis.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;

const Grid = @import("wavepacket.zig").Grid;
const Matrix = @import("tensor.zig").Matrix;
const Potential = @import("potential.zig").Potential;
const Vector = @import("tensor.zig").Vector;

const fft = @import("fourier_transform.zig").fft;
const fftn = @import("fourier_transform.zig").fftn;
const syevd = @import("linear_algebra.zig").syevd;

/// Holds all parameter configuration, input matrices, and grid references required for computing transition probabilities.
pub fn FluxAnalysis(comptime T: type) type {
    return struct {
        mass: []const T,
        gmma: []const T,

        e_min: T,
        e_max: T,
        initk: T,

        weight: T,

        flux_bounds: []const [2]T,

        npoint: usize,
        logarithmic: bool,

        dt: T,

        /// Initializes flux analysis parameters and allocates momentum and boundary vectors.
        pub fn init(opt: anytype, gpa: Allocator) !@This() {
            const flux_opt = opt.flux_analysis.?;

            const mass = try gpa.dupe(T, opt.mass);
            errdefer gpa.free(mass);

            const gammas = try gpa.dupe(T, opt.initial_conditions.gamma);
            errdefer gpa.free(gammas);

            const f_bounds = try gpa.dupe([2]T, flux_opt.flux_bounds);
            errdefer gpa.free(f_bounds);

            var r_perp: T = 1;

            for (1..gammas.len) |i| if (!(opt.grid.cylindrical and i == gammas.len - 1)) {
                r_perp *= std.math.sqrt(gammas[i] / std.math.pi);
            };

            const weight = if (opt.grid.cylindrical) (std.math.pi / gammas[gammas.len - 1]) / r_perp else 1 / r_perp;

            var ctx: @This() = undefined;

            ctx.e_min = @floatCast(flux_opt.e_min);
            ctx.e_max = @floatCast(flux_opt.e_max);

            ctx.npoint, ctx.logarithmic = .{ flux_opt.npoint, flux_opt.logarithmic };

            ctx.dt = @floatCast(opt.time_step);

            ctx.gmma, ctx.weight, ctx.mass, ctx.flux_bounds = .{ gammas, weight, mass, f_bounds };

            ctx.initk = @abs(opt.initial_conditions.momentum[0]);

            return ctx;
        }

        /// Deallocates the dynamically allocated slices in the context.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            gpa.free(self.mass);
            gpa.free(self.gmma);

            gpa.free(self.flux_bounds);
        }

        /// Computes the collision energy for the given index.
        pub fn getCollisionEnergy(opt: anytype, i: usize) T {
            if (opt.npoint <= 1) return @floatCast(opt.e_min);

            const fi = @as(T, @floatFromInt(i));

            const fn_minus_1 = @as(T, @floatFromInt(opt.npoint - 1));

            if (opt.logarithmic) {
                const log_min = std.math.log10(@as(T, @floatCast(opt.e_min)));
                const log_max = std.math.log10(@as(T, @floatCast(opt.e_max)));

                return std.math.pow(T, 10, log_min + fi * (log_max - log_min) / fn_minus_1);
            }

            const e_min: T, const e_max: T = .{ @floatCast(opt.e_min), @floatCast(opt.e_max) };

            return e_min + fi * (e_max - e_min) / fn_minus_1;
        }

        /// Computes the reactant asymptotic threshold energy (potential energy plus transversal zero-point energy).
        pub fn getThreshold(opt: anytype, pot: Potential(T), gpa: Allocator) !T {
            var V = try Matrix(T).init(pot.nstate(), pot.nstate(), gpa);
            defer V.deinit(gpa);

            try pot.eval(T, V.data, opt.initial_conditions.position, 0);

            if (opt.j_quantum_number > 0) {
                const j = @as(T, @floatFromInt(opt.j_quantum_number));

                const m, const r = .{ opt.mass[0], opt.initial_conditions.position[0] };

                for (0..pot.nstate()) |s| {
                    V.ptr(s, s).* += if (r != 0) j * (j + 1) / (2 * m * r * r) else 0;
                }
            }

            var Vreact: T = 0;

            if (opt.initial_conditions.adiabatic) {
                var U = try Matrix(T).init(pot.nstate(), pot.nstate(), gpa);
                defer U.deinit(gpa);

                var W = try Vector(T).init(pot.nstate(), gpa);
                defer W.deinit(gpa);

                try syevd(T, &W, &U, V);

                Vreact = W.at(opt.initial_conditions.state);
            }

            if (!opt.initial_conditions.adiabatic) {
                Vreact = V.at(opt.initial_conditions.state, opt.initial_conditions.state);
            }

            var e_perp: T = 0;

            for (1..opt.mass.len) |i| {
                const f_zpe = if (opt.grid.cylindrical and i == opt.mass.len - 1) @as(T, 1) else @as(T, 0.5);

                e_perp += f_zpe * opt.initial_conditions.gamma[i] / opt.mass[i];
            }

            return Vreact + e_perp;
        }

        /// Computes transition probabilities by integrating flux of energy-resolved wavefunctions at a dividing surface.
        pub fn run(self: @This(), grid: Grid(T), wfn_init: Matrix(Complex(T)), flux_acc: Matrix(Complex(T)), gpa: Allocator) !Matrix(T) {
            var npoint: usize = 1;

            while (try std.math.powi(usize, npoint, grid.ncol()) != grid.nrow()) {
                npoint += 1;
            }

            const ne = self.npoint;

            var sigma = try Matrix(T).initZero(ne, wfn_init.nrow(), gpa);
            errdefer sigma.deinit(gpa);

            var temp_phi = try gpa.alloc(Complex(T), grid.nrow());
            defer gpa.free(temp_phi);

            const shape = try gpa.alloc(usize, grid.ncol());
            defer gpa.free(shape);

            @memset(shape, npoint);

            const scratch = try gpa.alloc(Complex(T), if (grid.ncol() > 1) npoint else 0);
            defer gpa.free(scratch);

            for (0..grid.ncol()) |d| {
                const s_d = std.math.pow(usize, npoint, grid.ncol() - 1 - d);

                const n_min_f = (self.flux_bounds[d][0] - grid.getR(0, d)) / (grid.getR(s_d, d) - grid.getR(0, d));
                const n_max_f = (self.flux_bounds[d][1] - grid.getR(0, d)) / (grid.getR(s_d, d) - grid.getR(0, d));

                if (n_min_f < -0.5 or n_max_f >= @as(T, @floatFromInt(npoint)) - 0.5) {
                    std.log.err("FLUX BOUNDS MUST LIE WITHIN THE GRID BOUNDS", .{});

                    return error.InvalidInput;
                }

                const n_min: usize = @intFromFloat(@round(n_min_f));
                const n_max: usize = @intFromFloat(@round(n_max_f));

                const dx_d, const m = .{ grid.getR(s_d, d) - grid.getR(0, d), self.mass[0] };

                for (0..ne) |ei| {
                    const E = self.getCollisionEnergy(ei);

                    if (E <= 0) {
                        continue;
                    }

                    const k_inc = std.math.sqrt(2 * m * E);

                    const exp_arg = -std.math.pow(T, k_inc - self.initk, @as(T, 2)) / (self.gmma[0]);
                    const ak = std.math.sqrt(4 * std.math.pi / self.gmma[0]) * std.math.exp(exp_arg);

                    if (ak / (dx_d * dx_d) < 1e-3) {
                        continue;
                    }

                    for (0..wfn_init.nrow()) |f| {
                        const row, var sum: T = .{ ei * wfn_init.nrow() + f, 0 };

                        for (0..grid.nrow()) |j| {
                            temp_phi[j] = flux_acc.at(row, j);
                        }

                        if (grid.ncol() == 1) {
                            fft(T, temp_phi, temp_phi, -1);
                        }

                        if (grid.ncol() > 1) {
                            fftn(T, temp_phi, temp_phi, shape, -1, scratch);
                        }

                        for (0..grid.nrow()) |m_idx| {
                            temp_phi[m_idx] = temp_phi[m_idx].mul(Complex(T).init(0, grid.getK(m_idx, d)));
                        }

                        if (grid.ncol() == 1) {
                            fft(T, temp_phi, temp_phi, 1);
                        }

                        if (grid.ncol() > 1) {
                            fftn(T, temp_phi, temp_phi, shape, 1, scratch);
                        }

                        const factor = grid.dr * k_inc * self.dt * self.dt / (m * self.mass[d] * dx_d * ak);

                        for (0..grid.nrow()) |i| {
                            var in_bounds = true;

                            for (0..grid.ncol()) |k| {
                                if (k == d) continue;

                                const s_k = std.math.pow(usize, npoint, grid.ncol() - 1 - k);

                                const n_min_k_f = (self.flux_bounds[k][0] - grid.getR(0, k)) / (grid.getR(s_k, k) - grid.getR(0, k));
                                const n_max_k_f = (self.flux_bounds[k][1] - grid.getR(0, k)) / (grid.getR(s_k, k) - grid.getR(0, k));

                                const n_min_k: usize = @intFromFloat(@round(n_min_k_f));
                                const n_max_k: usize = @intFromFloat(@round(n_max_k_f));

                                if ((i / s_k) % npoint < n_min_k or (i / s_k) % npoint > n_max_k) {
                                    in_bounds = false;

                                    break;
                                }
                            }

                            if (!in_bounds) continue;

                            const val = flux_acc.at(row, i).conjugate().mul(temp_phi[i]).im;

                            if ((i / s_d) % npoint == n_max) sum += factor * val * self.weight;
                            if ((i / s_d) % npoint == n_min) sum -= factor * val * self.weight;
                        }

                        sigma.ptr(ei, f).* += sum;
                    }
                }
            }

            return sigma;
        }
    };
}
