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
pub fn FluxAnalysis(comptime T: type) type {
    return struct {
        mass: []const T,
        gmma: []const T,

        e_min: T,
        e_max: T,
        initk: T,

        e_thrs: T,
        e_step: T,
        weight: T,

        flux_bounds: []const [2]T,

        dt: T,

        /// Extracts necessary parameters and computes reactant potential energy directly during initialization.
        pub fn init(opt: anytype, pot: Potential(T), gpa: Allocator) !@This() {
            const flux_opt = opt.flux_analysis.?;

            const V_arr = try gpa.alloc(T, pot.nstate() * pot.nstate());
            defer gpa.free(V_arr);

            pot.eval(T, V_arr, opt.initial_conditions.position, 0);

            if (opt.j_quantum_number > 0) {
                const j = @as(T, @floatFromInt(opt.j_quantum_number));

                const m, const r = .{ opt.mass[0], opt.initial_conditions.position[0] };

                for (0..pot.nstate()) |s| {
                    V_arr[s * pot.nstate() + s] += if (r != 0) j * (j + 1) / (2 * m * r * r) else 0;
                }
            }

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

            const mass = try gpa.dupe(T, opt.mass);
            errdefer gpa.free(mass);

            const gammas = try gpa.dupe(T, opt.initial_conditions.gamma);
            errdefer gpa.free(gammas);

            const f_bounds = try gpa.dupe([2]T, flux_opt.flux_bounds);
            errdefer gpa.free(f_bounds);

            var e_perp: T, var r_perp: T = .{ 0, 1 };

            for (1..opt.mass.len) |i| {
                const f_zpe = if (opt.grid.cylindrical and i == opt.mass.len - 1) @as(T, 1) else @as(T, 0.5);

                e_perp += f_zpe * opt.initial_conditions.gamma[i] / mass[i];
            }

            for (1..gammas.len) |i| if (!(opt.grid.cylindrical and i == gammas.len - 1)) {
                r_perp *= std.math.sqrt(gammas[i] / std.math.pi);
            };

            const weight = if (opt.grid.cylindrical) (std.math.pi / gammas[gammas.len - 1]) / r_perp else 1 / r_perp;

            var ctx: @This() = undefined;

            ctx.e_min = flux_opt.e_min;
            ctx.e_max = flux_opt.e_max;

            ctx.e_step, ctx.dt, ctx.e_thrs = .{ flux_opt.e_step, opt.time_step, Vreact + e_perp };
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

        /// Computes transition probabilities by integrating flux of energy-resolved wavefunctions at a dividing surface.
        pub fn run(self: @This(), grid: Grid(T), wfn_init: Matrix(Complex(T)), flux_acc: Matrix(Complex(T)), gpa: Allocator) !Matrix(T) {
            var npoint: usize = 1;

            while (try std.math.powi(usize, npoint, grid.ncol()) != grid.nrow()) {
                npoint += 1;
            }

            const ne = @as(usize, @intFromFloat(@round((self.e_max - self.e_min) / self.e_step))) + 1;

            var sigma = try Matrix(T).initZero(ne, wfn_init.nrow(), gpa);
            errdefer sigma.deinit(gpa);

            var temp_phi = try gpa.alloc(Complex(T), grid.nrow());
            defer gpa.free(temp_phi);

            const shape = try gpa.alloc(i32, grid.ncol());
            defer gpa.free(shape);

            for (0..grid.ncol()) |i| {
                shape[i] = @as(i32, @intCast(npoint));
            }

            const ffft_plan = try FftPlan(Complex(T)).init(temp_phi, shape, -1, fftw.FFTW_ESTIMATE);
            defer ffft_plan.deinit();

            const ifft_plan = try FftPlan(Complex(T)).init(temp_phi, shape, 1, fftw.FFTW_ESTIMATE);
            defer ifft_plan.deinit();

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
                    const E = self.e_min + @as(T, @floatFromInt(ei)) * self.e_step;

                    if (E <= self.e_thrs) {
                        continue;
                    }

                    const k_inc = std.math.sqrt(2 * m * (E - self.e_thrs));

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

                        ffft_plan.execute(temp_phi);

                        for (0..grid.nrow()) |m_idx| {
                            temp_phi[m_idx] = temp_phi[m_idx].mul(Complex(T).init(0, grid.getK(m_idx, d)));
                        }

                        ifft_plan.execute(temp_phi);

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
