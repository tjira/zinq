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

            const mass = try gpa.alloc(T, opt.mass.len);
            errdefer gpa.free(mass);

            for (opt.mass, 0..) |m, i| {
                mass[i] = m;
            }

            const gamma = try gpa.alloc(T, opt.initial_conditions.gamma.len);
            errdefer gpa.free(gamma);

            for (opt.initial_conditions.gamma, 0..) |g, i| {
                gamma[i] = g;
            }

            const flux_bounds = try gpa.alloc([2]T, flux_opt.flux_bounds.len);
            errdefer gpa.free(flux_bounds);

            for (flux_opt.flux_bounds, 0..) |b, i| {
                flux_bounds[i] = .{ b[0], b[1] };
            }

            var e_perp: T, var r_perp: T = .{ 0, 1 };

            for (1..opt.mass.len) |i| {
                e_perp += opt.initial_conditions.gamma[i] / (2 * mass[i]);
            }

            for (1..opt.initial_conditions.gamma.len) |i| {
                if (opt.cylindrical and i == opt.initial_conditions.gamma.len - 1) {
                    continue;
                }

                r_perp *= std.math.sqrt(opt.initial_conditions.gamma[i] / std.math.pi);
            }

            const weight = if (opt.cylindrical) (std.math.pi / gamma[gamma.len - 1]) / r_perp else 1 / r_perp;

            return .{
                .flux_bounds = flux_bounds,
                .e_min = flux_opt.e_min,
                .e_max = flux_opt.e_max,
                .e_step = flux_opt.e_step,
                .mass = mass,
                .dt = opt.time_step,
                .e_thrs = Vreact + e_perp,
                .initk = @abs(opt.initial_conditions.momentum[0]),
                .gmma = gamma,
                .weight = weight,
            };
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

            while (try std.math.powi(usize, npoint, grid.r.ncol()) != grid.r.nrow()) {
                npoint += 1;
            }

            const ne = @as(usize, @intFromFloat(@round((self.e_max - self.e_min) / self.e_step))) + 1;

            var sigma = try Matrix(T).initZero(ne, wfn_init.nrow(), gpa);
            errdefer sigma.deinit(gpa);

            var temp_phi = try gpa.alloc(Complex(T), grid.r.nrow());
            defer gpa.free(temp_phi);

            const shape = try gpa.alloc(i32, grid.r.ncol());
            defer gpa.free(shape);

            for (0..grid.r.ncol()) |i| {
                shape[i] = @as(i32, @intCast(npoint));
            }

            const ffft_plan = try FftPlan(Complex(T)).init(temp_phi, shape, -1, fftw.FFTW_ESTIMATE);
            defer ffft_plan.deinit();

            const ifft_plan = try FftPlan(Complex(T)).init(temp_phi, shape, 1, fftw.FFTW_ESTIMATE);
            defer ifft_plan.deinit();

            for (0..grid.r.ncol()) |d| {
                const s_d, const r = .{ std.math.pow(usize, npoint, grid.r.ncol() - 1 - d), grid.r };

                const n_min_f = (self.flux_bounds[d][0] - r.at(0, d)) / (r.at(s_d, d) - r.at(0, d));
                const n_max_f = (self.flux_bounds[d][1] - r.at(0, d)) / (r.at(s_d, d) - r.at(0, d));

                if (n_min_f < -0.5 or n_max_f >= @as(T, @floatFromInt(npoint)) - 0.5) {
                    std.log.err("FLUX BOUNDS MUST LIE WITHIN THE GRID BOUNDS", .{});

                    return error.InvalidInput;
                }

                const n_min: usize = @intFromFloat(@round(n_min_f));
                const n_max: usize = @intFromFloat(@round(n_max_f));

                const dx_d, const m = .{ r.at(s_d, d) - r.at(0, d), self.mass[0] };

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

                        for (0..grid.r.nrow()) |j| {
                            temp_phi[j] = flux_acc.at(row, j);
                        }

                        ffft_plan.execute(temp_phi);

                        for (0..grid.r.nrow()) |m_idx| {
                            temp_phi[m_idx] = temp_phi[m_idx].mul(Complex(T).init(0, grid.k.at(m_idx, d)));
                        }

                        ifft_plan.execute(temp_phi);

                        const factor = grid.dr * k_inc * self.dt * self.dt / (m * self.mass[d] * dx_d * ak);

                        for (0..grid.r.nrow()) |i| {
                            var in_bounds = true;

                            for (0..grid.r.ncol()) |k| {
                                if (k == d) continue;

                                const s_k = std.math.pow(usize, npoint, grid.r.ncol() - 1 - k);

                                const n_min_k_f = (self.flux_bounds[k][0] - r.at(0, k)) / (r.at(s_k, k) - r.at(0, k));
                                const n_max_k_f = (self.flux_bounds[k][1] - r.at(0, k)) / (r.at(s_k, k) - r.at(0, k));

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
