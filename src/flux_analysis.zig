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
        mass: []const T,

        e_min: T,
        e_max: T,
        gamma: T,
        initk: T,

        e_step: T,
        Vreact: T,

        istate: usize,
        grid: Grid(T),

        flux_bounds: []const [2]T,
        prop_direction: []const T,

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

            const mass = try gpa.alloc(T, opt.mass.len);
            errdefer gpa.free(mass);

            for (opt.mass, 0..) |m, i| {
                mass[i] = m;
            }

            const flux_bounds = try gpa.alloc([2]T, flux_opt.flux_bounds.len);
            errdefer gpa.free(flux_bounds);

            for (flux_opt.flux_bounds, 0..) |b, i| {
                flux_bounds[i] = .{ b[0], b[1] };
            }

            const prop_direction = try gpa.alloc(T, opt.initial_conditions.momentum.len);
            errdefer gpa.free(prop_direction);

            var p_norm: T = 0;

            for (opt.initial_conditions.momentum) |p| {
                p_norm += p * p;
            }

            p_norm = std.math.sqrt(p_norm);

            for (opt.initial_conditions.momentum, 0..) |p, i| {
                prop_direction[i] = if (p_norm > 0) p / p_norm else 0;
            }

            var sum_inv_gamma: T = 0;

            for (opt.initial_conditions.gamma, 0..) |g, i| {
                sum_inv_gamma += (prop_direction[i] * prop_direction[i]) / g;
            }

            const gamma = if (sum_inv_gamma > 0) 1 / sum_inv_gamma else 0;

            return .{
                .flux_bounds = flux_bounds,
                .e_min = flux_opt.e_min,
                .e_max = flux_opt.e_max,
                .e_step = flux_opt.e_step,
                .mass = mass,
                .dt = opt.time_step,
                .istate = opt.initial_conditions.state,
                .Vreact = Vreact,
                .grid = grid,
                .prop_direction = prop_direction,
                .initk = p_norm,
                .gamma = gamma,
            };
        }

        /// Deallocates the dynamically allocated slices in the context.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            gpa.free(self.mass);
            gpa.free(self.flux_bounds);
            gpa.free(self.prop_direction);
        }

        /// Computes transition probabilities by integrating flux of energy-resolved wavefunctions at a dividing surface.
        pub fn analyze(self: @This(), wfn_init: Matrix(Complex(T)), flux_acc: Matrix(Complex(T)), gpa: Allocator) !Matrix(T) {
            var npoint: usize, var inv_m_eff: T = .{ 1, 0 };

            while (try std.math.powi(usize, npoint, self.grid.r.ncol()) != self.grid.r.nrow()) {
                npoint += 1;
            }

            const ne = @as(usize, @intFromFloat(@round((self.e_max - self.e_min) / self.e_step))) + 1;

            for (self.prop_direction, 0..) |e_p, i| if (self.mass[i] > 0) {
                inv_m_eff += (e_p * e_p) / self.mass[i];
            };

            const m_eff = if (inv_m_eff > 0) 1 / inv_m_eff else 0;

            var prob_matrix = try Matrix(T).initZero(ne, wfn_init.nrow(), gpa);
            errdefer prob_matrix.deinit(gpa);

            var temp_phi = try gpa.alloc(Complex(T), self.grid.r.nrow());
            defer gpa.free(temp_phi);

            const shape = try gpa.alloc(i32, self.grid.r.ncol());
            defer gpa.free(shape);

            for (0..self.grid.r.ncol()) |i| {
                shape[i] = @as(i32, @intCast(npoint));
            }

            const ffft_plan = try FftPlan(Complex(T)).init(temp_phi, shape, -1, fftw.FFTW_ESTIMATE);
            defer ffft_plan.deinit();

            const ifft_plan = try FftPlan(Complex(T)).init(temp_phi, shape, 1, fftw.FFTW_ESTIMATE);
            defer ifft_plan.deinit();

            for (0..self.grid.r.ncol()) |d| {
                const s_d, const r = .{ std.math.pow(usize, npoint, self.grid.r.ncol() - 1 - d), self.grid.r };

                const n_min_f = (self.flux_bounds[d][0] - r.at(0, d)) / (r.at(s_d, d) - r.at(0, d));
                const n_max_f = (self.flux_bounds[d][1] - r.at(0, d)) / (r.at(s_d, d) - r.at(0, d));

                const n_min: usize = @intFromFloat(@round(n_min_f));
                const n_max: usize = @intFromFloat(@round(n_max_f));

                if (n_min == 0 or n_max >= npoint - 1) {
                    std.log.err("FLUX BOUNDS MUST LIE WITHIN THE GRID INTERIOR", .{});

                    return error.InvalidInput;
                }

                const dx_d, const dr = .{ r.at(s_d, d) - r.at(0, d), self.grid.dr };

                for (0..ne) |ei| {
                    const E = self.e_min + @as(T, @floatFromInt(ei)) * self.e_step;

                    if (E <= self.Vreact) {
                        continue;
                    }

                    const k_inc = std.math.sqrt(2 * m_eff * (E - self.Vreact));

                    const exp_arg = -std.math.pow(T, k_inc - self.initk, @as(T, 2)) / (self.gamma);
                    const ak = std.math.sqrt(4 * std.math.pi / self.gamma) * std.math.exp(exp_arg);

                    if (ak / (dr * dr) < 1e-6) {
                        continue;
                    }

                    for (0..wfn_init.nrow()) |f| {
                        const row = ei * wfn_init.nrow() + f;

                        var sum: T = 0;

                        for (0..r.nrow()) |j| {
                            temp_phi[j] = flux_acc.at(row, j);
                        }

                        ffft_plan.execute(temp_phi);

                        for (0..r.nrow()) |m| {
                            temp_phi[m] = temp_phi[m].mul(Complex(T).init(0, self.grid.k.at(m, d)));
                        }

                        ifft_plan.execute(temp_phi);

                        const factor = dr * k_inc * self.dt * self.dt / (m_eff * self.mass[d] * dx_d * ak);

                        for (0..r.nrow()) |i| {
                            var in_bounds = true;

                            for (0..r.ncol()) |k| {
                                if (k == d) continue;

                                const s_k = std.math.pow(usize, npoint, r.ncol() - 1 - k);

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

                            if ((i / s_d) % npoint == n_max) sum += factor * val;
                            if ((i / s_d) % npoint == n_min) sum -= factor * val;
                        }

                        prob_matrix.ptr(ei, f).* += sum;
                    }
                }
            }

            return prob_matrix;
        }
    };
}
