const std = @import("std");
const libint = @import("libint");

const libxc = @cImport(@cInclude("xc.h"));

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Value = @import("value.zig").Value;
const Vector = @import("tensor.zig").Vector;

const getLebedevGrid = @import("lebedev_quadrature_nodes.zig").getLebedevGrid;
const primType = @import("value.zig").primType;

const A2BOHR = @import("constant.zig").A2BOHR;

// DFT STRUCTURE =======================================================================================================

pub fn DftPotential(comptime T: type) type {
    return struct {
        exch: ?libxc.xc_func_type,
        corr: ?libxc.xc_func_type,
        exco: ?libxc.xc_func_type,

        pts: Matrix(T),
        wgh: Vector(T),

        Vxc: Matrix(T),
        Exc: T = 0.000,

        polarized: bool,

        pub fn init(sys: MolecularSystem(T), funcs: anytype, n_rad: usize, n_leb: usize, polarized: bool, gpa: Allocator) !@This() {
            var grid = try (Becke(T){ .n_rad = n_rad, .n_leb = n_leb }).get(sys, gpa);

            errdefer grid[0].deinit(gpa);
            errdefer grid[1].deinit(gpa);

            const nbf = if (polarized) 2 * sys.nbf else sys.nbf;

            var exch_func: ?libxc.xc_func_type = null;
            var corr_func: ?libxc.xc_func_type = null;
            var exco_func: ?libxc.xc_func_type = null;

            const mode = if (polarized) libxc.XC_POLARIZED else libxc.XC_UNPOLARIZED;

            if (funcs[0]) |name| {
                const exch_nt = try gpa.dupeSentinel(u8, name, 0);
                defer gpa.free(exch_nt);

                const id = libxc.xc_functional_get_number(exch_nt.ptr);

                if (id < 0) return error.ExchFuncNotFound;

                var func: libxc.xc_func_type = undefined;

                if (libxc.xc_func_init(&func, id, mode) != 0) {
                    return error.LibxcInitFailed;
                }

                if (func.info.*.kind != libxc.XC_EXCHANGE) {
                    libxc.xc_func_end(&func);

                    return error.InvalidExchangeFunctional;
                }

                exch_func = func;
            }

            errdefer {
                if (exch_func) |*f| libxc.xc_func_end(f);
            }

            if (funcs[1]) |name| {
                const corr_nt = try gpa.dupeSentinel(u8, name, 0);
                defer gpa.free(corr_nt);

                const id = libxc.xc_functional_get_number(corr_nt.ptr);

                if (id < 0) return error.CorrFuncNotFound;

                var func: libxc.xc_func_type = undefined;

                if (libxc.xc_func_init(&func, id, mode) != 0) {
                    return error.LibxcInitFailed;
                }

                if (func.info.*.kind != libxc.XC_CORRELATION) {
                    libxc.xc_func_end(&func);

                    return error.InvalidCorrelationFunctional;
                }

                corr_func = func;
            }

            errdefer {
                if (corr_func) |*f| libxc.xc_func_end(f);
            }

            if (funcs[2]) |name| {
                const exco_nt = try gpa.dupeSentinel(u8, name, 0);
                defer gpa.free(exco_nt);

                const id = libxc.xc_functional_get_number(exco_nt.ptr);

                if (id < 0) return error.XcFuncNotFound;

                var func: libxc.xc_func_type = undefined;

                if (libxc.xc_func_init(&func, id, mode) != 0) {
                    return error.LibxcInitFailed;
                }

                if (func.info.*.kind != libxc.XC_EXCHANGE_CORRELATION) {
                    libxc.xc_func_end(&func);

                    return error.InvalidExchangeCorrelationFunctional;
                }

                exco_func = func;
            }

            errdefer {
                if (exco_func) |*f| libxc.xc_func_end(f);
            }

            if (exco_func != null and (exch_func != null or corr_func != null)) {
                return error.InvalidFunctionalCombination;
            }

            if (exch_func == null and corr_func == null and exco_func == null) {
                return error.NoFunctionalSpecified;
            }

            var Vxc = try Matrix(T).init(nbf, nbf, gpa);
            errdefer Vxc.deinit(gpa);

            return .{ .exch = exch_func, .corr = corr_func, .exco = exco_func, .pts = grid[0], .wgh = grid[1], .Vxc = Vxc, .polarized = polarized };
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            if (self.exch) |*func| libxc.xc_func_end(func);
            if (self.corr) |*func| libxc.xc_func_end(func);
            if (self.exco) |*func| libxc.xc_func_end(func);

            self.pts.deinit(gpa);
            self.wgh.deinit(gpa);
            self.Vxc.deinit(gpa);
        }

        pub fn evaluate(self: *@This(), sys: MolecularSystem(T), P: Matrix(T), gpa: Allocator) !void {
            self.Vxc.zero();

            const is_exch_sgga = self.exch != null and self.exch.?.info.*.family == libxc.XC_FAMILY_GGA;
            const is_corr_sgga = self.corr != null and self.corr.?.info.*.family == libxc.XC_FAMILY_GGA;
            const is_exco_sgga = self.exco != null and self.exco.?.info.*.family == libxc.XC_FAMILY_GGA;

            const is_exch_mgga = self.exch != null and self.exch.?.info.*.family == libxc.XC_FAMILY_MGGA;
            const is_corr_mgga = self.corr != null and self.corr.?.info.*.family == libxc.XC_FAMILY_MGGA;
            const is_exco_mgga = self.exco != null and self.exco.?.info.*.family == libxc.XC_FAMILY_MGGA;

            const has_sgga = is_exch_sgga or is_corr_sgga or is_exco_sgga;
            const has_mgga = is_exch_mgga or is_corr_mgga or is_exco_mgga;

            const size_factor: usize = if (self.polarized) 2 else 1;

            var exc = try Vector(T).initZero(self.pts.shape[0], gpa);
            defer exc.deinit(gpa);

            var rho_val = try Matrix(T).initZero(self.pts.shape[0], size_factor, gpa);
            defer rho_val.deinit(gpa);

            var rho_pot = try Matrix(T).initZero(self.pts.shape[0], size_factor, gpa);
            defer rho_pot.deinit(gpa);

            var sig_val: ?Matrix(T) = null;
            var sig_pot: ?Matrix(T) = null;

            var tau_val: ?Matrix(T) = null;
            var tau_pot: ?Matrix(T) = null;

            var del_rho: ?Matrix(T) = null;

            var dphi_dx: ?Matrix(T) = null;
            var dphi_dy: ?Matrix(T) = null;
            var dphi_dz: ?Matrix(T) = null;

            var X_a: ?Vector(T) = null;
            var X_b: ?Vector(T) = null;

            const needs_deriv = has_sgga or has_mgga;

            if (needs_deriv) {
                const sigma_cols: usize = if (self.polarized) 3 else 1;

                sig_val = try Matrix(T).initZero(self.pts.shape[0], sigma_cols, gpa);
                sig_pot = try Matrix(T).initZero(self.pts.shape[0], sigma_cols, gpa);

                const grad_cols: usize = if (self.polarized) 6 else 3;

                del_rho = try Matrix(T).initZero(self.pts.shape[0], grad_cols, gpa);

                dphi_dx = try Matrix(T).init(self.pts.shape[0], sys.nbf, gpa);
                dphi_dy = try Matrix(T).init(self.pts.shape[0], sys.nbf, gpa);
                dphi_dz = try Matrix(T).init(self.pts.shape[0], sys.nbf, gpa);

                X_a = try Vector(T).init(sys.nbf, gpa);

                if (self.polarized) {
                    X_b = try Vector(T).init(sys.nbf, gpa);
                }
            }

            if (has_mgga) {
                tau_val = try Matrix(T).initZero(self.pts.shape[0], size_factor, gpa);
                tau_pot = try Matrix(T).initZero(self.pts.shape[0], size_factor, gpa);
            }

            defer {
                if (sig_val) |*s| s.deinit(gpa);
                if (sig_pot) |*s| s.deinit(gpa);
                if (del_rho) |*g| g.deinit(gpa);
                if (tau_val) |*t| t.deinit(gpa);
                if (tau_pot) |*t| t.deinit(gpa);
                if (dphi_dx) |*d| d.deinit(gpa);
                if (dphi_dy) |*d| d.deinit(gpa);
                if (dphi_dz) |*d| d.deinit(gpa);

                if (X_a) |*x| x.deinit(gpa);
                if (X_b) |*x| x.deinit(gpa);
            }

            var phi = try Matrix(T).init(self.pts.shape[0], sys.nbf, gpa);
            defer phi.deinit(gpa);

            for (0..self.pts.shape[0]) |i| {
                const x = self.pts.at(i, 0);
                const y = self.pts.at(i, 1);
                const z = self.pts.at(i, 2);

                if (needs_deriv) {
                    const dx_row = dphi_dx.?.rowSlice(i);
                    const dy_row = dphi_dy.?.rowSlice(i);
                    const dz_row = dphi_dz.?.rowSlice(i);

                    libint.evaluate_basis_derivative(phi.rowSlice(i).ptr, dx_row.ptr, dy_row.ptr, dz_row.ptr, x, y, z, sys.ptr);
                }

                if (!needs_deriv) {
                    libint.evaluate_basis(phi.rowSlice(i).ptr, x, y, z, sys.ptr);
                }

                if (self.polarized) {
                    var rho_a: T = 0;
                    var rho_b: T = 0;

                    var del_rho_a_x: T = 0;
                    var del_rho_a_y: T = 0;
                    var del_rho_a_z: T = 0;
                    var del_rho_b_x: T = 0;
                    var del_rho_b_y: T = 0;
                    var del_rho_b_z: T = 0;

                    var tau_a: T = 0;
                    var tau_b: T = 0;

                    for (0..sys.nbf) |mu| for (0..sys.nbf) |nu| {
                        const b = sys.nbf;

                        const p_val_a = P.at(mu + 0, nu + 0);
                        const p_val_b = P.at(mu + b, nu + b);

                        const phi_mu = phi.at(i, mu);
                        const phi_nu = phi.at(i, nu);

                        rho_a += p_val_a * phi_mu * phi_nu;
                        rho_b += p_val_b * phi_mu * phi_nu;

                        if (needs_deriv) {
                            const dphi_dx_mu = dphi_dx.?.at(i, mu);
                            const dphi_dy_mu = dphi_dy.?.at(i, mu);
                            const dphi_dz_mu = dphi_dz.?.at(i, mu);

                            const term_a = p_val_a * phi_nu;
                            const term_b = p_val_b * phi_nu;

                            del_rho_a_x += term_a * dphi_dx_mu;
                            del_rho_a_y += term_a * dphi_dy_mu;
                            del_rho_a_z += term_a * dphi_dz_mu;

                            del_rho_b_x += term_b * dphi_dx_mu;
                            del_rho_b_y += term_b * dphi_dy_mu;
                            del_rho_b_z += term_b * dphi_dz_mu;

                            if (has_mgga) {
                                const dphi_dx_nu = dphi_dx.?.at(i, nu);
                                const dphi_dy_nu = dphi_dy.?.at(i, nu);
                                const dphi_dz_nu = dphi_dz.?.at(i, nu);

                                const dphi_dot = dphi_dx_mu * dphi_dx_nu + dphi_dy_mu * dphi_dy_nu + dphi_dz_mu * dphi_dz_nu;

                                tau_a += p_val_a * dphi_dot;
                                tau_b += p_val_b * dphi_dot;
                            }
                        }
                    };

                    rho_val.ptr(i, 0).* = rho_a;
                    rho_val.ptr(i, 1).* = rho_b;

                    if (needs_deriv) {
                        del_rho_a_x *= 2;
                        del_rho_a_y *= 2;
                        del_rho_a_z *= 2;
                        del_rho_b_x *= 2;
                        del_rho_b_y *= 2;
                        del_rho_b_z *= 2;

                        del_rho.?.ptr(i, 0).* = del_rho_a_x;
                        del_rho.?.ptr(i, 1).* = del_rho_a_y;
                        del_rho.?.ptr(i, 2).* = del_rho_a_z;
                        del_rho.?.ptr(i, 3).* = del_rho_b_x;
                        del_rho.?.ptr(i, 4).* = del_rho_b_y;
                        del_rho.?.ptr(i, 5).* = del_rho_b_z;

                        sig_val.?.ptr(i, 0).* = del_rho_a_x * del_rho_a_x + del_rho_a_y * del_rho_a_y + del_rho_a_z * del_rho_a_z;
                        sig_val.?.ptr(i, 1).* = del_rho_a_x * del_rho_b_x + del_rho_a_y * del_rho_b_y + del_rho_a_z * del_rho_b_z;
                        sig_val.?.ptr(i, 2).* = del_rho_b_x * del_rho_b_x + del_rho_b_y * del_rho_b_y + del_rho_b_z * del_rho_b_z;
                    }

                    if (has_mgga) {
                        tau_val.?.ptr(i, 0).* = 0.5 * tau_a;
                        tau_val.?.ptr(i, 1).* = 0.5 * tau_b;
                    }
                }

                if (!self.polarized) {
                    var rho_tot: T = 0;
                    var tau_tot: T = 0;

                    var del_rho_x: T = 0;
                    var del_rho_y: T = 0;
                    var del_rho_z: T = 0;

                    for (0..sys.nbf) |mu| for (0..sys.nbf) |nu| {
                        const phi_mu = phi.at(i, mu);
                        const phi_nu = phi.at(i, nu);

                        rho_tot += P.at(mu, nu) * phi_mu * phi_nu;

                        if (needs_deriv) {
                            const dphi_dx_mu = dphi_dx.?.at(i, mu);
                            const dphi_dy_mu = dphi_dy.?.at(i, mu);
                            const dphi_dz_mu = dphi_dz.?.at(i, mu);

                            const term = P.at(mu, nu) * phi_nu;

                            del_rho_x += term * dphi_dx_mu;
                            del_rho_y += term * dphi_dy_mu;
                            del_rho_z += term * dphi_dz_mu;

                            if (has_mgga) {
                                const dphi_dx_nu = dphi_dx.?.at(i, nu);
                                const dphi_dy_nu = dphi_dy.?.at(i, nu);
                                const dphi_dz_nu = dphi_dz.?.at(i, nu);

                                const dphi_dot = dphi_dx_mu * dphi_dx_nu + dphi_dy_mu * dphi_dy_nu + dphi_dz_mu * dphi_dz_nu;

                                tau_tot += P.at(mu, nu) * dphi_dot;
                            }
                        }
                    };

                    rho_val.ptr(i, 0).* = rho_tot;

                    if (needs_deriv) {
                        del_rho_x *= 2;
                        del_rho_y *= 2;
                        del_rho_z *= 2;

                        del_rho.?.ptr(i, 0).* = del_rho_x;
                        del_rho.?.ptr(i, 1).* = del_rho_y;
                        del_rho.?.ptr(i, 2).* = del_rho_z;

                        sig_val.?.ptr(i, 0).* = del_rho_x * del_rho_x + del_rho_y * del_rho_y + del_rho_z * del_rho_z;
                    }

                    if (has_mgga) {
                        tau_val.?.ptr(i, 0).* = 0.5 * tau_tot;
                    }
                }
            }

            const inp = XCInput(T){
                .rho_val = rho_val.data,

                .sig_val = if (sig_val) |s| s.data else null,
                .tau_val = if (tau_val) |t| t.data else null,
            };

            const out = XCOutput(T){
                .exc = exc.data,

                .sig_pot = if (sig_pot) |s| s.data else null,
                .tau_pot = if (tau_pot) |t| t.data else null,

                .rho_pot = rho_pot.data,
            };

            if (self.exch) |func| evaluateXCFunctional(T, out, func, self.polarized, inp);
            if (self.corr) |func| evaluateXCFunctional(T, out, func, self.polarized, inp);
            if (self.exco) |func| evaluateXCFunctional(T, out, func, self.polarized, inp);

            var energy_exc: T = 0;

            for (0..self.pts.shape[0]) |i| {
                const w = self.wgh.at(i);

                if (self.polarized) {
                    const rho_tot = rho_val.at(i, 0) + rho_val.at(i, 1);

                    energy_exc += rho_tot * exc.at(i) * w;

                    const v_a = rho_pot.at(i, 0);
                    const v_b = rho_pot.at(i, 1);

                    if (needs_deriv) {
                        const g_a_x = del_rho.?.at(i, 0);
                        const g_a_y = del_rho.?.at(i, 1);
                        const g_a_z = del_rho.?.at(i, 2);

                        const g_b_x = del_rho.?.at(i, 3);
                        const g_b_y = del_rho.?.at(i, 4);
                        const g_b_z = del_rho.?.at(i, 5);

                        const v_sig_aa = sig_pot.?.at(i, 0);
                        const v_sig_ab = sig_pot.?.at(i, 1);
                        const v_sig_bb = sig_pot.?.at(i, 2);

                        const w_a_x = 2 * v_sig_aa * g_a_x + v_sig_ab * g_b_x;
                        const w_a_y = 2 * v_sig_aa * g_a_y + v_sig_ab * g_b_y;
                        const w_a_z = 2 * v_sig_aa * g_a_z + v_sig_ab * g_b_z;

                        const w_b_x = 2 * v_sig_bb * g_b_x + v_sig_ab * g_a_x;
                        const w_b_y = 2 * v_sig_bb * g_b_y + v_sig_ab * g_a_y;
                        const w_b_z = 2 * v_sig_bb * g_b_z + v_sig_ab * g_a_z;

                        for (0..sys.nbf) |mu| {
                            const d_dx = dphi_dx.?.at(i, mu);
                            const d_dy = dphi_dy.?.at(i, mu);
                            const d_dz = dphi_dz.?.at(i, mu);

                            X_a.?.data[mu] = w_a_x * d_dx + w_a_y * d_dy + w_a_z * d_dz;
                            X_b.?.data[mu] = w_b_x * d_dx + w_b_y * d_dy + w_b_z * d_dz;
                        }

                        const v_tau_a = if (has_mgga) tau_pot.?.at(i, 0) else 0.0;
                        const v_tau_b = if (has_mgga) tau_pot.?.at(i, 1) else 0.0;

                        for (0..sys.nbf) |mu| for (0..sys.nbf) |nu| {
                            const phi_mu = phi.at(i, mu);
                            const phi_nu = phi.at(i, nu);

                            const b = sys.nbf;

                            var val_a = v_a * phi_mu * phi_nu + X_a.?.data[mu] * phi_nu + phi_mu * X_a.?.data[nu];
                            var val_b = v_b * phi_mu * phi_nu + X_b.?.data[mu] * phi_nu + phi_mu * X_b.?.data[nu];

                            if (has_mgga) {
                                const dx_mu_nu = dphi_dx.?.at(i, mu) * dphi_dx.?.at(i, nu);
                                const dy_mu_nu = dphi_dy.?.at(i, mu) * dphi_dy.?.at(i, nu);
                                const dz_mu_nu = dphi_dz.?.at(i, mu) * dphi_dz.?.at(i, nu);

                                const dphi_dot = dx_mu_nu + dy_mu_nu + dz_mu_nu;

                                val_a += 0.5 * v_tau_a * dphi_dot;
                                val_b += 0.5 * v_tau_b * dphi_dot;
                            }

                            self.Vxc.ptr(mu + 0, nu + 0).* += val_a * w;
                            self.Vxc.ptr(mu + b, nu + b).* += val_b * w;
                        };
                    }

                    if (!needs_deriv) {
                        for (0..sys.nbf) |mu| for (0..sys.nbf) |nu| {
                            const phi_mu = phi.at(i, mu);
                            const phi_nu = phi.at(i, nu);

                            const b = sys.nbf;

                            self.Vxc.ptr(mu + 0, nu + 0).* += v_a * phi_mu * phi_nu * w;
                            self.Vxc.ptr(mu + b, nu + b).* += v_b * phi_mu * phi_nu * w;
                        };
                    }
                }

                if (!self.polarized) {
                    energy_exc += rho_val.at(i, 0) * exc.at(i) * w;

                    if (needs_deriv) {
                        const v_sigma = sig_pot.?.at(i, 0);

                        const w_x = 2 * v_sigma * del_rho.?.at(i, 0);
                        const w_y = 2 * v_sigma * del_rho.?.at(i, 1);
                        const w_z = 2 * v_sigma * del_rho.?.at(i, 2);

                        for (0..sys.nbf) |mu| {
                            X_a.?.data[mu] = w_x * dphi_dx.?.at(i, mu) + w_y * dphi_dy.?.at(i, mu) + w_z * dphi_dz.?.at(i, mu);
                        }

                        const v_tau = if (has_mgga) tau_pot.?.at(i, 0) else 0.0;

                        for (0..sys.nbf) |mu| for (0..sys.nbf) |nu| {
                            const phi_mu = phi.at(i, mu);
                            const phi_nu = phi.at(i, nu);

                            var val = rho_pot.at(i, 0) * phi_mu * phi_nu + X_a.?.data[mu] * phi_nu + phi_mu * X_a.?.data[nu];

                            if (has_mgga) {
                                const dx_mu_nu = dphi_dx.?.at(i, mu) * dphi_dx.?.at(i, nu);
                                const dy_mu_nu = dphi_dy.?.at(i, mu) * dphi_dy.?.at(i, nu);
                                const dz_mu_nu = dphi_dz.?.at(i, mu) * dphi_dz.?.at(i, nu);

                                val += 0.5 * v_tau * (dx_mu_nu + dy_mu_nu + dz_mu_nu);
                            }

                            self.Vxc.ptr(mu, nu).* += val * w;
                        };
                    }

                    if (!needs_deriv) for (0..sys.nbf) |mu| for (0..sys.nbf) |nu| {
                        self.Vxc.ptr(mu, nu).* += rho_pot.at(i, 0) * phi.at(i, mu) * phi.at(i, nu) * w;
                    };
                }
            }

            self.Exc = @floatCast(energy_exc);
        }
    };
}

// GRID STRUCTURES =====================================================================================================

pub fn Becke(comptime T: type) type {
    return struct {
        n_rad: usize,
        n_leb: usize,

        pub fn get(self: @This(), sys: MolecularSystem(T), gpa: Allocator) !struct { Matrix(T), Vector(T) } {
            var centers = std.ArrayList([3]T).empty;
            defer centers.deinit(gpa);

            var atomic_numbers = std.ArrayList(usize).empty;
            defer atomic_numbers.deinit(gpa);

            const lebedev_points = try getLebedevGrid(self.n_leb);

            for (0..sys.atoms.len) |i| {
                try atomic_numbers.append(gpa, @intCast(sys.atoms[i]));

                const c = [3]T{
                    @floatCast(sys.coors[3 * i + 0]),
                    @floatCast(sys.coors[3 * i + 1]),
                    @floatCast(sys.coors[3 * i + 2]),
                };

                try centers.append(gpa, c);
            }

            const total_pts = centers.items.len * self.n_rad * self.n_leb;

            var pts = try Matrix(T).init(total_pts, 3, gpa);
            errdefer pts.deinit(gpa);

            var wgh = try Vector(T).init(total_pts, gpa);
            errdefer wgh.deinit(gpa);

            var W_vals = try gpa.alloc(T, centers.items.len);
            defer gpa.free(W_vals);

            var R_AB = try gpa.alloc(T, centers.items.len * centers.items.len);
            defer gpa.free(R_AB);

            var a_AB = try gpa.alloc(T, centers.items.len * centers.items.len);
            defer gpa.free(a_AB);

            for (0..centers.items.len) |A| for (0..centers.items.len) |B| {
                const dx = centers.items[A][0] - centers.items[B][0];
                const dy = centers.items[A][1] - centers.items[B][1];
                const dz = centers.items[A][2] - centers.items[B][2];

                R_AB[A * centers.items.len + B] = std.math.sqrt(dx * dx + dy * dy + dz * dz);

                if (A == B) a_AB[A * centers.items.len + B] = 0;

                if (A != B) {
                    const radius_A = getBraggSlaterRadius(T, atomic_numbers.items[A]);
                    const radius_B = getBraggSlaterRadius(T, atomic_numbers.items[B]);

                    const u = (radius_A - radius_B) / (radius_A + radius_B);

                    a_AB[A * centers.items.len + B] = std.math.clamp(u / (u * u - 1), -0.5, 0.5);
                }
            };

            var idx: usize = 0;

            for (0..centers.items.len) |A| for (0..self.n_rad) |ir| {
                const r, const w_r = self.getRadialNode(ir, getBraggSlaterRadius(T, atomic_numbers.items[A]));

                for (lebedev_points.nodes, lebedev_points.nweights) |lebedev_node, lebedev_weight| {
                    pts.ptr(idx, 0).* = centers.items[A][0] + r * @as(T, @floatCast(lebedev_node[0]));
                    pts.ptr(idx, 1).* = centers.items[A][1] + r * @as(T, @floatCast(lebedev_node[1]));
                    pts.ptr(idx, 2).* = centers.items[A][2] + r * @as(T, @floatCast(lebedev_node[2]));

                    var P_sum: T = 0;

                    for (0..centers.items.len) |I| {
                        var W_I: T = 1;

                        const px = pts.at(idx, 0);
                        const py = pts.at(idx, 1);
                        const pz = pts.at(idx, 2);

                        const dx_I = px - centers.items[I][0];
                        const dy_I = py - centers.items[I][1];
                        const dz_I = pz - centers.items[I][2];

                        const r_I = std.math.sqrt(dx_I * dx_I + dy_I * dy_I + dz_I * dz_I);

                        for (0..centers.items.len) |J| {
                            if (I == J) continue;

                            const dx_J = px - centers.items[J][0];
                            const dy_J = py - centers.items[J][1];
                            const dz_J = pz - centers.items[J][2];

                            const r_J = std.math.sqrt(dx_J * dx_J + dy_J * dy_J + dz_J * dz_J);

                            const mu = (r_I - r_J) / R_AB[I * centers.items.len + J];

                            W_I *= beckeStep(mu + a_AB[I * centers.items.len + J] * (1 - mu * mu));
                        }

                        W_vals[I], P_sum = .{ W_I, P_sum + W_I };
                    }

                    wgh.ptr(idx).*, idx = .{ 4 * std.math.pi * w_r * lebedev_weight * W_vals[A] / P_sum, idx + 1 };
                }
            };

            return .{ pts, wgh };
        }

        fn beckeStep(mu: T) T {
            var p = 1.5 * mu - 0.5 * mu * mu * mu;

            p = 1.5 * p - 0.5 * p * p * p;
            p = 1.5 * p - 0.5 * p * p * p;

            return 0.5 * (1.0 - p);
        }

        fn getRadialNode(self: @This(), ir: usize, br: T) struct { T, T } {
            const nr = self.n_rad;

            const f_ir = @as(T, @floatFromInt(ir));
            const f_nr = @as(T, @floatFromInt(nr));

            const x_cheb = std.math.cos(std.math.pi / f_nr * (f_ir + 0.5));
            const w_cheb = std.math.pi / f_nr * @sqrt(1 - x_cheb * x_cheb);

            const r = br * (1 + x_cheb) / (1 - x_cheb);

            return .{ r, 2 * w_cheb * br * r * r / ((1 - x_cheb) * (1 - x_cheb)) };
        }
    };
}

// HELPER FUNCTIONS ====================================================================================================

fn getBraggSlaterRadius(comptime T: type, Z: usize) T {
    const r: T = switch (Z) {
        1 => 0.35,
        2 => 0.32,
        3 => 1.45,
        4 => 1.05,
        5 => 0.85,
        6 => 0.65,
        7 => 0.60,
        8 => 0.60,
        9 => 0.50,

        10 => 0.50,
        11 => 1.80,
        12 => 1.50,
        13 => 1.25,
        14 => 1.15,
        15 => 1.00,
        16 => 1.00,
        17 => 1.00,
        18 => 1.00,

        else => 1.50,
    };

    return r * A2BOHR;
}

fn XCInput(comptime T: type) type {
    return struct {
        rho_val: []const T,

        sig_val: ?[]const T = null,
        tau_val: ?[]const T = null,
    };
}

fn XCOutput(comptime T: type) type {
    return struct {
        exc: []T,

        rho_pot: []T,

        sig_pot: ?[]T = null,
        tau_pot: ?[]T = null,
    };
}

fn evaluateXCFunctional(comptime T: type, out: XCOutput(T), func: libxc.xc_func_type, polarized: bool, inp: XCInput(T)) void {
    if (primType(T) != f64) @compileError("FUNCTIONAL EVALUATION NOW ONLY SUPPORTS F64 NUMBERS");

    if (func.info.*.family == libxc.XC_FAMILY_LDA) {
        return evaluateXCLDA(T, out, func, polarized, inp);
    }

    if (func.info.*.family == libxc.XC_FAMILY_GGA) {
        return evaluateXCGGA(T, out, func, polarized, inp);
    }

    if (func.info.*.family == libxc.XC_FAMILY_MGGA) {
        return evaluateXCMGGA(T, out, func, polarized, inp);
    }

    unreachable;
}

fn evaluateXCLDA(comptime T: type, out: XCOutput(T), func: libxc.xc_func_type, polarized: bool, inp: XCInput(T)) void {
    const CHUNK_SIZE = 512;

    var e_xc_tmp: [1 * CHUNK_SIZE]T = undefined;
    var vrho_tmp: [2 * CHUNK_SIZE]T = undefined;

    var i: usize = 0;

    while (i < out.exc.len) : (i += @min(CHUNK_SIZE, out.exc.len - i)) {
        const batch_len = @min(CHUNK_SIZE, out.exc.len - i);

        const exc_size, const vrh_size = .{ batch_len, batch_len * (if (polarized) @as(usize, 2) else 1) };

        const exc_batch = e_xc_tmp[0..exc_size];
        const vrh_batch = vrho_tmp[0..vrh_size];

        @memset(exc_batch, 0);
        @memset(vrh_batch, 0);

        if (polarized) {
            libxc.xc_lda_exc_vxc(&func, batch_len, inp.rho_val[2 * i ..].ptr, exc_batch.ptr, vrh_batch.ptr);

            for (0..batch_len) |j| {
                out.exc[i + j] += exc_batch[j];

                out.rho_pot[2 * (i + j) + 0] += vrh_batch[2 * j + 0];
                out.rho_pot[2 * (i + j) + 1] += vrh_batch[2 * j + 1];
            }
        }

        if (!polarized) {
            libxc.xc_lda_exc_vxc(&func, batch_len, inp.rho_val[i..].ptr, exc_batch.ptr, vrh_batch.ptr);

            for (0..batch_len) |j| {
                out.exc[i + j] += exc_batch[j];

                out.rho_pot[i + j] += vrh_batch[j];
            }
        }
    }
}

fn evaluateXCGGA(comptime T: type, out: XCOutput(T), func: libxc.xc_func_type, polarized: bool, inp: XCInput(T)) void {
    const CHUNK_SIZE = 512;

    const sig_val_unwrapped = inp.sig_val.?;
    const sig_pot_unwrapped = out.sig_pot.?;

    var e_xc_tmp: [1 * CHUNK_SIZE]T = undefined;
    var vrho_tmp: [2 * CHUNK_SIZE]T = undefined;
    var vsig_tmp: [3 * CHUNK_SIZE]T = undefined;

    var i: usize = 0;

    while (i < out.exc.len) : (i += @min(CHUNK_SIZE, out.exc.len - i)) {
        const batch_len = @min(CHUNK_SIZE, out.exc.len - i);

        const vrho_factor: usize = if (polarized) 2 else 1;
        const vsig_factor: usize = if (polarized) 3 else 1;

        const e_xc_size, const vrho_size = .{ batch_len, batch_len * vrho_factor };

        const e_xc_batch = e_xc_tmp[0..e_xc_size];
        const vrho_batch = vrho_tmp[0..vrho_size];

        @memset(e_xc_batch, 0);
        @memset(vrho_batch, 0);

        const vsig_batch = vsig_tmp[0 .. batch_len * vsig_factor];

        @memset(vsig_batch, 0);

        if (polarized) {
            const p1, const p2 = .{ inp.rho_val[2 * i ..].ptr, sig_val_unwrapped[3 * i ..].ptr };

            const p3 = e_xc_batch.ptr;
            const p4 = vrho_batch.ptr;
            const p5 = vsig_batch.ptr;

            libxc.xc_gga_exc_vxc(&func, batch_len, p1, p2, p3, p4, p5);

            for (0..batch_len) |j| {
                out.exc[i + j] += e_xc_batch[j];

                out.rho_pot[2 * (i + j) + 0] += vrho_batch[2 * j + 0];
                out.rho_pot[2 * (i + j) + 1] += vrho_batch[2 * j + 1];

                sig_pot_unwrapped[3 * (i + j) + 0] += vsig_batch[3 * j + 0];
                sig_pot_unwrapped[3 * (i + j) + 1] += vsig_batch[3 * j + 1];
                sig_pot_unwrapped[3 * (i + j) + 2] += vsig_batch[3 * j + 2];
            }
        }

        if (!polarized) {
            const p1, const p2 = .{ inp.rho_val[i..].ptr, sig_val_unwrapped[i..].ptr };

            const p3 = e_xc_batch.ptr;
            const p4 = vrho_batch.ptr;
            const p5 = vsig_batch.ptr;

            libxc.xc_gga_exc_vxc(&func, batch_len, p1, p2, p3, p4, p5);

            for (0..batch_len) |j| {
                out.exc[i + j], out.rho_pot[i + j] = .{ out.exc[i + j] + e_xc_batch[j], out.rho_pot[i + j] + vrho_batch[j] };

                sig_pot_unwrapped[i + j] += vsig_batch[j];
            }
        }
    }
}

fn evaluateXCMGGA(comptime T: type, out: XCOutput(T), func: libxc.xc_func_type, polarized: bool, inp: XCInput(T)) void {
    const CHUNK_SIZE = 512;

    const sig_val_unwrapped = inp.sig_val.?;
    const sig_pot_unwrapped = out.sig_pot.?;
    const tau_val_unwrapped = inp.tau_val.?;
    const tau_pot_unwrapped = out.tau_pot.?;

    var e_xc_tmp: [1 * CHUNK_SIZE]T = undefined;
    var vrho_tmp: [2 * CHUNK_SIZE]T = undefined;
    var vsig_tmp: [3 * CHUNK_SIZE]T = undefined;
    var vtau_tmp: [2 * CHUNK_SIZE]T = undefined;
    var lapl_tmp: [2 * CHUNK_SIZE]T = undefined;
    var vlap_tmp: [2 * CHUNK_SIZE]T = undefined;

    var i: usize = 0;

    while (i < out.exc.len) : (i += @min(CHUNK_SIZE, out.exc.len - i)) {
        const batch_len = @min(CHUNK_SIZE, out.exc.len - i);

        const vrho_factor: usize = if (polarized) 2 else 1;
        const vtau_factor: usize = if (polarized) 2 else 1;
        const vsig_factor: usize = if (polarized) 3 else 1;

        const e_xc_size, const vrho_size = .{ batch_len, batch_len * vrho_factor };

        const vsig_size = batch_len * vsig_factor;
        const vtau_size = batch_len * vtau_factor;

        const e_xc_batch = e_xc_tmp[0..e_xc_size];
        const vrho_batch = vrho_tmp[0..vrho_size];

        @memset(e_xc_batch, 0);
        @memset(vrho_batch, 0);

        const vsig_batch = vsig_tmp[0..vsig_size];
        const vtau_batch = vtau_tmp[0..vtau_size];

        @memset(vsig_batch, 0);
        @memset(vtau_batch, 0);

        const lapl_batch = lapl_tmp[0..vrho_size];
        const vlap_batch = vlap_tmp[0..vrho_size];

        @memset(lapl_batch, 0);
        @memset(vlap_batch, 0);

        if (polarized) {
            const p1 = inp.rho_val[2 * i ..].ptr;

            const p2 = sig_val_unwrapped[3 * i ..].ptr;
            const p4 = tau_val_unwrapped[2 * i ..].ptr;

            const p3 = lapl_batch.ptr;
            const p5 = e_xc_batch.ptr;
            const p6 = vrho_batch.ptr;
            const p7 = vsig_batch.ptr;
            const p8 = vlap_batch.ptr;
            const p9 = vtau_batch.ptr;

            libxc.xc_mgga_exc_vxc(&func, batch_len, p1, p2, p3, p4, p5, p6, p7, p8, p9);

            for (0..batch_len) |j| {
                out.exc[i + j] += e_xc_batch[j];

                out.rho_pot[2 * (i + j) + 0] += vrho_batch[2 * j + 0];
                out.rho_pot[2 * (i + j) + 1] += vrho_batch[2 * j + 1];

                sig_pot_unwrapped[3 * (i + j) + 0] += vsig_batch[3 * j + 0];
                sig_pot_unwrapped[3 * (i + j) + 1] += vsig_batch[3 * j + 1];
                sig_pot_unwrapped[3 * (i + j) + 2] += vsig_batch[3 * j + 2];

                tau_pot_unwrapped[2 * (i + j) + 0] += vtau_batch[2 * j + 0];
                tau_pot_unwrapped[2 * (i + j) + 1] += vtau_batch[2 * j + 1];
            }
        }

        if (!polarized) {
            const p1 = inp.rho_val[i..].ptr;

            const p2 = sig_val_unwrapped[i..].ptr;
            const p4 = tau_val_unwrapped[i..].ptr;

            const p3 = lapl_batch.ptr;
            const p5 = e_xc_batch.ptr;
            const p6 = vrho_batch.ptr;
            const p7 = vsig_batch.ptr;
            const p8 = vlap_batch.ptr;
            const p9 = vtau_batch.ptr;

            libxc.xc_mgga_exc_vxc(&func, batch_len, p1, p2, p3, p4, p5, p6, p7, p8, p9);

            for (0..batch_len) |j| {
                out.exc[i + j], out.rho_pot[i + j] = .{ out.exc[i + j] + e_xc_batch[j], out.rho_pot[i + j] + vrho_batch[j] };

                sig_pot_unwrapped[i + j] += vsig_batch[j];
                tau_pot_unwrapped[i + j] += vtau_batch[j];
            }
        }
    }
}
