const std = @import("std");
const libint = @import("libint");

const libxc = @cImport(@cInclude("xc.h"));

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Vector = @import("tensor.zig").Vector;

const getLebedevGrid = @import("lebedev_quadrature_nodes.zig").getLebedevGrid;
const primType = @import("value.zig").primType;

const A2BOHR = @import("constant.zig").A2BOHR;

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
        exx_coef: T = 0,

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
                exch_func = try initFunctional(name, mode, libxc.XC_EXCHANGE, gpa);
            }

            errdefer {
                if (exch_func) |*f| libxc.xc_func_end(f);
            }

            if (funcs[1]) |name| {
                corr_func = try initFunctional(name, mode, libxc.XC_CORRELATION, gpa);
            }

            errdefer {
                if (corr_func) |*f| libxc.xc_func_end(f);
            }

            if (funcs[2]) |name| {
                exco_func = try initFunctional(name, mode, libxc.XC_EXCHANGE_CORRELATION, gpa);
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

            var exx: T = 0.0;

            if (exch_func) |*func| {
                exx = @floatCast(libxc.xc_hyb_exx_coef(func));
            }

            if (exco_func) |*func| {
                exx = @floatCast(libxc.xc_hyb_exx_coef(func));
            }

            var Vxc = try Matrix(T).init(nbf, nbf, gpa);
            errdefer Vxc.deinit(gpa);

            const fcts = .{ exch_func, corr_func, exco_func };

            return .{ .exch = fcts[0], .corr = fcts[1], .exco = fcts[2], .pts = grid[0], .wgh = grid[1], .Vxc = Vxc, .polarized = polarized, .exx_coef = exx };
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

            const size_factor: usize, const flags = .{ if (self.polarized) 2 else 1, self.getFamilyFlags() };

            var basis = try BasisGrid(T).init(self.pts.shape[0], sys.nbf, flags.needs_deriv, flags.has_mgga, gpa);
            defer basis.deinit(gpa);

            var density = try DensityGrid(T).init(self.pts.shape[0], size_factor, flags.needs_deriv, flags.has_mgga, self.polarized, gpa);
            defer density.deinit(gpa);

            var potential = try PotentialGrid(T).init(self.pts.shape[0], size_factor, flags.needs_deriv, flags.has_mgga, self.polarized, gpa);
            defer potential.deinit(gpa);

            self.evaluateBasis(sys, basis);

            var ctx = DftEvaluationContext(T){
                .sys = sys,

                .basisgrid = basis,
                .density = density,
                .potgd = potential,

                .has_sgga = flags.has_sgga,
                .has_mgga = flags.has_mgga,

                .needs_deriv = flags.needs_deriv,
            };

            self.computeDensity(P, &ctx);

            const inp = XCInput(T){
                .rho_val = ctx.density.rho_val.data,

                .sig_val = if (ctx.density.sig_val) |s| s.data else null,
                .tau_val = if (ctx.density.tau_val) |t| t.data else null,
                .lap_val = if (ctx.density.lap_val) |l| l.data else null,
            };

            const out = XCOutput(T){
                .exc = ctx.potgd.exc.data,

                .sig_pot = if (ctx.potgd.sig_pot) |s| s.data else null,
                .tau_pot = if (ctx.potgd.tau_pot) |t| t.data else null,
                .lap_pot = if (ctx.potgd.lap_pot) |l| l.data else null,

                .rho_pot = ctx.potgd.rho_pot.data,
            };

            if (self.exch) |func| evaluateXCFunctional(T, out, func, self.polarized, inp);
            if (self.corr) |func| evaluateXCFunctional(T, out, func, self.polarized, inp);
            if (self.exco) |func| evaluateXCFunctional(T, out, func, self.polarized, inp);

            try self.integratePotential(ctx, gpa);
        }

        pub fn getFunctionalNames(self: @This(), gpa: Allocator) ![]const u8 {
            if (self.exco) |func| {
                return try gpa.dupe(u8, std.mem.span(func.info.*.name));
            }

            var parts = std.ArrayList([]const u8).empty;
            defer parts.deinit(gpa);

            if (self.exch) |func| {
                try parts.append(gpa, std.mem.span(func.info.*.name));
            }

            if (self.corr) |func| {
                try parts.append(gpa, std.mem.span(func.info.*.name));
            }

            return try std.mem.join(gpa, " + ", parts.items);
        }

        fn computeDensity(self: @This(), P: Matrix(T), ctx: *DftEvaluationContext(T)) void {
            if (self.polarized) {
                self.computeDensityPolarized(P, ctx);
            } else {
                self.computeDensityUnpolarized(P, ctx);
            }
        }

        fn computeDensityPolarized(self: @This(), P: Matrix(T), ctx: *DftEvaluationContext(T)) void {
            for (0..self.pts.shape[0]) |i| {
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

                var lapl_sum_a: T = 0;
                var lapl_sum_b: T = 0;

                for (0..ctx.sys.nbf) |mu| for (0..ctx.sys.nbf) |nu| {
                    const b = ctx.sys.nbf;

                    const p_val_a = P.at(mu + 0, nu + 0);
                    const p_val_b = P.at(mu + b, nu + b);

                    const phi_mu = ctx.basisgrid.phi.at(i, mu);
                    const phi_nu = ctx.basisgrid.phi.at(i, nu);

                    rho_a += p_val_a * phi_mu * phi_nu;
                    rho_b += p_val_b * phi_mu * phi_nu;

                    if (ctx.needs_deriv) {
                        const dphi_dx_mu = ctx.basisgrid.dphi_dx.?.at(i, mu);
                        const dphi_dy_mu = ctx.basisgrid.dphi_dy.?.at(i, mu);
                        const dphi_dz_mu = ctx.basisgrid.dphi_dz.?.at(i, mu);

                        const term_a = p_val_a * phi_nu;
                        const term_b = p_val_b * phi_nu;

                        del_rho_a_x += term_a * dphi_dx_mu;
                        del_rho_a_y += term_a * dphi_dy_mu;
                        del_rho_a_z += term_a * dphi_dz_mu;

                        del_rho_b_x += term_b * dphi_dx_mu;
                        del_rho_b_y += term_b * dphi_dy_mu;
                        del_rho_b_z += term_b * dphi_dz_mu;

                        if (ctx.has_mgga) {
                            const dphi_dx_nu = ctx.basisgrid.dphi_dx.?.at(i, nu);
                            const dphi_dy_nu = ctx.basisgrid.dphi_dy.?.at(i, nu);
                            const dphi_dz_nu = ctx.basisgrid.dphi_dz.?.at(i, nu);

                            const dphi_dot = dphi_dx_mu * dphi_dx_nu + dphi_dy_mu * dphi_dy_nu + dphi_dz_mu * dphi_dz_nu;

                            tau_a += p_val_a * dphi_dot;
                            tau_b += p_val_b * dphi_dot;

                            const lapl_mu = ctx.basisgrid.lapl.?.at(i, mu);

                            lapl_sum_a += p_val_a * lapl_mu * phi_nu;
                            lapl_sum_b += p_val_b * lapl_mu * phi_nu;
                        }
                    }
                };

                ctx.density.rho_val.ptr(i, 0).* = rho_a;
                ctx.density.rho_val.ptr(i, 1).* = rho_b;

                if (ctx.needs_deriv) {
                    del_rho_a_x *= 2;
                    del_rho_a_y *= 2;
                    del_rho_a_z *= 2;
                    del_rho_b_x *= 2;
                    del_rho_b_y *= 2;
                    del_rho_b_z *= 2;

                    ctx.density.del_rho.?.ptr(i, 0).* = del_rho_a_x;
                    ctx.density.del_rho.?.ptr(i, 1).* = del_rho_a_y;
                    ctx.density.del_rho.?.ptr(i, 2).* = del_rho_a_z;
                    ctx.density.del_rho.?.ptr(i, 3).* = del_rho_b_x;
                    ctx.density.del_rho.?.ptr(i, 4).* = del_rho_b_y;
                    ctx.density.del_rho.?.ptr(i, 5).* = del_rho_b_z;

                    ctx.density.sig_val.?.ptr(i, 0).* = del_rho_a_x * del_rho_a_x + del_rho_a_y * del_rho_a_y + del_rho_a_z * del_rho_a_z;
                    ctx.density.sig_val.?.ptr(i, 1).* = del_rho_a_x * del_rho_b_x + del_rho_a_y * del_rho_b_y + del_rho_a_z * del_rho_b_z;
                    ctx.density.sig_val.?.ptr(i, 2).* = del_rho_b_x * del_rho_b_x + del_rho_b_y * del_rho_b_y + del_rho_b_z * del_rho_b_z;
                }

                if (ctx.has_mgga) {
                    ctx.density.tau_val.?.ptr(i, 0).* = 0.5 * tau_a;
                    ctx.density.tau_val.?.ptr(i, 1).* = 0.5 * tau_b;

                    ctx.density.lap_val.?.ptr(i, 0).* = 2 * (lapl_sum_a + tau_a);
                    ctx.density.lap_val.?.ptr(i, 1).* = 2 * (lapl_sum_b + tau_b);
                }
            }
        }

        fn computeDensityUnpolarized(self: @This(), P: Matrix(T), ctx: *DftEvaluationContext(T)) void {
            for (0..self.pts.shape[0]) |i| {
                var rho_tot: T = 0;
                var tau_tot: T = 0;

                var lapl_sum_tot: T = 0;

                var del_rho_x: T = 0;
                var del_rho_y: T = 0;
                var del_rho_z: T = 0;

                for (0..ctx.sys.nbf) |mu| for (0..ctx.sys.nbf) |nu| {
                    const phi_mu = ctx.basisgrid.phi.at(i, mu);
                    const phi_nu = ctx.basisgrid.phi.at(i, nu);

                    rho_tot += P.at(mu, nu) * phi_mu * phi_nu;

                    if (ctx.needs_deriv) {
                        const dphi_dx_mu = ctx.basisgrid.dphi_dx.?.at(i, mu);
                        const dphi_dy_mu = ctx.basisgrid.dphi_dy.?.at(i, mu);
                        const dphi_dz_mu = ctx.basisgrid.dphi_dz.?.at(i, mu);

                        const term = P.at(mu, nu) * phi_nu;

                        del_rho_x += term * dphi_dx_mu;
                        del_rho_y += term * dphi_dy_mu;
                        del_rho_z += term * dphi_dz_mu;

                        if (ctx.has_mgga) {
                            const dphi_dx_nu = ctx.basisgrid.dphi_dx.?.at(i, nu);
                            const dphi_dy_nu = ctx.basisgrid.dphi_dy.?.at(i, nu);
                            const dphi_dz_nu = ctx.basisgrid.dphi_dz.?.at(i, nu);

                            const dphi_dot = dphi_dx_mu * dphi_dx_nu + dphi_dy_mu * dphi_dy_nu + dphi_dz_mu * dphi_dz_nu;

                            tau_tot += P.at(mu, nu) * dphi_dot;

                            lapl_sum_tot += P.at(mu, nu) * ctx.basisgrid.lapl.?.at(i, mu) * phi_nu;
                        }
                    }
                };

                ctx.density.rho_val.ptr(i, 0).* = rho_tot;

                if (ctx.needs_deriv) {
                    del_rho_x *= 2;
                    del_rho_y *= 2;
                    del_rho_z *= 2;

                    ctx.density.del_rho.?.ptr(i, 0).* = del_rho_x;
                    ctx.density.del_rho.?.ptr(i, 1).* = del_rho_y;
                    ctx.density.del_rho.?.ptr(i, 2).* = del_rho_z;

                    ctx.density.sig_val.?.ptr(i, 0).* = del_rho_x * del_rho_x + del_rho_y * del_rho_y + del_rho_z * del_rho_z;
                }

                if (ctx.has_mgga) {
                    ctx.density.tau_val.?.ptr(i, 0).* = 0.5 * tau_tot;

                    ctx.density.lap_val.?.ptr(i, 0).* = 2 * (lapl_sum_tot + tau_tot);
                }
            }
        }

        fn evaluateBasis(self: @This(), sys: MolecularSystem(T), basis: BasisGrid(T)) void {
            const needs_d1, const needs_d2 = .{ basis.dphi_dx != null, basis.lapl != null };

            for (0..self.pts.shape[0]) |i| {
                const x = self.pts.at(i, 0);
                const y = self.pts.at(i, 1);
                const z = self.pts.at(i, 2);

                if (needs_d2) {
                    const dx_row = basis.dphi_dx.?.rowSlice(i);
                    const dy_row = basis.dphi_dy.?.rowSlice(i);
                    const dz_row = basis.dphi_dz.?.rowSlice(i);

                    const lapl_row = basis.lapl.?.rowSlice(i);

                    libint.evaluate_basis_vgl(basis.phi.rowSlice(i).ptr, dx_row.ptr, dy_row.ptr, dz_row.ptr, lapl_row.ptr, x, y, z, sys.ptr);
                }

                if (needs_d1 and !needs_d2) {
                    const dx_row = basis.dphi_dx.?.rowSlice(i);
                    const dy_row = basis.dphi_dy.?.rowSlice(i);
                    const dz_row = basis.dphi_dz.?.rowSlice(i);

                    libint.evaluate_basis_vg(basis.phi.rowSlice(i).ptr, dx_row.ptr, dy_row.ptr, dz_row.ptr, x, y, z, sys.ptr);
                }

                if (!needs_d1 and !needs_d2) {
                    libint.evaluate_basis_v(basis.phi.rowSlice(i).ptr, x, y, z, sys.ptr);
                }
            }
        }

        fn getFamilyFlags(self: @This()) struct { has_sgga: bool, has_mgga: bool, needs_deriv: bool } {
            const exch_fam = if (self.exch) |func| func.info.*.family else null;
            const corr_fam = if (self.corr) |func| func.info.*.family else null;
            const exco_fam = if (self.exco) |func| func.info.*.family else null;

            const is_exch_sgga = self.exch != null and (exch_fam == libxc.XC_FAMILY_GGA or exch_fam == libxc.XC_FAMILY_HYB_GGA);
            const is_corr_sgga = self.corr != null and (corr_fam == libxc.XC_FAMILY_GGA or corr_fam == libxc.XC_FAMILY_HYB_GGA);
            const is_exco_sgga = self.exco != null and (exco_fam == libxc.XC_FAMILY_GGA or exco_fam == libxc.XC_FAMILY_HYB_GGA);

            const is_exch_mgga = self.exch != null and (exch_fam == libxc.XC_FAMILY_MGGA or exch_fam == libxc.XC_FAMILY_HYB_MGGA);
            const is_corr_mgga = self.corr != null and (corr_fam == libxc.XC_FAMILY_MGGA or corr_fam == libxc.XC_FAMILY_HYB_MGGA);
            const is_exco_mgga = self.exco != null and (exco_fam == libxc.XC_FAMILY_MGGA or exco_fam == libxc.XC_FAMILY_HYB_MGGA);

            const has_sgga = is_exch_sgga or is_corr_sgga or is_exco_sgga;
            const has_mgga = is_exch_mgga or is_corr_mgga or is_exco_mgga;

            return .{ .has_sgga = has_sgga, .has_mgga = has_mgga, .needs_deriv = has_sgga or has_mgga };
        }

        fn initFunctional(name: []const u8, mode: c_int, expected_kind: c_int, gpa: Allocator) !libxc.xc_func_type {
            const name_nt = try gpa.dupeSentinel(u8, name, 0);
            defer gpa.free(name_nt);

            const id = libxc.xc_functional_get_number(name_nt.ptr);
            if (id < 0) {
                return switch (expected_kind) {
                    libxc.XC_EXCHANGE => error.ExchFuncNotFound,
                    libxc.XC_CORRELATION => error.CorrFuncNotFound,
                    else => error.XcFuncNotFound,
                };
            }

            var func: libxc.xc_func_type = undefined;
            if (libxc.xc_func_init(&func, id, mode) != 0) {
                return error.LibxcInitFailed;
            }

            if (func.info.*.kind != expected_kind) {
                libxc.xc_func_end(&func);
                return switch (expected_kind) {
                    libxc.XC_EXCHANGE => error.InvalidExchangeFunctional,
                    libxc.XC_CORRELATION => error.InvalidCorrelationFunctional,
                    else => error.InvalidExchangeCorrelationFunctional,
                };
            }

            return func;
        }

        fn integratePotential(self: *@This(), ctx: DftEvaluationContext(T), gpa: Allocator) !void {
            if (self.polarized) {
                self.Exc = try self.integratePotentialPolarized(ctx, gpa);
            }

            if (!self.polarized) {
                self.Exc = try self.integratePotentialUnpolarized(ctx, gpa);
            }
        }

        fn integratePotentialPolarized(self: *@This(), ctx: DftEvaluationContext(T), gpa: Allocator) !T {
            var energy_exc: T = 0;

            var X_a: ?Vector(T) = null;
            var X_b: ?Vector(T) = null;

            if (ctx.needs_deriv) {
                X_a = try Vector(T).init(ctx.sys.nbf, gpa);
                X_b = try Vector(T).init(ctx.sys.nbf, gpa);
            }

            defer {
                if (X_a) |*x| x.deinit(gpa);
                if (X_b) |*x| x.deinit(gpa);
            }

            for (0..self.pts.shape[0]) |i| {
                const w = self.wgh.at(i);

                energy_exc += (ctx.density.rho_val.at(i, 0) + ctx.density.rho_val.at(i, 1)) * ctx.potgd.exc.at(i) * w;

                const v_a = ctx.potgd.rho_pot.at(i, 0);
                const v_b = ctx.potgd.rho_pot.at(i, 1);

                if (ctx.needs_deriv) {
                    const g_a_x = ctx.density.del_rho.?.at(i, 0);
                    const g_a_y = ctx.density.del_rho.?.at(i, 1);
                    const g_a_z = ctx.density.del_rho.?.at(i, 2);

                    const g_b_x = ctx.density.del_rho.?.at(i, 3);
                    const g_b_y = ctx.density.del_rho.?.at(i, 4);
                    const g_b_z = ctx.density.del_rho.?.at(i, 5);

                    const v_sig_aa = ctx.potgd.sig_pot.?.at(i, 0);
                    const v_sig_ab = ctx.potgd.sig_pot.?.at(i, 1);
                    const v_sig_bb = ctx.potgd.sig_pot.?.at(i, 2);

                    const w_a_x = 2 * v_sig_aa * g_a_x + v_sig_ab * g_b_x;
                    const w_a_y = 2 * v_sig_aa * g_a_y + v_sig_ab * g_b_y;
                    const w_a_z = 2 * v_sig_aa * g_a_z + v_sig_ab * g_b_z;

                    const w_b_x = 2 * v_sig_bb * g_b_x + v_sig_ab * g_a_x;
                    const w_b_y = 2 * v_sig_bb * g_b_y + v_sig_ab * g_a_y;
                    const w_b_z = 2 * v_sig_bb * g_b_z + v_sig_ab * g_a_z;

                    for (0..ctx.sys.nbf) |mu| {
                        const d_dx = ctx.basisgrid.dphi_dx.?.at(i, mu);
                        const d_dy = ctx.basisgrid.dphi_dy.?.at(i, mu);
                        const d_dz = ctx.basisgrid.dphi_dz.?.at(i, mu);

                        X_a.?.data[mu] = w_a_x * d_dx + w_a_y * d_dy + w_a_z * d_dz;
                        X_b.?.data[mu] = w_b_x * d_dx + w_b_y * d_dy + w_b_z * d_dz;
                    }

                    const v_tau_a = if (ctx.has_mgga) ctx.potgd.tau_pot.?.at(i, 0) else 0.0;
                    const v_tau_b = if (ctx.has_mgga) ctx.potgd.tau_pot.?.at(i, 1) else 0.0;
                    const v_lap_a = if (ctx.has_mgga) ctx.potgd.lap_pot.?.at(i, 0) else 0.0;
                    const v_lap_b = if (ctx.has_mgga) ctx.potgd.lap_pot.?.at(i, 1) else 0.0;

                    for (0..ctx.sys.nbf) |mu| for (0..ctx.sys.nbf) |nu| {
                        const phi_mu = ctx.basisgrid.phi.at(i, mu);
                        const phi_nu = ctx.basisgrid.phi.at(i, nu);

                        const b = ctx.sys.nbf;

                        var val_a = v_a * phi_mu * phi_nu + X_a.?.data[mu] * phi_nu + phi_mu * X_a.?.data[nu];
                        var val_b = v_b * phi_mu * phi_nu + X_b.?.data[mu] * phi_nu + phi_mu * X_b.?.data[nu];

                        if (ctx.has_mgga) {
                            const dx_mu_nu = ctx.basisgrid.dphi_dx.?.at(i, mu) * ctx.basisgrid.dphi_dx.?.at(i, nu);
                            const dy_mu_nu = ctx.basisgrid.dphi_dy.?.at(i, mu) * ctx.basisgrid.dphi_dy.?.at(i, nu);
                            const dz_mu_nu = ctx.basisgrid.dphi_dz.?.at(i, mu) * ctx.basisgrid.dphi_dz.?.at(i, nu);

                            const dphi_dot = dx_mu_nu + dy_mu_nu + dz_mu_nu;

                            val_a += 0.5 * v_tau_a * dphi_dot;
                            val_b += 0.5 * v_tau_b * dphi_dot;

                            const lapl_mu = ctx.basisgrid.lapl.?.at(i, mu);
                            const lapl_nu = ctx.basisgrid.lapl.?.at(i, nu);

                            const lapl_term = lapl_mu * phi_nu + 2 * dphi_dot + phi_mu * lapl_nu;

                            val_a += v_lap_a * lapl_term;
                            val_b += v_lap_b * lapl_term;
                        }

                        self.Vxc.ptr(mu + 0, nu + 0).* += val_a * w;
                        self.Vxc.ptr(mu + b, nu + b).* += val_b * w;
                    };
                }

                if (!ctx.needs_deriv) {
                    for (0..ctx.sys.nbf) |mu| for (0..ctx.sys.nbf) |nu| {
                        const phi_mu = ctx.basisgrid.phi.at(i, mu);
                        const phi_nu = ctx.basisgrid.phi.at(i, nu);

                        const b = ctx.sys.nbf;

                        self.Vxc.ptr(mu + 0, nu + 0).* += v_a * phi_mu * phi_nu * w;
                        self.Vxc.ptr(mu + b, nu + b).* += v_b * phi_mu * phi_nu * w;
                    };
                }
            }
            return energy_exc;
        }

        fn integratePotentialUnpolarized(self: *@This(), ctx: DftEvaluationContext(T), gpa: Allocator) !T {
            var energy_exc: T, var X_a: ?Vector(T) = .{ 0, null };

            if (ctx.needs_deriv) {
                X_a = try Vector(T).init(ctx.sys.nbf, gpa);
            }

            defer {
                if (X_a) |*x| x.deinit(gpa);
            }

            for (0..self.pts.shape[0]) |i| {
                const w = self.wgh.at(i);

                energy_exc += ctx.density.rho_val.at(i, 0) * ctx.potgd.exc.at(i) * w;

                if (ctx.needs_deriv) {
                    const v_sigma = ctx.potgd.sig_pot.?.at(i, 0);

                    const w_x = 2 * v_sigma * ctx.density.del_rho.?.at(i, 0);
                    const w_y = 2 * v_sigma * ctx.density.del_rho.?.at(i, 1);
                    const w_z = 2 * v_sigma * ctx.density.del_rho.?.at(i, 2);

                    for (0..ctx.sys.nbf) |mu| {
                        const dx_mu = ctx.basisgrid.dphi_dx.?.at(i, mu);
                        const dy_mu = ctx.basisgrid.dphi_dy.?.at(i, mu);
                        const dz_mu = ctx.basisgrid.dphi_dz.?.at(i, mu);

                        X_a.?.data[mu] = w_x * dx_mu + w_y * dy_mu + w_z * dz_mu;
                    }

                    const v_tau = if (ctx.has_mgga) ctx.potgd.tau_pot.?.at(i, 0) else 0.0;
                    const v_lap = if (ctx.has_mgga) ctx.potgd.lap_pot.?.at(i, 0) else 0.0;

                    for (0..ctx.sys.nbf) |mu| for (0..ctx.sys.nbf) |nu| {
                        const phi_mu = ctx.basisgrid.phi.at(i, mu);
                        const phi_nu = ctx.basisgrid.phi.at(i, nu);

                        var val = ctx.potgd.rho_pot.at(i, 0) * phi_mu * phi_nu + X_a.?.data[mu] * phi_nu + phi_mu * X_a.?.data[nu];

                        if (ctx.has_mgga) {
                            const dx_mu_nu = ctx.basisgrid.dphi_dx.?.at(i, mu) * ctx.basisgrid.dphi_dx.?.at(i, nu);
                            const dy_mu_nu = ctx.basisgrid.dphi_dy.?.at(i, mu) * ctx.basisgrid.dphi_dy.?.at(i, nu);
                            const dz_mu_nu = ctx.basisgrid.dphi_dz.?.at(i, mu) * ctx.basisgrid.dphi_dz.?.at(i, nu);

                            const dphi_dot = dx_mu_nu + dy_mu_nu + dz_mu_nu;

                            val += 0.5 * v_tau * dphi_dot;

                            const lapl_mu = ctx.basisgrid.lapl.?.at(i, mu);
                            const lapl_nu = ctx.basisgrid.lapl.?.at(i, nu);
                            const lapl_term = lapl_mu * phi_nu + 2 * dphi_dot + phi_mu * lapl_nu;

                            val += v_lap * lapl_term;
                        }

                        self.Vxc.ptr(mu, nu).* += val * w;
                    };
                }

                if (!ctx.needs_deriv) {
                    for (0..ctx.sys.nbf) |mu| for (0..ctx.sys.nbf) |nu| {
                        self.Vxc.ptr(mu, nu).* += ctx.potgd.rho_pot.at(i, 0) * ctx.basisgrid.phi.at(i, mu) * ctx.basisgrid.phi.at(i, nu) * w;
                    };
                }
            }
            return energy_exc;
        }
    };
}

fn BasisGrid(comptime T: type) type {
    return struct {
        phi: Matrix(T),

        dphi_dx: ?Matrix(T) = null,
        dphi_dy: ?Matrix(T) = null,
        dphi_dz: ?Matrix(T) = null,

        lapl: ?Matrix(T) = null,

        fn init(n_grid: usize, nbf: usize, needs_deriv: bool, has_mgga: bool, gpa: Allocator) !@This() {
            var phi = try Matrix(T).init(n_grid, nbf, gpa);
            errdefer phi.deinit(gpa);

            var dphi_dx: ?Matrix(T) = null;
            var dphi_dy: ?Matrix(T) = null;
            var dphi_dz: ?Matrix(T) = null;

            var lapl: ?Matrix(T) = null;

            if (needs_deriv) {
                dphi_dx = try Matrix(T).init(n_grid, nbf, gpa);
                errdefer if (dphi_dx) |*d| d.deinit(gpa);

                dphi_dy = try Matrix(T).init(n_grid, nbf, gpa);
                errdefer if (dphi_dy) |*d| d.deinit(gpa);

                dphi_dz = try Matrix(T).init(n_grid, nbf, gpa);
                errdefer if (dphi_dz) |*d| d.deinit(gpa);
            }

            if (has_mgga) {
                lapl = try Matrix(T).init(n_grid, nbf, gpa);
            }

            return .{ .phi = phi, .dphi_dx = dphi_dx, .dphi_dy = dphi_dy, .dphi_dz = dphi_dz, .lapl = lapl };
        }

        fn deinit(self: *@This(), gpa: Allocator) void {
            self.phi.deinit(gpa);

            if (self.dphi_dx) |*d| d.deinit(gpa);
            if (self.dphi_dy) |*d| d.deinit(gpa);
            if (self.dphi_dz) |*d| d.deinit(gpa);

            if (self.lapl) |*d| d.deinit(gpa);
        }
    };
}

fn Becke(comptime T: type) type {
    return struct {
        n_rad: usize,
        n_leb: usize,

        fn beckeStep(mu: T) T {
            var p = 1.5 * mu - 0.5 * mu * mu * mu;

            p = 1.5 * p - 0.5 * p * p * p;
            p = 1.5 * p - 0.5 * p * p * p;

            return 0.5 * (1.0 - p);
        }

        fn get(self: @This(), sys: MolecularSystem(T), gpa: Allocator) !struct { Matrix(T), Vector(T) } {
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

fn DensityGrid(comptime T: type) type {
    return struct {
        rho_val: Matrix(T),

        del_rho: ?Matrix(T) = null,
        sig_val: ?Matrix(T) = null,
        tau_val: ?Matrix(T) = null,
        lap_val: ?Matrix(T) = null,

        fn init(n_g: usize, size_factor: usize, needs_deriv: bool, has_mgga: bool, pol: bool, gpa: Allocator) !@This() {
            var rho_val = try Matrix(T).initZero(n_g, size_factor, gpa);
            errdefer rho_val.deinit(gpa);

            var del_rho: ?Matrix(T) = null;
            var sig_val: ?Matrix(T) = null;
            var tau_val: ?Matrix(T) = null;
            var lap_val: ?Matrix(T) = null;

            if (needs_deriv) {
                sig_val = try Matrix(T).initZero(n_g, if (pol) 3 else 1, gpa);
                errdefer if (sig_val) |*s| s.deinit(gpa);

                del_rho = try Matrix(T).initZero(n_g, if (pol) 6 else 3, gpa);
                errdefer if (del_rho) |*g| g.deinit(gpa);
            }

            if (has_mgga) {
                tau_val = try Matrix(T).initZero(n_g, size_factor, gpa);
                errdefer if (tau_val) |*t| t.deinit(gpa);

                lap_val = try Matrix(T).initZero(n_g, size_factor, gpa);
                errdefer if (lap_val) |*l| l.deinit(gpa);
            }

            return .{ .rho_val = rho_val, .del_rho = del_rho, .sig_val = sig_val, .tau_val = tau_val, .lap_val = lap_val };
        }

        fn deinit(self: *@This(), gpa: Allocator) void {
            self.rho_val.deinit(gpa);

            if (self.del_rho) |*g| g.deinit(gpa);
            if (self.sig_val) |*s| s.deinit(gpa);
            if (self.tau_val) |*t| t.deinit(gpa);
            if (self.lap_val) |*l| l.deinit(gpa);
        }
    };
}

fn DftEvaluationContext(comptime T: type) type {
    return struct {
        sys: MolecularSystem(T),
        density: DensityGrid(T),
        basisgrid: BasisGrid(T),
        potgd: PotentialGrid(T),

        has_sgga: bool,
        has_mgga: bool,

        needs_deriv: bool,
    };
}

fn PotentialGrid(comptime T: type) type {
    return struct {
        exc: Vector(T),

        rho_pot: Matrix(T),

        sig_pot: ?Matrix(T) = null,
        tau_pot: ?Matrix(T) = null,
        lap_pot: ?Matrix(T) = null,

        fn init(n_g: usize, size_factor: usize, needs_deriv: bool, has_mgga: bool, pol: bool, gpa: Allocator) !@This() {
            var exc = try Vector(T).initZero(n_g, gpa);
            errdefer exc.deinit(gpa);

            var rho_pot = try Matrix(T).initZero(n_g, size_factor, gpa);
            errdefer rho_pot.deinit(gpa);

            var sig_pot: ?Matrix(T) = null;
            var tau_pot: ?Matrix(T) = null;
            var lap_pot: ?Matrix(T) = null;

            if (needs_deriv) {
                sig_pot = try Matrix(T).initZero(n_g, if (pol) 3 else 1, gpa);
                errdefer if (sig_pot) |*s| s.deinit(gpa);
            }

            if (has_mgga) {
                tau_pot = try Matrix(T).initZero(n_g, size_factor, gpa);
                errdefer if (tau_pot) |*t| t.deinit(gpa);

                lap_pot = try Matrix(T).initZero(n_g, size_factor, gpa);
                errdefer if (lap_pot) |*l| l.deinit(gpa);
            }

            return .{ .exc = exc, .rho_pot = rho_pot, .sig_pot = sig_pot, .tau_pot = tau_pot, .lap_pot = lap_pot };
        }

        fn deinit(self: *@This(), gpa: Allocator) void {
            self.exc.deinit(gpa);

            self.rho_pot.deinit(gpa);

            if (self.sig_pot) |*s| s.deinit(gpa);
            if (self.tau_pot) |*t| t.deinit(gpa);
            if (self.lap_pot) |*l| l.deinit(gpa);
        }
    };
}

fn XCInput(comptime T: type) type {
    return struct {
        rho_val: []const T,

        sig_val: ?[]const T = null,
        tau_val: ?[]const T = null,
        lap_val: ?[]const T = null,
    };
}

fn XCOutput(comptime T: type) type {
    return struct {
        exc: []T,

        rho_pot: []T,

        sig_pot: ?[]T = null,
        tau_pot: ?[]T = null,
        lap_pot: ?[]T = null,
    };
}

fn evaluateXCFunctional(comptime T: type, out: XCOutput(T), func: libxc.xc_func_type, polarized: bool, inp: XCInput(T)) void {
    if (primType(T) != f64) @compileError("FUNCTIONAL EVALUATION NOW ONLY SUPPORTS F64 NUMBERS");

    if (func.info.*.family == libxc.XC_FAMILY_LDA or func.info.*.family == libxc.XC_FAMILY_HYB_LDA) {
        return evaluateXCLDA(T, out, func, polarized, inp);
    }

    if (func.info.*.family == libxc.XC_FAMILY_GGA or func.info.*.family == libxc.XC_FAMILY_HYB_GGA) {
        return evaluateXCGGA(T, out, func, polarized, inp);
    }

    if (func.info.*.family == libxc.XC_FAMILY_MGGA or func.info.*.family == libxc.XC_FAMILY_HYB_MGGA) {
        return evaluateXCMGGA(T, out, func, polarized, inp);
    }

    unreachable;
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

fn evaluateXCMGGA(comptime T: type, out: XCOutput(T), func: libxc.xc_func_type, polarized: bool, inp: XCInput(T)) void {
    const CHUNK_SIZE = 512;

    const sig_val_unwrapped = inp.sig_val.?;
    const sig_pot_unwrapped = out.sig_pot.?;
    const tau_val_unwrapped = inp.tau_val.?;
    const tau_pot_unwrapped = out.tau_pot.?;
    const lap_val_unwrapped = inp.lap_val.?;
    const lap_pot_unwrapped = out.lap_pot.?;

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

        @memset(vlap_batch, 0);

        if (polarized) {
            const p1 = inp.rho_val[2 * i ..].ptr;

            const p2 = sig_val_unwrapped[3 * i ..].ptr;
            const p4 = tau_val_unwrapped[2 * i ..].ptr;

            @memcpy(lapl_batch, lap_val_unwrapped[2 * i .. 2 * i + vrho_size]);

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
                lap_pot_unwrapped[2 * (i + j) + 0] += vlap_batch[2 * j + 0];
                lap_pot_unwrapped[2 * (i + j) + 1] += vlap_batch[2 * j + 1];
            }
        }

        if (!polarized) {
            const p1 = inp.rho_val[i..].ptr;

            const p2 = sig_val_unwrapped[i..].ptr;
            const p4 = tau_val_unwrapped[i..].ptr;

            @memcpy(lapl_batch, lap_val_unwrapped[i .. i + vrho_size]);

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
                lap_pot_unwrapped[i + j] += vlap_batch[j];
            }
        }
    }
}

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
