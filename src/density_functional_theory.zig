//! Implements DFT numerical integration grids and exchange-correlation potential evaluation.

const std = @import("std");

const libint = @import("cimport.zig").libint;
const libxc = @import("cimport.zig").libxc;

const Allocator = std.mem.Allocator;

const BasisGrid = @import("molecular_grid.zig").BasisGrid;
const Becke = @import("molecular_grid.zig").Becke;
const DensityGrid = @import("molecular_grid.zig").DensityGrid;
const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const PotentialGrid = @import("molecular_grid.zig").PotentialGrid;
const Vector = @import("tensor.zig").Vector;
const XCInput = @import("xc_functional.zig").XCInput;
const XCOutput = @import("xc_functional.zig").XCOutput;

const evaluateXCFunctional = @import("xc_functional.zig").evaluateXCFunctional;

/// Constructs a struct representing the DFT exchange-correlation potential and numerical integration grid.
pub fn DftPotential(comptime T: type) type {
    return struct {
        exch: ?libxc.xc_func_type,
        corr: ?libxc.xc_func_type,
        exco: ?libxc.xc_func_type,

        a_x: T = 0,

        pts: Matrix(T),
        wgh: Vector(T),

        V_xc: Matrix(T),
        E_xc: T = 0.000,

        polarized: bool,

        /// Initializes the DFT potential, setting up the molecular grid and initializing libxc functionals.
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

            var V_xc = try Matrix(T).init(nbf, nbf, gpa);
            errdefer V_xc.deinit(gpa);

            var dft: @This() = undefined;

            dft.exch = exch_func;
            dft.corr = corr_func;
            dft.exco = exco_func;

            dft.pts = grid[0];
            dft.wgh = grid[1];

            dft.V_xc, dft.polarized, dft.a_x = .{ V_xc, polarized, exx };

            dft.E_xc = 0;

            return dft;
        }

        /// Frees all allocated grid arrays and libxc functional resources.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            if (self.exch) |*func| libxc.xc_func_end(func);
            if (self.corr) |*func| libxc.xc_func_end(func);
            if (self.exco) |*func| libxc.xc_func_end(func);

            self.pts.deinit(gpa);
            self.wgh.deinit(gpa);
            self.V_xc.deinit(gpa);
        }

        /// Evaluates the exchange-correlation potential matrix and energy contribution by grid integration.
        pub fn evaluate(self: *@This(), sys: MolecularSystem(T), P: Matrix(T), gpa: Allocator) !void {
            self.V_xc.zero();

            const size_factor: usize, const flags = .{ if (self.polarized) 2 else 1, self.getFamilyFlags() };

            const nd, const mgga = .{ flags.needs_deriv, flags.has_mgga };

            var basis = try BasisGrid(T).init(self.pts.shape[0], sys.nbf, nd, mgga, gpa);
            defer basis.deinit(gpa);

            var density = try DensityGrid(T).init(self.pts.shape[0], size_factor, nd, mgga, self.polarized, gpa);
            defer density.deinit(gpa);

            var potential = try PotentialGrid(T).init(self.pts.shape[0], size_factor, nd, mgga, self.polarized, gpa);
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

        /// Retrieves a formatted string of the active exchange and correlation functional names.
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

        /// Calculates the electronic density and its derivatives on the grid from the density matrix.
        fn computeDensity(self: @This(), P: Matrix(T), ctx: *DftEvaluationContext(T)) void {
            if (self.polarized) {
                self.computeDensityPolarized(P, ctx);
            } else {
                self.computeDensityUnpolarized(P, ctx);
            }
        }

        /// Calculates the spin-polarized electronic density, its gradient, and kinetic energy density on the grid.
        fn computeDensityPolarized(self: @This(), P: Matrix(T), ctx: *DftEvaluationContext(T)) void {
            const b = ctx.sys.nbf;

            for (0..self.pts.shape[0]) |g| {
                const rho_a = contractDensity(P, ctx.basisgrid, 0, g, b);
                const rho_b = contractDensity(P, ctx.basisgrid, b, g, b);

                ctx.density.rho_val.ptr(g, 0).* = rho_a;
                ctx.density.rho_val.ptr(g, 1).* = rho_b;

                if (ctx.needs_deriv) {
                    const del_rho_a = contractDensityGradient(P, ctx.basisgrid, 0, g, b);
                    const del_rho_b = contractDensityGradient(P, ctx.basisgrid, b, g, b);

                    ctx.density.del_rho.?.ptr(g, 0).* = del_rho_a[0];
                    ctx.density.del_rho.?.ptr(g, 1).* = del_rho_a[1];
                    ctx.density.del_rho.?.ptr(g, 2).* = del_rho_a[2];
                    ctx.density.del_rho.?.ptr(g, 3).* = del_rho_b[0];
                    ctx.density.del_rho.?.ptr(g, 4).* = del_rho_b[1];
                    ctx.density.del_rho.?.ptr(g, 5).* = del_rho_b[2];

                    const ga, const gb = .{ del_rho_a, del_rho_b };

                    const sig_aa = ga[0] * ga[0] + ga[1] * ga[1] + ga[2] * ga[2];
                    const sig_ab = ga[0] * gb[0] + ga[1] * gb[1] + ga[2] * gb[2];
                    const sig_bb = gb[0] * gb[0] + gb[1] * gb[1] + gb[2] * gb[2];

                    ctx.density.sig_val.?.ptr(g, 0).* = sig_aa;
                    ctx.density.sig_val.?.ptr(g, 1).* = sig_ab;
                    ctx.density.sig_val.?.ptr(g, 2).* = sig_bb;
                }

                if (ctx.has_mgga) {
                    const mgga_a = contractKineticEnergyDensity(P, ctx.basisgrid, 0, g, b);
                    const mgga_b = contractKineticEnergyDensity(P, ctx.basisgrid, b, g, b);

                    ctx.density.tau_val.?.ptr(g, 0).* = mgga_a.tau;
                    ctx.density.tau_val.?.ptr(g, 1).* = mgga_b.tau;

                    ctx.density.lap_val.?.ptr(g, 0).* = mgga_a.lapl;
                    ctx.density.lap_val.?.ptr(g, 1).* = mgga_b.lapl;
                }
            }
        }

        /// Calculates the spin-unpolarized electronic density, its gradient, and kinetic energy density on the grid.
        fn computeDensityUnpolarized(self: @This(), P: Matrix(T), ctx: *DftEvaluationContext(T)) void {
            for (0..self.pts.shape[0]) |g| {
                const rho_tot = contractDensity(P, ctx.basisgrid, 0, g, ctx.sys.nbf);

                ctx.density.rho_val.ptr(g, 0).* = rho_tot;

                if (ctx.needs_deriv) {
                    const del_rho = contractDensityGradient(P, ctx.basisgrid, 0, g, ctx.sys.nbf);

                    ctx.density.del_rho.?.ptr(g, 0).* = del_rho[0];
                    ctx.density.del_rho.?.ptr(g, 1).* = del_rho[1];
                    ctx.density.del_rho.?.ptr(g, 2).* = del_rho[2];

                    const sig_tot = del_rho[0] * del_rho[0] + del_rho[1] * del_rho[1] + del_rho[2] * del_rho[2];

                    ctx.density.sig_val.?.ptr(g, 0).* = sig_tot;
                }

                if (ctx.has_mgga) {
                    const mgga = contractKineticEnergyDensity(P, ctx.basisgrid, 0, g, ctx.sys.nbf);

                    ctx.density.tau_val.?.ptr(g, 0).* = mgga.tau;

                    ctx.density.lap_val.?.ptr(g, 0).* = mgga.lapl;
                }
            }
        }

        /// Computes the intermediate gradient vector X_mu = 2 * v_sigma * grad(rho) . grad(chi_mu) on the grid.
        fn computeIntermediateGradientVector(X: *Vector(T), del_rho: Matrix(T), basis: BasisGrid(T), v_sigma: T, g: usize, nbf: usize) void {
            const w_x = 2 * v_sigma * del_rho.at(g, 0);
            const w_y = 2 * v_sigma * del_rho.at(g, 1);
            const w_z = 2 * v_sigma * del_rho.at(g, 2);

            for (0..nbf) |mu| {
                const dx_mu = basis.dphi_dx.?.at(g, mu);
                const dy_mu = basis.dphi_dy.?.at(g, mu);
                const dz_mu = basis.dphi_dz.?.at(g, mu);

                X.data[mu] = w_x * dx_mu + w_y * dy_mu + w_z * dz_mu;
            }
        }

        /// Computes spin-resolved intermediate gradient vectors X_a and X_b on the numerical integration grid.
        fn computeIntermediateGradientVectorPolarized(X_a: *Vector(T), X_b: *Vector(T), del_rho: Matrix(T), basis: BasisGrid(T), sig_pot: Matrix(T), g: usize, nbf: usize) void {
            const g_a_x = del_rho.at(g, 0);
            const g_a_y = del_rho.at(g, 1);
            const g_a_z = del_rho.at(g, 2);

            const g_b_x = del_rho.at(g, 3);
            const g_b_y = del_rho.at(g, 4);
            const g_b_z = del_rho.at(g, 5);

            const v_sig_aa = sig_pot.at(g, 0);
            const v_sig_ab = sig_pot.at(g, 1);
            const v_sig_bb = sig_pot.at(g, 2);

            const w_a_x = 2 * v_sig_aa * g_a_x + v_sig_ab * g_b_x;
            const w_a_y = 2 * v_sig_aa * g_a_y + v_sig_ab * g_b_y;
            const w_a_z = 2 * v_sig_aa * g_a_z + v_sig_ab * g_b_z;

            const w_b_x = 2 * v_sig_bb * g_b_x + v_sig_ab * g_a_x;
            const w_b_y = 2 * v_sig_bb * g_b_y + v_sig_ab * g_a_y;
            const w_b_z = 2 * v_sig_bb * g_b_z + v_sig_ab * g_a_z;

            for (0..nbf) |mu| {
                const d_dx = basis.dphi_dx.?.at(g, mu);
                const d_dy = basis.dphi_dy.?.at(g, mu);
                const d_dz = basis.dphi_dz.?.at(g, mu);

                X_a.data[mu] = w_a_x * d_dx + w_a_y * d_dy + w_a_z * d_dz;
                X_b.data[mu] = w_b_x * d_dx + w_b_y * d_dy + w_b_z * d_dz;
            }
        }

        /// Evaluates the electron density rho at a grid point by contracting the density matrix with basis functions.
        fn contractDensity(P: Matrix(T), basis: BasisGrid(T), offset: usize, g: usize, nbf: usize) T {
            var rho: T = 0;

            for (0..nbf) |mu| for (0..nbf) |nu| {
                rho += P.at(mu + offset, nu + offset) * basis.phi.at(g, mu) * basis.phi.at(g, nu);
            };

            return rho;
        }

        /// Contracts atomic orbital gradients with the density matrix at a grid point to evaluate grad(rho).
        fn contractDensityGradient(P: Matrix(T), basis: BasisGrid(T), offset: usize, g: usize, nbf: usize) [3]T {
            var del_rho = [3]T{ 0, 0, 0 };

            for (0..nbf) |mu| {
                const dphi_dx_mu = basis.dphi_dx.?.at(g, mu);
                const dphi_dy_mu = basis.dphi_dy.?.at(g, mu);
                const dphi_dz_mu = basis.dphi_dz.?.at(g, mu);

                for (0..nbf) |nu| {
                    const term = P.at(mu + offset, nu + offset) * basis.phi.at(g, nu);

                    del_rho[0] += term * dphi_dx_mu;
                    del_rho[1] += term * dphi_dy_mu;
                    del_rho[2] += term * dphi_dz_mu;
                }
            }

            return .{ 2 * del_rho[0], 2 * del_rho[1], 2 * del_rho[2] };
        }

        /// Contracts atomic orbital gradients and Laplacians with the density matrix to evaluate tau and lapl(rho).
        fn contractKineticEnergyDensity(P: Matrix(T), basis: BasisGrid(T), offset: usize, g: usize, nbf: usize) struct { tau: T, lapl: T } {
            var tau: T, var lapl_sum: T = .{ 0, 0 };

            for (0..nbf) |mu| {
                const dphi_dx_mu = basis.dphi_dx.?.at(g, mu);
                const dphi_dy_mu = basis.dphi_dy.?.at(g, mu);
                const dphi_dz_mu = basis.dphi_dz.?.at(g, mu);

                const lapl_mu = basis.lapl.?.at(g, mu);

                for (0..nbf) |nu| {
                    const p_val = P.at(mu + offset, nu + offset);

                    const phi_nu = basis.phi.at(g, nu);

                    const dphi_dx_nu = basis.dphi_dx.?.at(g, nu);
                    const dphi_dy_nu = basis.dphi_dy.?.at(g, nu);
                    const dphi_dz_nu = basis.dphi_dz.?.at(g, nu);

                    const dphi_dot = dphi_dx_mu * dphi_dx_nu + dphi_dy_mu * dphi_dy_nu + dphi_dz_mu * dphi_dz_nu;

                    tau += p_val * dphi_dot;

                    lapl_sum += p_val * lapl_mu * phi_nu;
                }
            }

            return .{ .tau = 0.5 * tau, .lapl = 2 * (lapl_sum + tau) };
        }

        /// Evaluates atomic orbital basis functions and their spatial derivatives at all grid points.
        fn evaluateBasis(self: @This(), sys: MolecularSystem(T), basis: BasisGrid(T)) void {
            const needs_d1, const needs_d2 = .{ basis.dphi_dx != null, basis.lapl != null };

            for (0..self.pts.shape[0]) |g| {
                const x = self.pts.at(g, 0);
                const y = self.pts.at(g, 1);
                const z = self.pts.at(g, 2);

                if (needs_d2) {
                    const phi_row = basis.phi.rowSlice(g).ptr;

                    const dx_row = basis.dphi_dx.?.rowSlice(g).ptr;
                    const dy_row = basis.dphi_dy.?.rowSlice(g).ptr;
                    const dz_row = basis.dphi_dz.?.rowSlice(g).ptr;

                    const lapl_row = basis.lapl.?.rowSlice(g).ptr;

                    libint.libint_evaluate_basis_vgl(phi_row, dx_row, dy_row, dz_row, lapl_row, x, y, z, sys.ptr);
                }

                if (needs_d1 and !needs_d2) {
                    const phi_row = basis.phi.rowSlice(g).ptr;

                    const dx_row = basis.dphi_dx.?.rowSlice(g).ptr;
                    const dy_row = basis.dphi_dy.?.rowSlice(g).ptr;
                    const dz_row = basis.dphi_dz.?.rowSlice(g).ptr;

                    libint.libint_evaluate_basis_vg(phi_row, dx_row, dy_row, dz_row, x, y, z, sys.ptr);
                }

                if (!needs_d1 and !needs_d2) {
                    libint.libint_evaluate_basis_v(basis.phi.rowSlice(g).ptr, x, y, z, sys.ptr);
                }
            }
        }

        /// Analyzes the active functionals to determine GGA or meta-GGA classification and derivative requirements.
        fn getFamilyFlags(self: @This()) struct { has_sgga: bool, has_mgga: bool, needs_deriv: bool } {
            const exch_fam = if (self.exch) |func| func.info.*.family else null;
            const corr_fam = if (self.corr) |func| func.info.*.family else null;
            const exco_fam = if (self.exco) |func| func.info.*.family else null;

            const xc_gga, const xc_hyb_gga = .{ libxc.XC_FAMILY_GGA, libxc.XC_FAMILY_HYB_GGA };

            const is_exch_sgga = self.exch != null and (exch_fam == xc_gga or exch_fam == xc_hyb_gga);
            const is_corr_sgga = self.corr != null and (corr_fam == xc_gga or corr_fam == xc_hyb_gga);
            const is_exco_sgga = self.exco != null and (exco_fam == xc_gga or exco_fam == xc_hyb_gga);

            const xc_mgga, const xc_hyb_mgga = .{ libxc.XC_FAMILY_MGGA, libxc.XC_FAMILY_HYB_MGGA };

            const is_exch_mgga = self.exch != null and (exch_fam == xc_mgga or exch_fam == xc_hyb_mgga);
            const is_corr_mgga = self.corr != null and (corr_fam == xc_mgga or corr_fam == xc_hyb_mgga);
            const is_exco_mgga = self.exco != null and (exco_fam == xc_mgga or exco_fam == xc_hyb_mgga);

            const has_sgga = is_exch_sgga or is_corr_sgga or is_exco_sgga;
            const has_mgga = is_exch_mgga or is_corr_mgga or is_exco_mgga;

            return .{ .has_sgga = has_sgga, .has_mgga = has_mgga, .needs_deriv = has_sgga or has_mgga };
        }

        /// Allocates and initializes a libxc functional using its name and spin mode.
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

        /// Performs numerical integration of the exchange-correlation potential over grid points.
        fn integratePotential(self: *@This(), ctx: DftEvaluationContext(T), gpa: Allocator) !void {
            if (self.polarized) {
                self.E_xc = try self.integratePotentialPolarized(ctx, gpa);
            }

            if (!self.polarized) {
                self.E_xc = try self.integratePotentialUnpolarized(ctx, gpa);
            }
        }

        /// Integrates the spin-polarized exchange-correlation potential to yield the energy and matrix.
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

            for (0..self.pts.shape[0]) |g| {
                const w = self.wgh.at(g);

                energy_exc += (ctx.density.rho_val.at(g, 0) + ctx.density.rho_val.at(g, 1)) * ctx.potgd.exc.at(g) * w;

                const v_a = ctx.potgd.rho_pot.at(g, 0);
                const v_b = ctx.potgd.rho_pot.at(g, 1);

                if (ctx.needs_deriv) {
                    computeIntermediateGradientVectorPolarized(&X_a.?, &X_b.?, ctx.density.del_rho.?, ctx.basisgrid, ctx.potgd.sig_pot.?, g, ctx.sys.nbf);

                    const v_tau_a = if (ctx.has_mgga) ctx.potgd.tau_pot.?.at(g, 0) else 0.0;
                    const v_tau_b = if (ctx.has_mgga) ctx.potgd.tau_pot.?.at(g, 1) else 0.0;
                    const v_lap_a = if (ctx.has_mgga) ctx.potgd.lap_pot.?.at(g, 0) else 0.0;
                    const v_lap_b = if (ctx.has_mgga) ctx.potgd.lap_pot.?.at(g, 1) else 0.0;

                    for (0..ctx.sys.nbf) |mu| for (0..ctx.sys.nbf) |nu| {
                        const phi_mu = ctx.basisgrid.phi.at(g, mu);
                        const phi_nu = ctx.basisgrid.phi.at(g, nu);

                        const b = ctx.sys.nbf;

                        var val_a = v_a * phi_mu * phi_nu + X_a.?.data[mu] * phi_nu + phi_mu * X_a.?.data[nu];
                        var val_b = v_b * phi_mu * phi_nu + X_b.?.data[mu] * phi_nu + phi_mu * X_b.?.data[nu];

                        if (ctx.has_mgga) {
                            const dx_mu_nu = ctx.basisgrid.dphi_dx.?.at(g, mu) * ctx.basisgrid.dphi_dx.?.at(g, nu);
                            const dy_mu_nu = ctx.basisgrid.dphi_dy.?.at(g, mu) * ctx.basisgrid.dphi_dy.?.at(g, nu);
                            const dz_mu_nu = ctx.basisgrid.dphi_dz.?.at(g, mu) * ctx.basisgrid.dphi_dz.?.at(g, nu);

                            const dphi_dot = dx_mu_nu + dy_mu_nu + dz_mu_nu;

                            val_a += 0.5 * v_tau_a * dphi_dot;
                            val_b += 0.5 * v_tau_b * dphi_dot;

                            const lapl_mu = ctx.basisgrid.lapl.?.at(g, mu);
                            const lapl_nu = ctx.basisgrid.lapl.?.at(g, nu);

                            const lapl_term = lapl_mu * phi_nu + 2 * dphi_dot + phi_mu * lapl_nu;

                            val_a += v_lap_a * lapl_term;
                            val_b += v_lap_b * lapl_term;
                        }

                        self.V_xc.ptr(mu + 0, nu + 0).* += val_a * w;
                        self.V_xc.ptr(mu + b, nu + b).* += val_b * w;
                    };
                }

                if (!ctx.needs_deriv) {
                    for (0..ctx.sys.nbf) |mu| for (0..ctx.sys.nbf) |nu| {
                        const phi_mu = ctx.basisgrid.phi.at(g, mu);
                        const phi_nu = ctx.basisgrid.phi.at(g, nu);

                        const b = ctx.sys.nbf;

                        self.V_xc.ptr(mu + 0, nu + 0).* += v_a * phi_mu * phi_nu * w;
                        self.V_xc.ptr(mu + b, nu + b).* += v_b * phi_mu * phi_nu * w;
                    };
                }
            }
            return energy_exc;
        }

        /// Integrates the spin-unpolarized exchange-correlation potential to yield the energy and matrix.
        fn integratePotentialUnpolarized(self: *@This(), ctx: DftEvaluationContext(T), gpa: Allocator) !T {
            var energy_exc: T, var X_a: ?Vector(T) = .{ 0, null };

            if (ctx.needs_deriv) {
                X_a = try Vector(T).init(ctx.sys.nbf, gpa);
            }

            defer {
                if (X_a) |*x| x.deinit(gpa);
            }

            for (0..self.pts.shape[0]) |g| {
                const w = self.wgh.at(g);

                energy_exc += ctx.density.rho_val.at(g, 0) * ctx.potgd.exc.at(g) * w;

                if (ctx.needs_deriv) {
                    const v_sigma = ctx.potgd.sig_pot.?.at(g, 0);

                    computeIntermediateGradientVector(&X_a.?, ctx.density.del_rho.?, ctx.basisgrid, v_sigma, g, ctx.sys.nbf);

                    const v_tau = if (ctx.has_mgga) ctx.potgd.tau_pot.?.at(g, 0) else 0.0;
                    const v_lap = if (ctx.has_mgga) ctx.potgd.lap_pot.?.at(g, 0) else 0.0;

                    for (0..ctx.sys.nbf) |mu| for (0..ctx.sys.nbf) |nu| {
                        const v_rho = ctx.potgd.rho_pot.at(g, 0);

                        const phi_mu = ctx.basisgrid.phi.at(g, mu);
                        const phi_nu = ctx.basisgrid.phi.at(g, nu);

                        var val = v_rho * phi_mu * phi_nu + X_a.?.data[mu] * phi_nu + phi_mu * X_a.?.data[nu];

                        if (ctx.has_mgga) {
                            const dx_mu_nu = ctx.basisgrid.dphi_dx.?.at(g, mu) * ctx.basisgrid.dphi_dx.?.at(g, nu);
                            const dy_mu_nu = ctx.basisgrid.dphi_dy.?.at(g, mu) * ctx.basisgrid.dphi_dy.?.at(g, nu);
                            const dz_mu_nu = ctx.basisgrid.dphi_dz.?.at(g, mu) * ctx.basisgrid.dphi_dz.?.at(g, nu);

                            const dphi_dot = dx_mu_nu + dy_mu_nu + dz_mu_nu;

                            val += 0.5 * v_tau * dphi_dot;

                            const lapl_mu = ctx.basisgrid.lapl.?.at(g, mu);
                            const lapl_nu = ctx.basisgrid.lapl.?.at(g, nu);

                            const lapl_term = lapl_mu * phi_nu + 2 * dphi_dot + phi_mu * lapl_nu;

                            val += v_lap * lapl_term;
                        }

                        self.V_xc.ptr(mu, nu).* += val * w;
                    };
                }

                if (!ctx.needs_deriv) for (0..ctx.sys.nbf) |mu| for (0..ctx.sys.nbf) |nu| {
                    const phi_mu = ctx.basisgrid.phi.at(g, mu);
                    const phi_nu = ctx.basisgrid.phi.at(g, nu);

                    self.V_xc.ptr(mu, nu).* += ctx.potgd.rho_pot.at(g, 0) * phi_mu * phi_nu * w;
                };
            }

            return energy_exc;
        }
    };
}

/// Holds temporary arrays and grid metrics during DFT exchange-correlation potential evaluation.
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
