//! Represents a molecular grid, managing radial and angular quadrature integration grids for density functional theory.

const std = @import("std");
const libint = @import("libint");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Vector = @import("tensor.zig").Vector;

const getLebedevGrid = @import("lebedev_quadrature_nodes.zig").getLebedevGrid;

const A2BOHR = @import("constant.zig").A2BOHR;

/// Returns a generic type representing a molecular grid for basis function evaluation at grid points.
pub fn BasisGrid(comptime T: type) type {
    return struct {
        phi: Matrix(T),

        dphi_dx: ?Matrix(T) = null,
        dphi_dy: ?Matrix(T) = null,
        dphi_dz: ?Matrix(T) = null,

        lapl: ?Matrix(T) = null,

        /// Initializes a BasisGrid structure, allocating memory for basis functions and optional derivatives at grid points.
        pub fn init(n_grid: usize, nbf: usize, needs_deriv: bool, has_mgga: bool, gpa: Allocator) !@This() {
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

        /// Frees allocated memory for basis functions and derivative values at grid points.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.phi.deinit(gpa);

            if (self.dphi_dx) |*d| d.deinit(gpa);
            if (self.dphi_dy) |*d| d.deinit(gpa);
            if (self.dphi_dz) |*d| d.deinit(gpa);

            if (self.lapl) |*d| d.deinit(gpa);
        }
    };
}

/// Returns a generic type implementing Becke multicenter integration partition schemes for molecular grids.
pub fn Becke(comptime T: type) type {
    return struct {
        n_rad: usize,
        n_leb: usize,

        /// Computes the Becke partition coordinate transform step function for atomic cell smoothing.
        fn beckeStep(mu: T) T {
            var p = 1.5 * mu - 0.5 * mu * mu * mu;

            p = 1.5 * p - 0.5 * p * p * p;
            p = 1.5 * p - 0.5 * p * p * p;

            return 0.5 * (1.0 - p);
        }

        /// Generates Becke partitioned grid points and weights for multicenter molecular integration.
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

        /// Computes Chebyshev radial quadrature nodes and weights mapped to a physical atomic radius.
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

/// Returns a generic type representing density values, gradients, laplacians, and kinetic energy density on the grid.
pub fn DensityGrid(comptime T: type) type {
    return struct {
        rho_val: Matrix(T),

        del_rho: ?Matrix(T) = null,
        sig_val: ?Matrix(T) = null,
        tau_val: ?Matrix(T) = null,
        lap_val: ?Matrix(T) = null,

        /// Initializes a DensityGrid structure, allocating memory for density, gradient, and laplacian grids.
        pub fn init(n_g: usize, size_factor: usize, needs_deriv: bool, has_mgga: bool, pol: bool, gpa: Allocator) !@This() {
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

        /// Frees allocated memory for the grid density and derivative values.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.rho_val.deinit(gpa);

            if (self.del_rho) |*g| g.deinit(gpa);
            if (self.sig_val) |*s| s.deinit(gpa);
            if (self.tau_val) |*t| t.deinit(gpa);
            if (self.lap_val) |*l| l.deinit(gpa);
        }
    };
}

/// Returns a generic type representing DFT exchange-correlation potentials and their gradients on the grid.
pub fn PotentialGrid(comptime T: type) type {
    return struct {
        exc: Vector(T),

        rho_pot: Matrix(T),

        sig_pot: ?Matrix(T) = null,
        tau_pot: ?Matrix(T) = null,
        lap_pot: ?Matrix(T) = null,

        /// Initializes a PotentialGrid, allocating memory for exchange-correlation potential terms on grid points.
        pub fn init(n_g: usize, size_factor: usize, needs_deriv: bool, has_mgga: bool, pol: bool, gpa: Allocator) !@This() {
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

        /// Frees allocated memory for the grid potential terms.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.exc.deinit(gpa);

            self.rho_pot.deinit(gpa);

            if (self.sig_pot) |*s| s.deinit(gpa);
            if (self.tau_pot) |*t| t.deinit(gpa);
            if (self.lap_pot) |*l| l.deinit(gpa);
        }
    };
}

/// Retrieves the Bragg-Slater radius in Bohr units for a given atomic number Z.
pub fn getBraggSlaterRadius(comptime T: type, Z: usize) T {
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
