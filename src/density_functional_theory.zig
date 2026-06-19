const std = @import("std");
const libint = @import("libint");

const libxc = @cImport(@cInclude("xc.h"));

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Vector = @import("tensor.zig").Vector;

const getLebedevGrid = @import("lebedev_quadrature_nodes.zig").getLebedevGrid;

const A2BOHR = @import("constant.zig").A2BOHR;

// DFT STRUCTURE =======================================================================================================

pub fn DftPotential(comptime T: type) type {
    return struct {
        exch: libxc.xc_func_type,
        corr: libxc.xc_func_type,

        pts: Matrix(T),
        wgh: Vector(T),

        Vxc: Matrix(T),
        Exc: T = 0.000,

        pub fn init(sys: MolecularSystem(T), exch: []const u8, corr: []const u8, n_rad: usize, n_leb: usize, gpa: Allocator) !@This() {
            const grid = try (Becke(T){ .n_rad = n_rad, .n_leb = n_leb }).get(sys, gpa);

            const exch_nt = try gpa.dupeSentinel(u8, exch, 0);
            defer gpa.free(exch_nt);

            const corr_nt = try gpa.dupeSentinel(u8, corr, 0);
            defer gpa.free(corr_nt);

            const id_exch = libxc.xc_functional_get_number(exch_nt.ptr);
            const id_corr = libxc.xc_functional_get_number(corr_nt.ptr);

            if (id_exch < 0) return error.ExchFuncNotFound;
            if (id_corr < 0) return error.CorrFuncNotFound;

            var exch_func: libxc.xc_func_type = undefined;
            var corr_func: libxc.xc_func_type = undefined;

            if (libxc.xc_func_init(&exch_func, id_exch, libxc.XC_UNPOLARIZED) != 0) return error.LibxcInitFailed;
            errdefer libxc.xc_func_end(&exch_func);

            if (libxc.xc_func_init(&corr_func, id_corr, libxc.XC_UNPOLARIZED) != 0) return error.LibxcInitFailed;
            errdefer libxc.xc_func_end(&corr_func);

            var Vxc = try Matrix(T).init(sys.nbf, sys.nbf, gpa);
            errdefer Vxc.deinit(gpa);

            return .{ .exch = exch_func, .corr = corr_func, .pts = grid[0], .wgh = grid[1], .Vxc = Vxc };
        }

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            libxc.xc_func_end(&self.exch);
            libxc.xc_func_end(&self.corr);

            self.pts.deinit(gpa);
            self.wgh.deinit(gpa);
            self.Vxc.deinit(gpa);
        }

        pub fn evaluate(self: *@This(), sys: MolecularSystem(T), P: Matrix(T), gpa: Allocator) !void {
            self.Vxc.zero();

            const rho = try gpa.alloc(T, self.pts.shape[0]);
            defer gpa.free(rho);

            @memset(rho, 0.0);

            const exc = try gpa.alloc(T, self.pts.shape[0]);
            defer gpa.free(exc);

            @memset(exc, 0.0);

            const vrho = try gpa.alloc(T, self.pts.shape[0]);
            defer gpa.free(vrho);

            @memset(vrho, 0.0);

            var phi = try Matrix(T).init(self.pts.shape[0], sys.nbf, gpa);
            defer phi.deinit(gpa);

            for (0..self.pts.shape[0]) |i| {
                const x = self.pts.at(i, 0);
                const y = self.pts.at(i, 1);
                const z = self.pts.at(i, 2);

                libint.evaluate_basis(phi.rowSlice(i).ptr, x, y, z, sys.ptr);

                for (0..sys.nbf) |mu| for (0..sys.nbf) |nu| {
                    rho[i] += P.at(mu, nu) * phi.at(i, mu) * phi.at(i, nu);
                };
            }

            inline for (.{ self.exch, self.corr }) |name| {
                evaluateXCFunctional(name, rho, exc, vrho);
            }

            self.Exc = 0.0;

            for (0..self.pts.shape[0]) |i| {
                const w = self.wgh.at(i);

                self.Exc += rho[i] * exc[i] * w;

                for (0..sys.nbf) |mu| for (0..sys.nbf) |nu| {
                    self.Vxc.ptr(mu, nu).* += vrho[i] * phi.at(i, mu) * phi.at(i, nu) * w;
                };
            }
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
    const r_f64: T = switch (Z) {
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

    return @as(T, @floatCast(r_f64 * A2BOHR));
}

fn evaluateXCFunctional(func: libxc.xc_func_type, rho: []const f64, exc: []f64, vrh: []f64) void {
    const CHUNK_SIZE = 512;

    var exc_tmp: [CHUNK_SIZE]f64 = undefined;
    var vrh_tmp: [CHUNK_SIZE]f64 = undefined;

    var i: usize = 0;

    while (i < rho.len) {
        const batch_len = @min(CHUNK_SIZE, rho.len - i);

        const exc_batch = exc_tmp[0..batch_len];
        const vrh_batch = vrh_tmp[0..batch_len];

        @memset(exc_batch, 0.0);
        @memset(vrh_batch, 0.0);

        libxc.xc_lda_exc_vxc(&func, batch_len, rho[i .. i + batch_len].ptr, exc_batch.ptr, vrh_batch.ptr);

        for (0..batch_len) |j| {
            exc[i + j] += exc_batch[j];
            vrh[i + j] += vrh_batch[j];
        }

        i += batch_len;
    }
}
