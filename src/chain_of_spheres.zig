//! Implements the Chain-of-Spheres (COSX) semi-numerical exchange approximation for Hartree–Fock and DFT.

const std = @import("std");

const Allocator = std.mem.Allocator;

const BasisGrid = @import("molecular_grid.zig").BasisGrid;
const Becke = @import("molecular_grid.zig").Becke;
const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Vector = @import("tensor.zig").Vector;

const gemm = @import("linear_algebra.zig").gemm;

/// Angular and radial grid specifications for Chain-of-Spheres (COSX) numerical integration.
pub const GridOptions = struct {
    angular: usize = 50,
    radial: usize = 25,
};

/// Options for the Chain-of-Spheres (COSX) semi-numerical exchange approximation.
pub const Options = struct {
    grid: GridOptions = .{},
};

/// Numerical integration grid and potential kernel for evaluating the COSX exchange matrix.
pub fn ChainOfSpheres(comptime T: type) type {
    return struct {
        kernel: Matrix(T),

        phi: Matrix(T),
        wgh: Vector(T),

        /// Initializes the COSX integration grid, evaluates atomic basis functions, and constructs the Coulomb kernel.
        pub fn init(sys: MolecularSystem(T), opt: Options, gpa: Allocator) !@This() {
            var grid = try (Becke(T){ .n_rad = opt.grid.radial, .n_leb = opt.grid.angular }).get(sys, gpa);
            defer grid[0].deinit(gpa);

            errdefer grid[1].deinit(gpa);

            const n_grid = grid[0].shape[0];

            var basis_grid = try BasisGrid(T).init(n_grid, sys.nbf, false, false, gpa);
            defer basis_grid.deinit(gpa);

            basis_grid.evaluate(sys, grid[0]);

            var phi = try basis_grid.phi.clone(gpa);
            errdefer phi.deinit(gpa);

            var kernel = try Matrix(T).init(n_grid, n_grid, gpa);
            errdefer kernel.deinit(gpa);

            for (0..n_grid) |a| {
                const xa = grid[0].at(a, 0);
                const ya = grid[0].at(a, 1);
                const za = grid[0].at(a, 2);

                const wa = grid[1].at(a);

                const r_self = 2.0 * std.math.pi * std.math.pow(T, 0.75 * wa / std.math.pi, 2.0 / 3.0);

                for (0..a) |b| {
                    const xb = grid[0].at(b, 0);
                    const yb = grid[0].at(b, 1);
                    const zb = grid[0].at(b, 2);

                    const dx = xa - xb;
                    const dy = ya - yb;
                    const dz = za - zb;

                    const r = @sqrt(dx * dx + dy * dy + dz * dz);

                    kernel.ptr(a, b).*, kernel.ptr(b, a).* = .{ grid[1].at(b) / r, wa / r };
                }

                kernel.ptr(a, a).* = r_self;
            }

            return .{ .kernel = kernel, .phi = phi, .wgh = grid[1] };
        }

        /// Deallocates grid weights, basis grid values, and the Coulomb kernel matrix.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.kernel.deinit(gpa);

            self.phi.deinit(gpa);
            self.wgh.deinit(gpa);
        }

        /// Evaluates the semi-numerical exchange matrix using occupied MO coefficients and accumulates into F.
        pub fn contractExchange(self: *@This(), F: *Matrix(T), C: Matrix(T), nocc: usize, exch_factor: T, gpa: Allocator) !void {
            const n_grid, const nbf = .{ self.phi.shape[0], self.phi.shape[1] };

            const scale: T = if (C.shape[0] != nbf) -exch_factor else -2.0 * exch_factor;

            var K_tot = try Matrix(T).initZero(nbf, nbf, gpa);
            defer K_tot.deinit(gpa);

            var Psi = try Matrix(T).init(n_grid, nocc, gpa);
            defer Psi.deinit(gpa);

            const C_spatial = if (C.shape[0] != nbf) C.sliceRows(0, nbf) else C;

            var C_occ = try Matrix(T).init(nbf, nocc, gpa);
            defer C_occ.deinit(gpa);

            for (0..nbf) |mu| for (0..nocc) |i| {
                C_occ.ptr(mu, i).* = C_spatial.at(mu, i);
            };

            gemm(T, &Psi, self.phi, C_occ, 1, 0, false, false);

            var rho_i = try Matrix(T).init(n_grid, nbf, gpa);
            defer rho_i.deinit(gpa);

            var A_i = try Matrix(T).init(n_grid, nbf, gpa);
            defer A_i.deinit(gpa);

            var Q_i = try Matrix(T).init(nbf, n_grid, gpa);
            defer Q_i.deinit(gpa);

            var K_i = try Matrix(T).init(nbf, nbf, gpa);
            defer K_i.deinit(gpa);

            for (0..nocc) |i| {
                for (0..n_grid) |g| {
                    const psi_val = Psi.at(g, i);

                    const w_psi = self.wgh.at(g) * psi_val;

                    for (0..nbf) |nu| {
                        const phi_val = self.phi.at(g, nu);

                        rho_i.ptr(g, nu).*, Q_i.ptr(nu, g).* = .{ phi_val * psi_val, phi_val * w_psi };
                    }
                }

                gemm(T, &A_i, self.kernel, rho_i, 1, 0, false, false);

                gemm(T, &K_i, Q_i, A_i, 1, 0, false, false);

                for (0..nbf) |mu| for (0..nbf) |nu| {
                    K_tot.ptr(mu, nu).* += K_i.at(mu, nu);
                };
            }

            for (0..nbf) |mu| for (0..nbf) |nu| {
                const k_sym = 0.5 * (K_tot.at(mu, nu) + K_tot.at(nu, mu));

                F.ptr(mu, nu).* += scale * k_sym;

                if (C.shape[0] != nbf) {
                    F.ptr(mu + nbf, nu + nbf).* += scale * k_sym;
                }
            };
        }
    };
}
