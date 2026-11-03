//! Computes Møller-Plesset correlation energy corrections to a Hartree-Fock reference wavefunction.

const std = @import("std");

const cblas = @import("cimport.zig").cblas;

const Allocator = std.mem.Allocator;

const FrequencyOptions = @import("frequency_analysis.zig").Options;
const HartreeFockOptions = @import("hartree_fock.zig").Options;
const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const ScalarDual = @import("dual.zig").ScalarDual;
const Tensor = @import("tensor.zig").Tensor;
const Value = @import("value.zig").Value;
const Vector = @import("tensor.zig").Vector;

const ao2mo_pppp = @import("integral_transform.zig").ao2mo_pppp;
const ao2so_coef = @import("integral_transform.zig").ao2so_coef;
const ao2so_pppp = @import("integral_transform.zig").ao2so_pppp;
const generateDets = @import("configuration_interaction.zig").generateDets;
const printf = @import("read_write.zig").printf;
const slater = @import("configuration_interaction.zig").slater;

const AU2CM = @import("constant.zig").AU2CM;

/// Parameters governing the Møller-Plesset perturbation theory calculation and its derivatives.
pub const Options = struct {
    frequency: ?FrequencyOptions = null,
    gradient: ?GradientOptions = null,
    hartree_fock: HartreeFockOptions,
    hessian: ?HessianOptions = null,
    nthreads: u32 = 1,
    optimize: ?OptimizeOptions = null,
    order: u32 = 2,
    write: Write = .{},
};

/// Options for computing Moller-Plesset energy gradients analytically or numerically.
pub const GradientOptions = union(enum) {
    analytic: GradientAnalyticOptions,
    numeric: GradientNumericOptions,
};

/// Tagged union specifying Moller-Plesset nuclear Hessian calculation methods.
pub const HessianOptions = union(enum) {
    numeric: HessianNumericOptions,
};

/// Tagged union specifying geometry optimization algorithms.
pub const OptimizeOptions = union(enum) {
    bfgs: BfgsOptions,
    steepest_descent: SteepestDescentOptions,
};

/// Configuration options for BFGS quasi-Newton molecular geometry optimization.
const BfgsOptions = struct {
    gradient: GradientOptions = .analytic,
    iterations: u32 = 100,
    step: f64 = 1,
    threshold: f64 = 1e-4,
};

/// Analytical nuclear gradient evaluation settings.
const GradientAnalyticOptions = struct {};

/// Finite difference displacement settings for numerical gradient evaluation.
const GradientNumericOptions = struct {
    step: f64 = 1e-5,
};

/// Finite difference displacement step for numerical Hessian evaluation.
const HessianNumericOptions = struct {
    step: f64 = 1e-5,
};

/// Configuration options for steepest descent molecular geometry optimization.
const SteepestDescentOptions = struct {
    gradient: GradientOptions = .analytic,
    iterations: u32 = 100,
    step: f64 = 1e-1,
    threshold: f64 = 1e-4,
};

/// Destination paths for outputting Moller-Plesset geometries, gradients, and Hessians.
const Write = struct {
    geometry: ?[]const u8 = null,
    gradient: ?[]const u8 = null,
    hessian: ?[]const u8 = null,
};

/// Holds Moller-Plesset calculation outputs: HF reference result, perturbation energies, and derivatives.
pub fn Result(comptime T: type) type {
    return struct {
        hartree_fock: HartreeFockResult(T),

        energy: []T,

        grad: []Matrix(T) = &.{},
        hess: []Matrix(T) = &.{},

        /// Frees all allocated resources associated with the Møller-Plesset result.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.hartree_fock.deinit(gpa);

            gpa.free(self.energy);

            for (0..self.grad.len) |i| {
                self.grad[i].deinit(gpa);
            }

            if (self.grad.len > 0) gpa.free(self.grad);

            for (0..self.hess.len) |i| {
                self.hess[i].deinit(gpa);
            }

            if (self.hess.len > 0) gpa.free(self.hess);
        }
    };
}

/// Computes the Møller-Plesset perturbation series for electronic correlation energy corrections.
pub fn solve(comptime T: type, io: std.Io, opt: Options, sys: *MolecularSystem(T), ref: HartreeFockResult(T), log: bool, gpa: Allocator) !Result(T) {
    _ = sys;

    var hfres = ref;
    errdefer hfres.deinit(gpa);

    cblas.openblas_set_num_threads(@intCast(opt.nthreads));

    var energy = try gpa.alloc(T, 1);
    errdefer gpa.free(energy);

    var grad: []Matrix(T) = &.{};

    errdefer {
        for (0..grad.len) |k| grad[k].deinit(gpa);

        if (grad.len > 0) gpa.free(grad);
    }

    var grads: ?[]Matrix(T) = null;

    defer if (grads) |g| {
        for (0..(opt.order + 1)) |i| g[i].deinit(gpa);

        gpa.free(g);
    };

    if (opt.gradient) |g| switch (g) {
        .analytic => {
            var g0 = try hfres.grad[0].clone(gpa);
            errdefer g0.deinit(gpa);

            grad = try gpa.alloc(Matrix(T), 1);
            grad[0] = g0;

            grads = try gradient(T, opt.order, hfres, gpa);
        },
        .numeric => {},
    };

    energy[0] = hfres.energy[0];

    const generalized = opt.hartree_fock.generalized;

    const nocc = if (generalized) hfres.ints.sys.nel else hfres.ints.sys.nel / 2;

    const energies = try mp(T, opt.order, hfres.ints.g.?, hfres.C_mo, hfres.e_mo, nocc, generalized, gpa);
    defer gpa.free(energies);

    for (2..opt.order + 1) |i| {
        energy[0] += energies[i];

        if (log) {
            try printf(io, "\nMP{d:02} TOTAL ENERGY: {d:.14} Eh", .{ i, energy[0] });

            if (grads != null) {
                try printf(io, "\n", .{});
            }
        }

        if (grads) |g| {
            for (0..grad[0].nrow()) |j| for (0..grad[0].ncol()) |k| {
                grad[0].ptr(j, k).* += g[i].at(j, k);
            };

            if (log) {
                try printf(io, "\nMP{d:02} ANALYTICAL NUCLEAR ENERGY GRADIENT (Eh/a0)\n", .{i});

                for (0..grad[0].nrow()) |j| for (0..grad[0].ncol()) |k| {
                    try printf(io, "{d:20.14}{s}", .{ grad[0].at(j, k), if (k == 2) "\n" else " " });
                };
            }
        }

        if (log and grads == null and i == opt.order) {
            try printf(io, "\n", .{});
        }
    }

    return Result(T){ .hartree_fock = hfres, .energy = energy, .grad = grad, .hess = &.{} };
}

/// Evaluates Hamiltonian matrix elements over all Slater determinant pairs using Slater-Condon rules.
fn buildCiMatrix(comptime T: type, dets: []const []const usize, H_MS: Matrix(T), g_MS: Tensor(T, 4), gpa: Allocator) !Matrix(T) {
    var H_CI = try Matrix(T).init(dets.len, dets.len, gpa);
    errdefer H_CI.deinit(gpa);

    for (0..dets.len) |i| for (i..dets.len) |j| {
        const val = slater(T, dets[i], dets[j], H_MS, g_MS);

        H_CI.ptr(i, j).* = val;
        H_CI.ptr(j, i).* = val;
    };

    return H_CI;
}

/// Computes the nuclear gradient of the Møller-Plesset energy corrections using dual number differentiation.
fn gradient(comptime T: type, order: usize, hfres: HartreeFockResult(T), gpa: Allocator) ![]Matrix(T) {
    const nbf, const generalized = .{ hfres.C_mo.shape[0], hfres.ints.sys.nbf != hfres.C_mo.shape[0] };

    const nocc = if (generalized) hfres.ints.sys.nel else hfres.ints.sys.nel / 2;

    const grads = try gpa.alloc(Matrix(T), order + 1);
    errdefer gpa.free(grads);

    for (0..order + 1) |k| {
        grads[k] = Matrix(T).initZero(hfres.ints.sys.atoms.len, 3, gpa) catch |err| {
            for (0..k) |i| {
                grads[i].deinit(gpa);
            }

            return err;
        };
    }

    for (0..3 * hfres.ints.sys.atoms.len) |i| {
        var C_d = try Matrix(ScalarDual(T)).init(nbf, nbf, gpa);
        defer C_d.deinit(gpa);

        for (0..nbf) |mu| for (0..nbf) |p| {
            C_d.ptr(mu, p).* = ScalarDual(T).init(hfres.C_mo.at(mu, p), hfres.dC_mo.?.at(.{ i, mu, p }));
        };

        var e_d = try Vector(ScalarDual(T)).init(nbf, gpa);
        defer e_d.deinit(gpa);

        for (0..nbf) |j| {
            e_d.ptr(j).* = ScalarDual(T).init(hfres.e_mo.at(j), hfres.de_mo.?.at(i, j));
        }

        var g_d = try Tensor(ScalarDual(T), 4).init(hfres.ints.g.?.shape, gpa);
        defer g_d.deinit(gpa);

        for (0..nbf) |mu| for (0..nbf) |lambda| for (0..nbf) |nu| for (0..nbf) |sigma| {
            const dg_x = hfres.ints.dg.?.at(.{ i, mu, lambda, nu, sigma });

            g_d.ptr(.{ mu, lambda, nu, sigma }).* = ScalarDual(T).init(hfres.ints.g.?.at(.{ mu, lambda, nu, sigma }), dg_x);
        };

        const corr_energies_d = try mp(ScalarDual(T), order, g_d, C_d, e_d, nocc, generalized, gpa);
        defer gpa.free(corr_energies_d);

        for (2..order + 1) |k| {
            grads[k].ptr(i / 3, i % 3).* = corr_energies_d[k].der;
        }
    }

    return grads;
}

/// Computes Møller-Plesset energy corrections up to the specified perturbation order.
fn mp(comptime T: type, order: usize, g: Tensor(T, 4), C: Matrix(T), e: Vector(T), nocc: usize, generalized: bool, gpa: Allocator) ![]T {
    const nsp, const nel = if (generalized) .{ C.shape[0], nocc } else .{ 2 * C.shape[0], 2 * nocc };

    var C_SO = if (generalized) C else try Matrix(T).initZero(nsp, nsp, gpa);
    defer if (!generalized) C_SO.deinit(gpa);

    if (!generalized) {
        ao2so_coef(T, &C_SO, C);
    }

    var g_SO = if (generalized) g else try Tensor(T, 4).initZero(.{ nsp, nsp, nsp, nsp }, gpa);
    defer if (!generalized) g_SO.deinit(gpa);

    if (!generalized) {
        ao2so_pppp(T, &g_SO, g);
    }

    var g_MS = try Tensor(T, 4).init(.{ nsp, nsp, nsp, nsp }, gpa);
    defer g_MS.deinit(gpa);

    try ao2mo_pppp(T, &g_MS, g_SO, C_SO, gpa);

    var H_MS = try Matrix(T).initZero(nsp, nsp, gpa);
    defer H_MS.deinit(gpa);

    var e_SO = try Vector(T).init(nsp, gpa);
    defer e_SO.deinit(gpa);

    if (generalized) for (0..nsp) |p| {
        e_SO.ptr(p).* = e.at(p);
    };

    if (!generalized) for (0..C.shape[0]) |i| {
        e_SO.ptr(2 * i + 0).* = e.at(i);
        e_SO.ptr(2 * i + 1).* = e.at(i);
    };

    for (0..nsp) |p| for (0..nsp) |q| {
        var sum = Value(T).fromFloat(0);

        for (0..nel) |i| {
            const term1 = Value(T).init(g_MS.at(.{ p, i, q, i }));
            const term2 = Value(T).init(g_MS.at(.{ p, i, i, q }));

            sum = sum.add(term1.sub(term2));
        }

        H_MS.ptr(p, q).* = (if (p == q) Value(T).init(e_SO.at(p)) else Value(T).fromFloat(0)).sub(sum).val;
    };

    var excitations = try gpa.alloc(u32, @min(order, @min(nel, nsp - nel)));
    defer gpa.free(excitations);

    for (0..excitations.len) |i| {
        excitations[i] = @intCast(i + 1);
    }

    var dets = try generateDets(nel, nsp, excitations, gpa);

    defer {
        for (0..dets.items.len) |i| gpa.free(dets.items[i]);

        dets.deinit(gpa);
    }

    var H_CI = try buildCiMatrix(T, dets.items, H_MS, g_MS, gpa);
    defer H_CI.deinit(gpa);

    var E0 = try gpa.alloc(Value(T), dets.items.len);
    defer gpa.free(E0);

    for (0..dets.items.len) |i| {
        var sum = Value(T).fromFloat(0);

        for (dets.items[i]) |p| {
            sum = sum.add(Value(T).init(e_SO.at(p)));
        }

        E0[i] = sum;
    }

    return try solveRayleighSchrodinger(T, H_CI, E0, order, gpa);
}

/// Solves Rayleigh-Schrödinger perturbation equations iteratively up to the requested perturbation order.
fn solveRayleighSchrodinger(comptime T: type, H_CI: Matrix(T), E0: []const Value(T), order: usize, gpa: Allocator) ![]T {
    const ndets = E0.len;

    var C_coeff = try gpa.alloc(Value(T), (order + 1) * ndets);
    defer gpa.free(C_coeff);

    @memset(C_coeff, Value(T).fromFloat(0));

    var E = try gpa.alloc(Value(T), order + 1);
    defer gpa.free(E);

    E[0], C_coeff[0] = .{ E0[0], Value(T).fromFloat(1) };

    for (1..order + 1) |k| {
        if (k == 1) {
            E[1] = Value(T).init(H_CI.at(0, 0)).sub(E0[0]);
        }

        if (k > 1) {
            var sum = Value(T).fromFloat(0);

            for (1..ndets) |j| {
                sum = sum.add(Value(T).init(H_CI.at(0, j)).mul(C_coeff[(k - 1) * ndets + j]));
            }

            E[k] = sum;
        }

        if (k == order) break;

        for (1..ndets) |i| {
            var term1 = Value(T).fromFloat(0);
            var term3 = Value(T).fromFloat(0);

            for (0..ndets) |j| {
                term1 = term1.add(Value(T).init(H_CI.at(i, j)).mul(C_coeff[(k - 1) * ndets + j]));
            }

            const term2 = E0[i].mul(C_coeff[(k - 1) * ndets + i]);

            for (1..k + 1) |j| {
                term3 = term3.add(E[j].mul(C_coeff[(k - j) * ndets + i]));
            }

            C_coeff[k * ndets + i] = term1.sub(term2).sub(term3).div(E0[0].sub(E0[i]));
        }
    }

    const energies = try gpa.alloc(T, order + 1);
    errdefer gpa.free(energies);

    for (0..order + 1) |k| {
        energies[k] = E[k].val;
    }

    return energies;
}
