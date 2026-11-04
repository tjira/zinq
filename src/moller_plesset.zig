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

const generateDets = @import("determinant.zig").generateDets;
const printf = @import("read_write.zig").printf;
const slater = @import("slater_condon.zig").slater;
const transformInts = @import("integral_transform.zig").transformInts;

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

/// Computes the nuclear gradient of Møller-Plesset correlation energies using dual number differentiation.
pub fn calculateNuclearGradient(comptime T: type, order: usize, hfres: HartreeFockResult(T), gpa: Allocator) ![]Matrix(T) {
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

    for (0..3 * hfres.ints.sys.atoms.len) |coord| {
        var H_d = try Matrix(ScalarDual(T)).init(nbf, nbf, gpa);
        defer H_d.deinit(gpa);

        for (0..nbf) |mu| for (0..nbf) |nu| {
            H_d.ptr(mu, nu).* = ScalarDual(T).init(hfres.ints.H.?.at(mu, nu), hfres.ints.dH.?.at(.{ coord, mu, nu }));
        };

        var C_d = try Matrix(ScalarDual(T)).init(nbf, nbf, gpa);
        defer C_d.deinit(gpa);

        for (0..nbf) |mu| for (0..nbf) |p| {
            C_d.ptr(mu, p).* = ScalarDual(T).init(hfres.C_mo.at(mu, p), hfres.dC_mo.?.at(.{ coord, mu, p }));
        };

        var e_d = try Vector(ScalarDual(T)).init(nbf, gpa);
        defer e_d.deinit(gpa);

        for (0..nbf) |p| {
            e_d.ptr(p).* = ScalarDual(T).init(hfres.e_mo.at(p), hfres.de_mo.?.at(coord, p));
        }

        var g_d = try Tensor(ScalarDual(T), 4).init(hfres.ints.g.?.shape, gpa);
        defer g_d.deinit(gpa);

        for (0..nbf) |mu| for (0..nbf) |lambda| for (0..nbf) |nu| for (0..nbf) |sigma| {
            const idx = .{ mu, lambda, nu, sigma };

            const dg_x = hfres.ints.dg.?.at(.{ coord, mu, lambda, nu, sigma });

            g_d.ptr(idx).* = ScalarDual(T).init(hfres.ints.g.?.at(idx), dg_x);
        };

        const corr = try calculateCorrelationEnergies(ScalarDual(T), order, H_d, g_d, C_d, e_d, nocc, generalized, gpa);
        defer gpa.free(corr);

        for (2..order + 1) |k| {
            grads[k].ptr(coord / 3, coord % 3).* = corr[k].der;
        }
    }

    return grads;
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

            grads = try calculateNuclearGradient(T, opt.order, hfres, gpa);
        },
        .numeric => {},
    };

    energy[0] = hfres.energy[0];

    const generalized = opt.hartree_fock.generalized;

    const nocc = if (generalized) hfres.ints.sys.nel else hfres.ints.sys.nel / 2;

    const H_ao, const g_ao, const C_mo, const e_mo = .{ hfres.ints.H.?, hfres.ints.g.?, hfres.C_mo, hfres.e_mo };

    const energies = try calculateCorrelationEnergies(T, opt.order, H_ao, g_ao, C_mo, e_mo, nocc, generalized, gpa);
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

/// Evaluates electronic Hamiltonian matrix elements across Slater determinants using Slater-Condon rules.
fn buildHamiltonianMatrix(comptime T: type, dets: []const []const usize, h_ms: Matrix(T), g_ms: Tensor(T, 4), gpa: Allocator) !Matrix(T) {
    var H_dets = try Matrix(T).init(dets.len, dets.len, gpa);
    errdefer H_dets.deinit(gpa);

    for (0..dets.len) |i| for (i..dets.len) |j| {
        const val = slater(T, dets[i], dets[j], h_ms, g_ms);

        H_dets.ptr(i, j).* = val;
        H_dets.ptr(j, i).* = val;
    };

    return H_dets;
}

/// Computes Møller-Plesset correlation energy corrections up to the requested perturbation order.
fn calculateCorrelationEnergies(comptime T: type, order: usize, H: Matrix(T), g: Tensor(T, 4), C: Matrix(T), e: Vector(T), nocc: usize, generalized: bool, gpa: Allocator) ![]T {
    const nsp, const nel = if (generalized) .{ C.shape[0], nocc } else .{ 2 * C.shape[0], 2 * nocc };

    var h_ms, var g_ms = try transformInts(T, C, H, g, generalized, gpa);

    defer {
        h_ms.deinit(gpa);
        g_ms.deinit(gpa);
    }

    var e_ms = try Vector(T).init(nsp, gpa);
    defer e_ms.deinit(gpa);

    if (generalized) for (0..nsp) |p| {
        e_ms.ptr(p).* = e.at(p);
    };

    if (!generalized) for (0..C.shape[0]) |p| {
        e_ms.ptr(2 * p + 0).* = e.at(p);
        e_ms.ptr(2 * p + 1).* = e.at(p);
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

    var H_dets = try buildHamiltonianMatrix(T, dets.items, h_ms, g_ms, gpa);
    defer H_dets.deinit(gpa);

    const E0 = try gpa.alloc(Value(T), dets.items.len);
    defer gpa.free(E0);

    calculateZerothOrderEnergies(T, E0, dets.items, e_ms);

    return try solveRayleighSchrodinger(T, H_dets, E0, order, gpa);
}

/// Evaluates unperturbed zeroth-order energies of Slater determinants under the unperturbed Fock Hamiltonian.
fn calculateZerothOrderEnergies(comptime T: type, E0: []Value(T), dets: []const []const usize, e_ms: Vector(T)) void {
    for (0..dets.len) |i| {
        var sum = Value(T).fromFloat(0);

        for (dets[i]) |p| {
            sum = sum.add(Value(T).init(e_ms.at(p)));
        }

        E0[i] = sum;
    }
}

/// Computes the m-th order Rayleigh-Schrödinger correlation energy correction from Hamiltonian couplings.
fn evaluateEnergyCorrection(comptime T: type, m: usize, H_dets: Matrix(T), c_dets: []const Value(T)) Value(T) {
    const ndets = H_dets.nrow();

    var sum = Value(T).fromFloat(0);

    for (1..ndets) |j| {
        sum = sum.add(Value(T).init(H_dets.at(0, j)).mul(c_dets[(m - 1) * ndets + j]));
    }

    return sum;
}

/// Solves Rayleigh-Schrödinger perturbation equations iteratively in the Slater determinant basis.
fn solveRayleighSchrodinger(comptime T: type, H_dets: Matrix(T), E0: []const Value(T), order: usize, gpa: Allocator) ![]T {
    const ndets = E0.len;

    var c_dets = try gpa.alloc(Value(T), (order + 1) * ndets);
    defer gpa.free(c_dets);

    @memset(c_dets, Value(T).fromFloat(0));

    var E = try gpa.alloc(Value(T), order + 1);
    defer gpa.free(E);

    E[0], c_dets[0] = .{ E0[0], Value(T).fromFloat(1) };

    for (1..order + 1) |m| {
        if (m == 1) {
            E[1] = Value(T).init(H_dets.at(0, 0)).sub(E0[0]);
        }

        if (m > 1) {
            E[m] = evaluateEnergyCorrection(T, m, H_dets, c_dets);
        }

        if (m == order) break;

        updateAmplitudes(T, m, H_dets, E0, E, c_dets);
    }

    const energies = try gpa.alloc(T, order + 1);
    errdefer gpa.free(energies);

    for (0..order + 1) |m| {
        energies[m] = E[m].val;
    }

    return energies;
}

/// Evaluates m-th order configuration mixing amplitudes using Rayleigh-Schrödinger perturbation theory.
fn updateAmplitudes(comptime T: type, m: usize, H_dets: Matrix(T), E0: []const Value(T), E: []const Value(T), c_dets: []Value(T)) void {
    const ndets = E0.len;

    for (1..ndets) |i| {
        var h_psi = Value(T).fromFloat(0);
        var e_psi = Value(T).fromFloat(0);

        for (0..ndets) |j| {
            h_psi = h_psi.add(Value(T).init(H_dets.at(i, j)).mul(c_dets[(m - 1) * ndets + j]));
        }

        const h0_psi = E0[i].mul(c_dets[(m - 1) * ndets + i]);

        for (1..m + 1) |j| {
            e_psi = e_psi.add(E[j].mul(c_dets[(m - j) * ndets + i]));
        }

        c_dets[m * ndets + i] = h_psi.sub(h0_psi).sub(e_psi).div(E0[0].sub(E0[i]));
    }
}