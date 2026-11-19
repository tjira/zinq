//! Computes Møller-Plesset correlation energy corrections to a Hartree-Fock reference wavefunction.

const std = @import("std");

const cblas = @import("cimport.zig").cblas;

const Allocator = std.mem.Allocator;

const DualMatrix = @import("dual.zig").DualMatrix;
const DualTensor = @import("dual.zig").DualTensor;
const DualVector = @import("dual.zig").DualVector;
const FittingDerivativeWorkspace = @import("density_fitting.zig").FittingDerivativeWorkspace;
const FrequencyOptions = @import("frequency_analysis.zig").Options;
const HartreeFockOptions = @import("hartree_fock.zig").Options;
const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const OptimizeOptions = @import("molecular_optimization.zig").OptimizeOptions(GradientOptions);
const ScalarDual = @import("dual.zig").ScalarDual;
const Tensor = @import("tensor.zig").Tensor;
const Value = @import("value.zig").Value;
const Vector = @import("tensor.zig").Vector;

const ao2mo_oovv = @import("integral_transform.zig").ao2mo_oovv;
const ao2mo_ov_3c = @import("integral_transform.zig").ao2mo_ov_3c;
const buildFittingIntegrals = @import("density_fitting.zig").buildFittingIntegrals;
const buildHamiltonianMatrix = @import("slater_condon.zig").buildHamiltonianMatrix;
const gemm = @import("linear_algebra.zig").gemm;
const generateDets = @import("determinant.zig").generateDets;
const isDual = @import("value.zig").isDual;
const isDualMatrix = @import("value.zig").isDualMatrix;
const isDualTensor = @import("value.zig").isDualTensor;
const printf = @import("read_write.zig").printf;
const transformInts = @import("integral_transform.zig").transformInts;

const AU2CM = @import("constants.zig").AU2CM;

/// Parameters governing the Møller-Plesset perturbation theory calculation and its derivatives.
pub const Options = struct {
    aux_basis: ?[]const u8 = null,
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
    const grads = try gpa.alloc(Matrix(T), order + 1);
    errdefer gpa.free(grads);

    var initialized: usize = 0;
    errdefer for (0..initialized) |k| grads[k].deinit(gpa);

    for (0..order + 1) |k| {
        grads[k] = try Matrix(T).initZero(hfres.ints.sys.atoms.len, 3, gpa);

        initialized += 1;
    }

    try calculateMpGradient(T, grads, order, hfres, gpa);

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

    const analytic_grad = if (opt.gradient) |g| g == .analytic else false;

    if (opt.aux_basis) |aux| {
        try buildFittingIntegrals(T, io, &hfres.ints, aux, analytic_grad, opt.nthreads, gpa);
    }

    if (opt.gradient) |g| switch (g) {
        .analytic => {
            var g0 = try hfres.grad[0].clone(gpa);
            errdefer g0.deinit(gpa);

            const g_arr = try gpa.alloc(Matrix(T), 1);
            errdefer gpa.free(g_arr);

            grads = try calculateNuclearGradient(T, opt.order, hfres, gpa);

            g_arr[0] = g0;

            grad = g_arr;
        },
        .numeric => {},
    };

    energy[0] = hfres.energy[0];

    try evaluateMp(T, io, opt, hfres, energy, grad, grads, log, gpa);

    return .{ .hartree_fock = hfres, .energy = energy, .grad = grad, .hess = &.{} };
}

/// Evaluates canonical second-order Møller-Plesset correlation energy from two-electron integrals.
pub fn calculateCanonicalMp2Energy(comptime T: type, g: anytype, C: anytype, e: anytype, nocc: usize, gen: bool, gpa: Allocator) !if (isDual(@TypeOf(g))) ScalarDual(T) else T {
    const is_dual = comptime isDualTensor(@TypeOf(g));

    const OutType = if (is_dual) ScalarDual(T) else T;

    const nvir = C.shape[0] - nocc;

    var g_mo = try @TypeOf(g).initZero(.{ nocc, nocc, nvir, nvir }, gpa);
    defer g_mo.deinit(gpa);

    try ao2mo_oovv(T, &g_mo, g, C, nocc, gpa);

    var e_corr = Value(OutType).fromFloat(0);

    for (0..nocc) |i| for (0..nocc) |j| {
        const val_ei = Value(OutType).init(e.at(i));
        const val_ej = Value(OutType).init(e.at(j));

        for (0..nvir) |a| {
            const val_ea = Value(OutType).init(e.at(nocc + a));

            for (0..nvir) |b| {
                const val_eb = Value(OutType).init(e.at(nocc + b));

                const denom = val_ei.add(val_ej).sub(val_ea).sub(val_eb);

                const g_ijab = Value(OutType).init(g_mo.at(.{ i, j, a, b }));
                const g_ijba = Value(OutType).init(g_mo.at(.{ i, j, b, a }));

                var num = Value(OutType).fromFloat(0);

                if (gen) {
                    const diff = g_ijab.sub(g_ijba);
                    num = diff.mul(diff).muls(0.25);
                }

                if (!gen) {
                    num = g_ijab.mul(g_ijab.muls(2).sub(g_ijba));
                }

                e_corr = e_corr.add(num.div(denom));
            }
        }
    };

    return e_corr.val;
}

/// Computes Møller-Plesset correlation energy corrections up to the requested perturbation order.
pub fn calculateCorrelationEnergies(comptime T: type, ord: usize, ints: anytype, mo: anytype, gpa: Allocator) !if (isDualTensor(@TypeOf(ints[1]))) []ScalarDual(T) else []T {
    const H, const g = ints;

    const C, const e, const nocc, const gen = mo;

    const is_dual = comptime isDualTensor(@TypeOf(g));

    const OutType = if (is_dual) ScalarDual(T) else T;

    if (ord == 2) {
        const energies = try gpa.alloc(OutType, 3);
        errdefer gpa.free(energies);

        energies[0] = Value(OutType).fromFloat(0).val;
        energies[1] = Value(OutType).fromFloat(0).val;

        if (comptime isDualTensor(@TypeOf(g))) {
            if (g.val.shape.len == 4) {
                energies[2] = try calculateCanonicalMp2Energy(T, g, C, e, nocc, gen, gpa);
            }

            if (g.val.shape.len == 3) {
                energies[2] = try calculateRiMp2Energy(T, g, C, e, nocc, gpa);
            }
        }

        if (comptime !isDualTensor(@TypeOf(g))) {
            if (g.shape.len == 4) {
                energies[2] = try calculateCanonicalMp2Energy(T, g, C, e, nocc, gen, gpa);
            }

            if (g.shape.len == 3) {
                energies[2] = try calculateRiMp2Energy(T, g, C, e, nocc, gpa);
            }
        }

        return energies;
    }

    const nsp, const nel = if (gen) .{ C.shape[0], nocc } else .{ 2 * C.shape[0], 2 * nocc };

    var h_ms, var g_ms = try transformInts(T, .{ C, gen }, H, g, gpa);

    defer {
        h_ms.deinit(gpa);
        g_ms.deinit(gpa);
    }

    return try solveHigherOrderRspt(T, OutType, ord, nel, nsp, gen, h_ms, g_ms, e, gpa);
}

/// Evaluates the second-order Møller-Plesset correlation energy using resolution-of-the-identity density fitting.
pub fn calculateRiMp2Energy(comptime T: type, B_ao: anytype, C: anytype, e: anytype, nocc: usize, gpa: Allocator) !if (isDualTensor(@TypeOf(B_ao))) ScalarDual(T) else T {
    const is_dual, const generalized = .{ comptime isDualTensor(@TypeOf(B_ao)), C.shape[0] != B_ao.shape[0] };

    const OutType = if (is_dual) ScalarDual(T) else T;

    const naux, const nvir = .{ B_ao.shape[2], C.shape[1] - nocc };

    var B_mo = try @TypeOf(B_ao).initZero(.{ nocc, nvir, naux }, gpa);
    defer B_mo.deinit(gpa);

    try ao2mo_ov_3c(T, &B_mo, B_ao, C, gpa);

    const MatType = if (is_dual) DualMatrix(T) else Matrix(T);

    var M_ij = try MatType.initZero(nvir, nvir, gpa);
    defer M_ij.deinit(gpa);

    var e_corr = Value(OutType).fromFloat(0);

    for (0..nocc) |i| {
        for (0..nocc) |j| {
            gemm(T, &M_ij, B_mo.slice(i).asMatrix(), B_mo.slice(j).asMatrix(), 1, 0, false, true);

            const val_ei = Value(OutType).init(e.at(i));
            const val_ej = Value(OutType).init(e.at(j));

            for (0..nvir) |a| {
                const val_ea = Value(OutType).init(e.at(nocc + a));

                for (0..nvir) |b| {
                    const val_eb = Value(OutType).init(e.at(nocc + b));

                    const denom = val_ei.add(val_ej).sub(val_ea).sub(val_eb);

                    const g_ijab = Value(OutType).init(M_ij.at(a, b));
                    const g_ijba = Value(OutType).init(M_ij.at(b, a));

                    var num = Value(OutType).fromFloat(0);

                    if (generalized) {
                        const diff = g_ijab.sub(g_ijba);
                        num = diff.mul(diff).muls(0.25);
                    }

                    if (!generalized) {
                        num = g_ijab.mul(g_ijab.muls(2).sub(g_ijba));
                    }

                    e_corr = e_corr.add(num.div(denom));
                }
            }
        }
    }

    return e_corr.val;
}

/// Evaluates unperturbed zeroth-order energies of Slater determinants under the unperturbed Fock Hamiltonian.
fn calculateZerothOrderEnergies(comptime OutType: type, E0: []Value(OutType), dets: []const []const usize, e_ms: anytype) void {
    for (0..dets.len) |i| {
        var sum = Value(OutType).fromFloat(0);

        for (dets[i]) |p| {
            sum = sum.add(Value(OutType).init(e_ms.at(p)));
        }

        E0[i] = sum;
    }
}

/// Evaluates Møller-Plesset correlation energies and nuclear gradients across perturbation orders.
fn evaluateMp(comptime T: type, io: std.Io, opt: Options, hfres: HartreeFockResult(T), energy: []T, grad: []Matrix(T), grads: ?[]Matrix(T), log: bool, gpa: Allocator) !void {
    const generalized = opt.hartree_fock.generalized;

    const nocc = if (generalized) hfres.ints.sys.nel else hfres.ints.sys.nel / 2;

    var energies: []T = undefined;

    const mo_ctx = .{ hfres.C_mo, hfres.e_mo, nocc, generalized };

    if (hfres.ints.B) |B| {
        energies = try calculateCorrelationEnergies(T, opt.order, .{ hfres.ints.H.?, B }, mo_ctx, gpa);
    }

    if (hfres.ints.B == null) {
        energies = try calculateCorrelationEnergies(T, opt.order, .{ hfres.ints.H.?, hfres.ints.g.? }, mo_ctx, gpa);
    }

    defer gpa.free(energies);

    for (2..opt.order + 1) |i| {
        energy[0] += energies[i];

        if (log) {
            try printf(io, "\nMP{d:02} TOTAL ENERGY: {d:.14} Eh\n", .{ i, energy[0] });
        }

        if (grads) |g| {
            for (0..grad[0].nrow()) |j| for (0..grad[0].ncol()) |k| {
                grad[0].ptr(j, k).* += g[i].at(j, k);
            };
        }
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

/// Solves higher-order Rayleigh-Schrödinger perturbation theory across Slater determinants.
fn solveHigherOrderRspt(comptime T: type, comptime Out: type, order: usize, nel: usize, nsp: usize, gen: bool, h_ms: anytype, g_ms: anytype, e: anytype, gpa: Allocator) ![]Out {
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

    const is_dual = comptime isDualMatrix(@TypeOf(h_ms));

    const VecType = if (is_dual) DualVector(T) else Vector(T);

    var e_ms = try VecType.init(nsp, gpa);
    defer e_ms.deinit(gpa);

    if (gen) for (0..nsp) |p| {
        e_ms.set(p, e.at(p));
    };

    if (!gen) for (0..nsp / 2) |p| {
        e_ms.set(2 * p + 0, e.at(p));
        e_ms.set(2 * p + 1, e.at(p));
    };

    var H_dets = try buildHamiltonianMatrix(T, dets.items, h_ms, g_ms, gpa);
    defer H_dets.deinit(gpa);

    const E0 = try gpa.alloc(Value(Out), dets.items.len);
    defer gpa.free(E0);

    calculateZerothOrderEnergies(Out, E0, dets.items, e_ms);

    return try solveRayleighSchrodinger(Out, H_dets, E0, order, gpa);
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

/// Evaluates Møller-Plesset analytical nuclear gradients across perturbation orders using dual numbers.
fn calculateMpGradient(comptime T: type, grads: []Matrix(T), order: usize, hfres: HartreeFockResult(T), gpa: Allocator) !void {
    const generalized = hfres.ints.sys.nbf != hfres.C_mo.shape[0];

    const nocc = if (generalized) hfres.ints.sys.nel else hfres.ints.sys.nel / 2;

    var deriv_ws: ?FittingDerivativeWorkspace(T) = null;

    if (hfres.ints.B) |B| {
        deriv_ws = try FittingDerivativeWorkspace(T).init(hfres.ints.J.?, B.shape, gpa);
    }

    defer if (deriv_ws) |*w| w.deinit(gpa);

    for (0..3 * hfres.ints.sys.atoms.len) |coord| {
        const C_d = DualMatrix(T).fromMatrices(hfres.C_mo, hfres.dC_mo.?.slice(coord).asMatrix());

        const e_d = DualVector(T).fromVectors(hfres.e_mo, hfres.de_mo.?.row(coord));

        if (order == 2) {
            var der: T = undefined;

            if (deriv_ws) |*w| {
                const B_d = try w.evaluateCoord(hfres.ints, coord, gpa);

                der = (try calculateRiMp2Energy(T, B_d, C_d, e_d, nocc, gpa)).der;
            }

            if (deriv_ws == null) {
                const g_d = DualTensor(T, 4).fromTensors(hfres.ints.g.?, hfres.ints.dg.?.slice(coord));

                der = (try calculateCanonicalMp2Energy(T, g_d, C_d, e_d, nocc, generalized, gpa)).der;
            }

            grads[2].ptr(coord / 3, coord % 3).* = der;
        }

        if (order > 2) {
            const H_d = DualMatrix(T).fromMatrices(hfres.ints.H.?, hfres.ints.dH.?.slice(coord).asMatrix());

            var corr: []ScalarDual(T) = undefined;

            if (deriv_ws) |*w| {
                const B_d = try w.evaluateCoord(hfres.ints, coord, gpa);

                corr = try calculateCorrelationEnergies(T, order, .{ H_d, B_d }, .{ C_d, e_d, nocc, generalized }, gpa);
            }

            if (deriv_ws == null) {
                const g_d = DualTensor(T, 4).fromTensors(hfres.ints.g.?, hfres.ints.dg.?.slice(coord));

                corr = try calculateCorrelationEnergies(T, order, .{ H_d, g_d }, .{ C_d, e_d, nocc, generalized }, gpa);
            }

            defer gpa.free(corr);

            for (2..order + 1) |k| {
                grads[k].ptr(coord / 3, coord % 3).* = corr[k].der;
            }
        }
    }
}
