//! Evaluates analytic nuclear gradients for Møller-Plesset perturbation theory across orders and RI formulations.

const std = @import("std");

const Allocator = std.mem.Allocator;

const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const MetricDecomposition = @import("density_fitting.zig").MetricDecomposition;
const ScalarDual = @import("dual.zig").ScalarDual;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

const buildMetricDerivKernel = @import("density_fitting.zig").buildMetricDerivKernel;
const calculateCorrelationEnergies = @import("moller_plesset.zig").calculateCorrelationEnergies;
const calculateCorrelationEnergiesRi = @import("moller_plesset.zig").calculateCorrelationEnergiesRi;
const calculateRiMp2Energy = @import("moller_plesset.zig").calculateRiMp2Energy;
const evaluateFittingTensorDeriv = @import("density_fitting.zig").evaluateFittingTensorDeriv;
const evaluateMetricDeriv = @import("density_fitting.zig").evaluateMetricDeriv;

/// Evaluates canonical Møller-Plesset analytical nuclear gradients across perturbation orders using dual numbers.
pub fn calculateCanonicalMpGradient(comptime T: type, grads: []Matrix(T), order: usize, hfres: HartreeFockResult(T), nocc: usize, generalized: bool, gpa: Allocator) !void {
    const nmo = hfres.C_mo.shape[0];

    var H_d = try Matrix(ScalarDual(T)).init(nmo, nmo, gpa);
    defer H_d.deinit(gpa);

    var C_d = try Matrix(ScalarDual(T)).init(nmo, nmo, gpa);
    defer C_d.deinit(gpa);

    var e_d = try Vector(ScalarDual(T)).init(nmo, gpa);
    defer e_d.deinit(gpa);

    var g_d = try Tensor(ScalarDual(T), 4).init(hfres.ints.g.?.shape, gpa);
    defer g_d.deinit(gpa);

    const nao = g_d.shape[0];

    for (0..3 * hfres.ints.sys.atoms.len) |coord| {
        for (0..nmo) |mu| for (0..nmo) |nu| {
            H_d.ptr(mu, nu).* = ScalarDual(T).init(hfres.ints.H.?.at(mu, nu), hfres.ints.dH.?.at(.{ coord, mu, nu }));
        };

        for (0..nmo) |mu| for (0..nmo) |p| {
            C_d.ptr(mu, p).* = ScalarDual(T).init(hfres.C_mo.at(mu, p), hfres.dC_mo.?.at(.{ coord, mu, p }));
        };

        for (0..nmo) |p| {
            e_d.ptr(p).* = ScalarDual(T).init(hfres.e_mo.at(p), hfres.de_mo.?.at(coord, p));
        }

        for (0..nao) |mu| for (0..nao) |lambda| for (0..nao) |nu| for (0..nao) |sigma| {
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
}

/// Evaluates resolution-of-the-identity MP nuclear gradients across perturbation orders using dual numbers.
pub fn calculateRiMpGradient(comptime T: type, grads: []Matrix(T), order: usize, hfres: HartreeFockResult(T), nocc: usize, generalized: bool, gpa: Allocator) !void {
    const B = hfres.ints.B.?;

    const J = hfres.ints.J orelse unreachable;

    const dI = hfres.ints.dI orelse unreachable;
    const dJ = hfres.ints.dJ orelse unreachable;

    const I_3c = hfres.ints.I orelse unreachable;

    const nbf, const naux, const nmo = .{ B.shape[0], B.shape[2], hfres.C_mo.shape[0] };

    var decomp = try MetricDecomposition(T).init(J, gpa);
    defer decomp.deinit(gpa);

    var M_kernel = try buildMetricDerivKernel(T, decomp.W, gpa);
    defer M_kernel.deinit(gpa);

    var dJ_inv_half = try Matrix(T).init(naux, naux, gpa);
    defer dJ_inv_half.deinit(gpa);

    var dB_ao = try Tensor(T, 3).initZero(B.shape, gpa);
    defer dB_ao.deinit(gpa);

    var B_d = try Tensor(ScalarDual(T), 3).init(B.shape, gpa);
    defer B_d.deinit(gpa);

    var C_d = try Matrix(ScalarDual(T)).init(nmo, nmo, gpa);
    defer C_d.deinit(gpa);

    var e_d = try Vector(ScalarDual(T)).init(nmo, gpa);
    defer e_d.deinit(gpa);

    const h_dim = hfres.ints.H.?.nrow();

    var H_d = if (order > 2) try Matrix(ScalarDual(T)).init(h_dim, h_dim, gpa) else undefined;
    defer if (order > 2) H_d.deinit(gpa);

    for (0..3 * hfres.ints.sys.atoms.len) |coord| {
        try evaluateMetricDeriv(T, &dJ_inv_half, dJ.slice(coord).asMatrix(), decomp.U, M_kernel, gpa);

        evaluateFittingTensorDeriv(T, &dB_ao, dI.slice(coord), I_3c, decomp.J_inv_half, dJ_inv_half);

        for (0..nbf) |mu| for (0..nbf) |nu| for (0..naux) |P| {
            B_d.ptr(.{ mu, nu, P }).* = ScalarDual(T).init(B.at(.{ mu, nu, P }), dB_ao.at(.{ mu, nu, P }));
        };

        for (0..nmo) |mu| for (0..nmo) |p| {
            C_d.ptr(mu, p).* = ScalarDual(T).init(hfres.C_mo.at(mu, p), hfres.dC_mo.?.at(.{ coord, mu, p }));
        };

        for (0..nmo) |p| {
            e_d.ptr(p).* = ScalarDual(T).init(hfres.e_mo.at(p), hfres.de_mo.?.at(coord, p));
        }

        if (order == 2) {
            const e_corr = try calculateRiMp2Energy(ScalarDual(T), B_d, C_d, e_d, nocc, generalized, gpa);

            grads[2].ptr(coord / 3, coord % 3).* = e_corr.der;
        }

        if (order > 2) {
            for (0..h_dim) |mu| for (0..h_dim) |nu| {
                H_d.ptr(mu, nu).* = ScalarDual(T).init(hfres.ints.H.?.at(mu, nu), hfres.ints.dH.?.at(.{ coord, mu, nu }));
            };

            const corr = try calculateCorrelationEnergiesRi(ScalarDual(T), order, H_d, B_d, C_d, e_d, nocc, generalized, gpa);
            defer gpa.free(corr);

            for (2..order + 1) |k| {
                grads[k].ptr(coord / 3, coord % 3).* = corr[k].der;
            }
        }
    }
}
