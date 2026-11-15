//! Evaluates analytic nuclear gradients for Møller-Plesset perturbation theory across orders and RI formulations.

const std = @import("std");

const Allocator = std.mem.Allocator;

const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const FittingDerivativeWorkspace = @import("density_fitting.zig").FittingDerivativeWorkspace;
const ScalarDual = @import("dual.zig").ScalarDual;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

const calculateCanonicalCorrelationEnergies = @import("moller_plesset.zig").calculateCanonicalCorrelationEnergies;
const calculateRiCorrelationEnergies = @import("moller_plesset.zig").calculateRiCorrelationEnergies;
const calculateRiMp2Energy = @import("moller_plesset.zig").calculateRiMp2Energy;
const packDualTensor3c = @import("density_fitting.zig").packDualTensor3c;

/// Evaluates canonical Møller-Plesset analytical nuclear gradients across perturbation orders using dual numbers.
pub fn calculateCanonicalMpGradient(comptime T: type, grads: []Matrix(T), order: usize, hfres: HartreeFockResult(T), gpa: Allocator) !void {
    const generalized = hfres.ints.sys.nbf != hfres.C_mo.shape[0];

    const nocc = if (generalized) hfres.ints.sys.nel else hfres.ints.sys.nel / 2;

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

        const corr = try calculateCanonicalCorrelationEnergies(ScalarDual(T), order, H_d, g_d, C_d, e_d, nocc, generalized, gpa);
        defer gpa.free(corr);

        for (2..order + 1) |k| {
            grads[k].ptr(coord / 3, coord % 3).* = corr[k].der;
        }
    }
}

/// Evaluates resolution-of-the-identity MP nuclear gradients across perturbation orders using dual numbers.
pub fn calculateRiMpGradient(comptime T: type, grads: []Matrix(T), order: usize, hfres: HartreeFockResult(T), gpa: Allocator) !void {
    const generalized = hfres.ints.sys.nbf != hfres.C_mo.shape[0];

    const nocc = if (generalized) hfres.ints.sys.nel else hfres.ints.sys.nel / 2;

    const B = hfres.ints.B.?;

    const J = hfres.ints.J orelse unreachable;

    const dI = hfres.ints.dI orelse unreachable;
    const dJ = hfres.ints.dJ orelse unreachable;

    const I_3c = hfres.ints.I orelse unreachable;

    const nmo = hfres.C_mo.shape[0];

    var deriv_ws = try FittingDerivativeWorkspace(T).init(J, B.shape, gpa);
    defer deriv_ws.deinit(gpa);

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
        try deriv_ws.evaluate(dJ.slice(coord).asMatrix(), dI.slice(coord), I_3c, gpa);

        packDualTensor3c(T, &B_d, B, deriv_ws.dB_ao);

        for (0..nmo) |mu| for (0..nmo) |p| {
            C_d.ptr(mu, p).* = ScalarDual(T).init(hfres.C_mo.at(mu, p), hfres.dC_mo.?.at(.{ coord, mu, p }));
        };

        for (0..nmo) |p| {
            e_d.ptr(p).* = ScalarDual(T).init(hfres.e_mo.at(p), hfres.de_mo.?.at(coord, p));
        }

        if (order == 2) {
            const e_corr = try calculateRiMp2Energy(ScalarDual(T), B_d, C_d, e_d, nocc, gpa);

            grads[2].ptr(coord / 3, coord % 3).* = e_corr.der;
        }

        if (order > 2) {
            for (0..h_dim) |mu| for (0..h_dim) |nu| {
                H_d.ptr(mu, nu).* = ScalarDual(T).init(hfres.ints.H.?.at(mu, nu), hfres.ints.dH.?.at(.{ coord, mu, nu }));
            };

            const corr = try calculateRiCorrelationEnergies(ScalarDual(T), order, H_d, B_d, C_d, e_d, nocc, generalized, gpa);
            defer gpa.free(corr);

            for (2..order + 1) |k| {
                grads[k].ptr(coord / 3, coord % 3).* = corr[k].der;
            }
        }
    }
}
