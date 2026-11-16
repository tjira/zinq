//! Evaluates analytic nuclear gradients for Møller-Plesset perturbation theory across orders and RI formulations.

const std = @import("std");

const Allocator = std.mem.Allocator;

const DualMatrix = @import("dual.zig").DualMatrix;
const DualTensor = @import("dual.zig").DualTensor;
const DualVector = @import("dual.zig").DualVector;
const FittingDerivativeWorkspace = @import("density_fitting.zig").FittingDerivativeWorkspace;
const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;

const calculateCanonicalCorrelationEnergies = @import("moller_plesset.zig").calculateCanonicalCorrelationEnergies;
const calculateCanonicalMp2Energy = @import("moller_plesset.zig").calculateCanonicalMp2Energy;
const calculateRiCorrelationEnergies = @import("moller_plesset.zig").calculateRiCorrelationEnergies;
const calculateRiMp2Energy = @import("moller_plesset.zig").calculateRiMp2Energy;

/// Evaluates canonical Møller-Plesset analytical nuclear gradients across perturbation orders using dual numbers.
pub fn calculateCanonicalMpGradient(comptime T: type, grads: []Matrix(T), order: usize, hfres: HartreeFockResult(T), gpa: Allocator) !void {
    const generalized = hfres.ints.sys.nbf != hfres.C_mo.shape[0];

    const nocc = if (generalized) hfres.ints.sys.nel else hfres.ints.sys.nel / 2;

    for (0..3 * hfres.ints.sys.atoms.len) |coord| {
        const C_d = DualMatrix(T).fromMatrices(hfres.C_mo, hfres.dC_mo.?.slice(coord).asMatrix());

        const g_d = DualTensor(T, 4).fromTensors(hfres.ints.g.?, hfres.ints.dg.?.slice(coord));

        const e_d = DualVector(T).fromVectors(hfres.e_mo, hfres.de_mo.?.row(coord));

        if (order == 2) {
            const e_corr = try calculateCanonicalMp2Energy(T, g_d, C_d, e_d, nocc, generalized, gpa);

            grads[2].ptr(coord / 3, coord % 3).* = e_corr.der;
        }

        if (order > 2) {
            const H_d = DualMatrix(T).fromMatrices(hfres.ints.H.?, hfres.ints.dH.?.slice(coord).asMatrix());

            const corr = try calculateCanonicalCorrelationEnergies(T, order, H_d, g_d, C_d, e_d, nocc, generalized, gpa);
            defer gpa.free(corr);

            for (2..order + 1) |k| {
                grads[k].ptr(coord / 3, coord % 3).* = corr[k].der;
            }
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

    var deriv_ws = try FittingDerivativeWorkspace(T).init(J, B.shape, gpa);
    defer deriv_ws.deinit(gpa);

    for (0..3 * hfres.ints.sys.atoms.len) |coord| {
        try deriv_ws.evaluate(dJ.slice(coord).asMatrix(), dI.slice(coord), I_3c, gpa);

        const B_dual = DualTensor(T, 3).fromTensors(B, deriv_ws.dB_ao);

        const C_dual = DualMatrix(T).fromMatrices(hfres.C_mo, hfres.dC_mo.?.slice(coord).asMatrix());

        const e_dual = DualVector(T).fromVectors(hfres.e_mo, hfres.de_mo.?.row(coord));

        if (order == 2) {
            const e_corr = try calculateRiMp2Energy(T, B_dual, C_dual, e_dual, nocc, gpa);

            grads[2].ptr(coord / 3, coord % 3).* = e_corr.der;
        }

        if (order > 2) {
            const H_dual = DualMatrix(T).fromMatrices(hfres.ints.H.?, hfres.ints.dH.?.slice(coord).asMatrix());

            const corr = try calculateRiCorrelationEnergies(T, order, H_dual, B_dual, C_dual, e_dual, nocc, generalized, gpa);
            defer gpa.free(corr);

            for (2..order + 1) |k| {
                grads[k].ptr(coord / 3, coord % 3).* = corr[k].der;
            }
        }
    }
}
