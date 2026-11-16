//! Evaluates analytic nuclear gradients for Møller-Plesset perturbation theory across orders and RI formulations.

const std = @import("std");

const Allocator = std.mem.Allocator;

const DualMatrix = @import("dual.zig").DualMatrix;
const DualTensor = @import("dual.zig").DualTensor;
const DualVector = @import("dual.zig").DualVector;
const ScalarDual = @import("dual.zig").ScalarDual;
const FittingDerivativeWorkspace = @import("density_fitting.zig").FittingDerivativeWorkspace;
const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;

const calculateCanonicalMp2Energy = @import("moller_plesset.zig").calculateCanonicalMp2Energy;
const calculateCorrelationEnergies = @import("moller_plesset.zig").calculateCorrelationEnergies;
const calculateRiMp2Energy = @import("moller_plesset.zig").calculateRiMp2Energy;

/// Evaluates Møller-Plesset analytical nuclear gradients across perturbation orders using dual numbers.
pub fn calculateMpGradient(comptime T: type, grads: []Matrix(T), order: usize, hfres: HartreeFockResult(T), gpa: Allocator) !void {
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
                try w.evaluate(hfres.ints.dJ.?.slice(coord).asMatrix(), hfres.ints.dI.?.slice(coord), hfres.ints.I.?, gpa);

                const B_d = DualTensor(T, 3).fromTensors(hfres.ints.B.?, w.dB_ao);

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
                try w.evaluate(hfres.ints.dJ.?.slice(coord).asMatrix(), hfres.ints.dI.?.slice(coord), hfres.ints.I.?, gpa);

                const B_d = DualTensor(T, 3).fromTensors(hfres.ints.B.?, w.dB_ao);

                corr = try calculateCorrelationEnergies(T, order, H_d, B_d, C_d, e_d, nocc, generalized, gpa);
            }

            if (deriv_ws == null) {
                const g_d = DualTensor(T, 4).fromTensors(hfres.ints.g.?, hfres.ints.dg.?.slice(coord));

                corr = try calculateCorrelationEnergies(T, order, H_d, g_d, C_d, e_d, nocc, generalized, gpa);
            }

            defer gpa.free(corr);

            for (2..order + 1) |k| {
                grads[k].ptr(coord / 3, coord % 3).* = corr[k].der;
            }
        }
    }
}
