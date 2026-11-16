//! Evaluates analytic nuclear gradients for Configuration Interaction and RI-CI electronic states.

const std = @import("std");

const Allocator = std.mem.Allocator;

const DualMatrix = @import("dual.zig").DualMatrix;
const DualTensor = @import("dual.zig").DualTensor;
const FittingDerivativeWorkspace = @import("density_fitting.zig").FittingDerivativeWorkspace;
const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;

const slater = @import("slater_condon.zig").slater;
const transformInts = @import("integral_transform.zig").transformInts;
const transformIntsRi = @import("integral_transform.zig").transformIntsRi;

/// Computes analytic nuclear gradients of CI state energies from two-electron integrals or fitting tensors.
pub fn calculateCiGradient(comptime T: type, hfres: HartreeFockResult(T), C_ci: Matrix(T), dets: []const []const usize, states: []const u32, gpa: Allocator) ![]Matrix(T) {
    const generalized = hfres.ints.sys.nbf != hfres.C_mo.shape[0];

    var grads = try gpa.alloc(Matrix(T), states.len);
    errdefer gpa.free(grads);

    var grads_done: usize = 0;
    errdefer for (0..grads_done) |k| grads[k].deinit(gpa);

    for (0..states.len) |idx| {
        grads[idx] = try hfres.ints.sys.nuclearRepulsionGradient(gpa);

        grads_done += 1;
    }

    var deriv_ws: ?FittingDerivativeWorkspace(T) = null;

    if (hfres.ints.B) |B| {
        deriv_ws = try FittingDerivativeWorkspace(T).init(hfres.ints.J.?, B.shape, gpa);
    }

    defer if (deriv_ws) |*w| w.deinit(gpa);

    for (0..3 * hfres.ints.sys.atoms.len) |coord| {
        const C_d = DualMatrix(T).fromMatrices(hfres.C_mo, hfres.dC_mo.?.slice(coord).asMatrix());

        const H_d = DualMatrix(T).fromMatrices(hfres.ints.H.?, hfres.ints.dH.?.slice(coord).asMatrix());

        var h_ms_d: DualMatrix(T), var g_ms_d: DualTensor(T, 4) = .{ undefined, undefined };

        if (deriv_ws) |*w| {
            try w.evaluate(hfres.ints.dJ.?.slice(coord).asMatrix(), hfres.ints.dI.?.slice(coord), hfres.ints.I.?, gpa);

            const B_d = DualTensor(T, 3).fromTensors(hfres.ints.B.?, w.dB_ao);

            h_ms_d, g_ms_d = try transformIntsRi(T, C_d, H_d, B_d, generalized, gpa);
        }

        if (deriv_ws == null) {
            const g_d = DualTensor(T, 4).fromTensors(hfres.ints.g.?, hfres.ints.dg.?.slice(coord));

            h_ms_d, g_ms_d = try transformInts(T, C_d, H_d, g_d, generalized, gpa);
        }

        defer {
            h_ms_d.deinit(gpa);
            g_ms_d.deinit(gpa);
        }

        for (0..dets.len) |m| for (m..dets.len) |n| {
            const factor: T = if (m == n) 1 else 2;

            const der = slater(T, dets[m], dets[n], h_ms_d.der, g_ms_d.der);

            for (states, 0..) |k, k_idx| {
                grads[k_idx].ptr(coord / 3, coord % 3).* += factor * C_ci.at(m, k) * C_ci.at(n, k) * der;
            }
        };
    }

    return grads;
}
