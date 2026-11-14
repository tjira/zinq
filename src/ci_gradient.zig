//! Evaluates analytic nuclear gradients for Configuration Interaction and RI-CI electronic states.

const std = @import("std");

const Allocator = std.mem.Allocator;

const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const MetricDecomposition = @import("density_fitting.zig").MetricDecomposition;
const ScalarDual = @import("dual.zig").ScalarDual;
const Tensor = @import("tensor.zig").Tensor;

const buildMetricDerivKernel = @import("density_fitting.zig").buildMetricDerivKernel;
const evaluateFittingTensorDeriv = @import("density_fitting.zig").evaluateFittingTensorDeriv;
const evaluateMetricDeriv = @import("density_fitting.zig").evaluateMetricDeriv;
const slater = @import("slater_condon.zig").slater;
const transformInts = @import("integral_transform.zig").transformInts;
const transformIntsRi = @import("integral_transform.zig").transformIntsRi;

/// Computes analytic nuclear gradients of CI state energies from four-center two-electron integrals.
pub fn calculateCanonicalCiGradient(comptime T: type, hfres: HartreeFockResult(T), C_ci: Matrix(T), dets: []const []const usize, states: []const u32, gpa: Allocator) ![]Matrix(T) {
    const generalized = hfres.ints.sys.nbf != hfres.C_mo.shape[0];

    var grads = try gpa.alloc(Matrix(T), states.len);
    errdefer gpa.free(grads);

    var grads_done: usize = 0;
    errdefer for (0..grads_done) |k| grads[k].deinit(gpa);

    for (0..states.len) |idx| {
        grads[idx] = try hfres.ints.sys.nuclearRepulsionGradient(gpa);

        grads_done += 1;
    }

    var H_d = try Matrix(ScalarDual(T)).init(hfres.ints.H.?.nrow(), hfres.ints.H.?.ncol(), gpa);
    defer H_d.deinit(gpa);

    var g_d = try Tensor(ScalarDual(T), 4).init(hfres.ints.g.?.shape, gpa);
    defer g_d.deinit(gpa);

    var C_d = try Matrix(ScalarDual(T)).init(hfres.C_mo.nrow(), hfres.C_mo.ncol(), gpa);
    defer C_d.deinit(gpa);

    for (0..3 * hfres.ints.sys.atoms.len) |coord| {
        for (0..H_d.nrow()) |mu| for (0..H_d.ncol()) |nu| {
            H_d.ptr(mu, nu).* = ScalarDual(T).init(hfres.ints.H.?.at(mu, nu), hfres.ints.dH.?.at(.{ coord, mu, nu }));
        };

        const nao = g_d.shape[0];

        for (0..nao) |mu| for (0..nao) |lambda| for (0..nao) |nu| for (0..nao) |sigma| {
            const indices = .{ mu, lambda, nu, sigma };

            const val = hfres.ints.g.?.at(indices);

            g_d.ptr(indices).* = ScalarDual(T).init(val, hfres.ints.dg.?.at(.{ coord, mu, lambda, nu, sigma }));
        };

        for (0..C_d.nrow()) |mu| for (0..C_d.ncol()) |p| {
            C_d.ptr(mu, p).* = ScalarDual(T).init(hfres.C_mo.at(mu, p), hfres.dC_mo.?.at(.{ coord, mu, p }));
        };

        var h_ms_d, var g_ms_d = try transformInts(ScalarDual(T), C_d, H_d, g_d, generalized, gpa);
        defer {
            h_ms_d.deinit(gpa);
            g_ms_d.deinit(gpa);
        }

        for (0..dets.len) |m| for (m..dets.len) |n| {
            const factor: T = if (m == n) 1 else 2;

            const val = slater(ScalarDual(T), dets[m], dets[n], h_ms_d, g_ms_d);

            for (states, 0..) |k, k_idx| {
                grads[k_idx].ptr(coord / 3, coord % 3).* += factor * C_ci.at(m, k) * C_ci.at(n, k) * val.der;
            }
        };
    }

    return grads;
}

/// Computes analytic nuclear gradients of CI state energies using density fitting three-center tensors.
pub fn calculateRiCiGradient(comptime T: type, hfres: HartreeFockResult(T), C_ci: Matrix(T), dets: []const []const usize, states: []const u32, gpa: Allocator) ![]Matrix(T) {
    const generalized = hfres.ints.sys.nbf != hfres.C_mo.shape[0];

    const B = hfres.ints.B.?;

    const J = hfres.ints.J orelse unreachable;

    const dI = hfres.ints.dI orelse unreachable;
    const dJ = hfres.ints.dJ orelse unreachable;

    const I_3c = hfres.ints.I orelse unreachable;

    const nbf, const naux = .{ B.shape[0], B.shape[2] };

    var grads = try gpa.alloc(Matrix(T), states.len);
    errdefer gpa.free(grads);

    var grads_done: usize = 0;
    errdefer for (0..grads_done) |k| grads[k].deinit(gpa);

    for (0..states.len) |idx| {
        grads[idx] = try hfres.ints.sys.nuclearRepulsionGradient(gpa);

        grads_done += 1;
    }

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

    var H_d = try Matrix(ScalarDual(T)).init(hfres.ints.H.?.nrow(), hfres.ints.H.?.ncol(), gpa);
    defer H_d.deinit(gpa);

    var C_d = try Matrix(ScalarDual(T)).init(hfres.C_mo.nrow(), hfres.C_mo.ncol(), gpa);
    defer C_d.deinit(gpa);

    for (0..3 * hfres.ints.sys.atoms.len) |coord| {
        try evaluateMetricDeriv(T, &dJ_inv_half, dJ.slice(coord).asMatrix(), decomp.U, M_kernel, gpa);

        evaluateFittingTensorDeriv(T, &dB_ao, dI.slice(coord), I_3c, decomp.J_inv_half, dJ_inv_half);

        for (0..nbf) |mu| for (0..nbf) |nu| for (0..naux) |P| {
            B_d.ptr(.{ mu, nu, P }).* = ScalarDual(T).init(B.at(.{ mu, nu, P }), dB_ao.at(.{ mu, nu, P }));
        };

        for (0..H_d.nrow()) |mu| for (0..H_d.ncol()) |nu| {
            H_d.ptr(mu, nu).* = ScalarDual(T).init(hfres.ints.H.?.at(mu, nu), hfres.ints.dH.?.at(.{ coord, mu, nu }));
        };

        for (0..C_d.nrow()) |mu| for (0..C_d.ncol()) |p| {
            C_d.ptr(mu, p).* = ScalarDual(T).init(hfres.C_mo.at(mu, p), hfres.dC_mo.?.at(.{ coord, mu, p }));
        };

        var h_ms_d, var g_ms_d = try transformIntsRi(ScalarDual(T), C_d, H_d, B_d, generalized, gpa);
        defer {
            h_ms_d.deinit(gpa);
            g_ms_d.deinit(gpa);
        }

        for (0..dets.len) |m| for (m..dets.len) |n| {
            const factor: T = if (m == n) 1 else 2;

            const val = slater(ScalarDual(T), dets[m], dets[n], h_ms_d, g_ms_d);

            for (states, 0..) |k, k_idx| {
                grads[k_idx].ptr(coord / 3, coord % 3).* += factor * C_ci.at(m, k) * C_ci.at(n, k) * val.der;
            }
        };
    }

    return grads;
}
