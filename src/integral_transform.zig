const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;

pub fn ao2mo_oovv(comptime T: type, g_xxxx: Tensor(T, 4), C: Matrix(T), nocc: usize, gpa: Allocator) !Tensor(T, 4) {
    const N = g_xxxx.shape[0];

    std.debug.assert(g_xxxx.shape[1] == N);
    std.debug.assert(g_xxxx.shape[2] == N);
    std.debug.assert(g_xxxx.shape[3] == N);

    std.debug.assert(C.shape[0] == N);
    std.debug.assert(C.shape[1] == N);

    std.debug.assert(nocc <= N);

    const nvir = N - nocc;

    var g_oxxx = try Tensor(T, 4).initZero(.{ nocc, N, N, N }, gpa);
    defer g_oxxx.deinit(gpa);

    var g_ooxx = try Tensor(T, 4).initZero(.{ nocc, nocc, N, N }, gpa);
    defer g_ooxx.deinit(gpa);

    var g_oovx = try Tensor(T, 4).initZero(.{ nocc, nocc, nvir, N }, gpa);
    defer g_oovx.deinit(gpa);

    var g_oovv = try Tensor(T, 4).initZero(.{ nocc, nocc, nvir, nvir }, gpa);
    errdefer g_oovv.deinit(gpa);

    for (0..nocc) |i| for (0..N) |lambda| for (0..N) |nu| for (0..N) |mu| for (0..N) |sigma| {
        g_oxxx.ptr(.{ i, lambda, nu, sigma }).* += C.at(mu, i) * g_xxxx.at(.{ mu, lambda, nu, sigma });
    };

    for (0..nocc) |i| for (0..nocc) |j| for (0..N) |nu| for (0..N) |lambda| for (0..N) |sigma| {
        g_ooxx.ptr(.{ i, j, nu, sigma }).* += C.at(lambda, j) * g_oxxx.at(.{ i, lambda, nu, sigma });
    };

    for (0..nocc) |i| for (0..nocc) |j| for (nocc..N) |a| for (0..N) |nu| for (0..N) |sigma| {
        g_oovx.ptr(.{ i, j, a - nocc, sigma }).* += C.at(nu, a) * g_ooxx.at(.{ i, j, nu, sigma });
    };

    for (0..nocc) |i| for (0..nocc) |j| for (0..nvir) |a| for (0..N) |sigma| for (nocc..N) |b| {
        g_oovv.ptr(.{ i, j, a, b - nocc }).* += C.at(sigma, b) * g_oovx.at(.{ i, j, a, sigma });
    };

    return g_oovv;
}
