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

pub fn ao2mo_pppp(comptime T: type, g_xxxx: Tensor(T, 4), C: Matrix(T), gpa: Allocator) !Tensor(T, 4) {
    const N = g_xxxx.shape[0];

    std.debug.assert(g_xxxx.shape[1] == N);
    std.debug.assert(g_xxxx.shape[2] == N);
    std.debug.assert(g_xxxx.shape[3] == N);

    std.debug.assert(C.shape[0] == N);
    std.debug.assert(C.shape[1] == N);

    var g_pxxx = try Tensor(T, 4).initZero(.{ N, N, N, N }, gpa);
    defer g_pxxx.deinit(gpa);

    var g_ppxx = try Tensor(T, 4).initZero(.{ N, N, N, N }, gpa);
    defer g_ppxx.deinit(gpa);

    var g_pppx = try Tensor(T, 4).initZero(.{ N, N, N, N }, gpa);
    defer g_pppx.deinit(gpa);

    var g_pppp = try Tensor(T, 4).initZero(.{ N, N, N, N }, gpa);
    errdefer g_pppp.deinit(gpa);

    for (0..N) |p| for (0..N) |lambda| for (0..N) |nu| for (0..N) |sigma| {
        var sum: T = 0;

        for (0..N) |mu| {
            sum += C.at(mu, p) * g_xxxx.at(.{ mu, lambda, nu, sigma });
        }

        g_pxxx.ptr(.{ p, lambda, nu, sigma }).* = sum;
    };

    for (0..N) |p| for (0..N) |q| for (0..N) |nu| for (0..N) |sigma| {
        var sum: T = 0;

        for (0..N) |lambda| {
            sum += C.at(lambda, q) * g_pxxx.at(.{ p, lambda, nu, sigma });
        }

        g_ppxx.ptr(.{ p, q, nu, sigma }).* = sum;
    };

    for (0..N) |p| for (0..N) |q| for (0..N) |r| for (0..N) |sigma| {
        var sum: T = 0;

        for (0..N) |nu| {
            sum += C.at(nu, r) * g_ppxx.at(.{ p, q, nu, sigma });
        }

        g_pppx.ptr(.{ p, q, r, sigma }).* = sum;
    };

    for (0..N) |p| for (0..N) |q| for (0..N) |r| for (0..N) |s| {
        var sum: T = 0;

        for (0..N) |sigma| {
            sum += C.at(sigma, s) * g_pppx.at(.{ p, q, r, sigma });
        }

        g_pppp.ptr(.{ p, q, r, s }).* = sum;
    };

    return g_pppp;
}

pub fn ao2mo_pp(comptime T: type, A_mo: *Matrix(T), A: Matrix(T), C: Matrix(T)) void {
    ao2mo_pp_general(T, A_mo, A, C, C);
}

pub fn ao2mo_pp_general(comptime T: type, A_mo: *Matrix(T), A: Matrix(T), C1: Matrix(T), C2: Matrix(T)) void {
    const N = A.nrow();

    std.debug.assert(C1.shape[0] == N);
    std.debug.assert(C1.shape[1] == N);
    std.debug.assert(C2.shape[0] == N);
    std.debug.assert(C2.shape[1] == N);
    std.debug.assert(A_mo.nrow() == N);
    std.debug.assert(A_mo.ncol() == N);

    for (0..N) |p| for (0..N) |q| {
        var sum: T = 0;

        for (0..N) |mu| for (0..N) |nu| {
            sum += C1.at(mu, p) * C2.at(nu, q) * A.at(mu, nu);
        };

        A_mo.ptr(p, q).* = sum;
    };
}

pub fn mo2ao_xx(comptime T: type, A_ao: *Matrix(T), A_mo: Matrix(T), C: Matrix(T)) void {
    const N = C.shape[0];

    std.debug.assert(A_mo.nrow() == N);
    std.debug.assert(A_mo.ncol() == N);
    std.debug.assert(A_ao.nrow() == N);
    std.debug.assert(A_ao.ncol() == N);

    for (0..N) |mu| for (0..N) |p| {
        var sum: T = 0;

        for (0..N) |q| {
            sum += C.at(mu, q) * A_mo.at(q, p);
        }

        A_ao.ptr(mu, p).* = sum;
    };
}
