const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;
const Value = @import("value.zig").Value;

pub fn ao2mo_pppp(comptime T: type, g_pppp: *Tensor(T, 4), g_xxxx: Tensor(T, 4), C: Matrix(T), gpa: Allocator) !void {
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

    g_pppp.zero();

    const addTo = struct {
        pub fn addTo(ptr: *T, c: T, val: T) void {
            const val_c = Value(T).init(c);
            const val_val = Value(T).init(val);
            ptr.* = Value(T).init(ptr.*).add(val_c.mul(val_val)).val;
        }
    }.addTo;

    for (0..N) |i| for (0..N) |lambda| for (0..N) |nu| for (0..N) |mu| for (0..N) |sigma| {
        addTo(g_pxxx.ptr(.{ i, lambda, nu, sigma }), C.at(mu, i), g_xxxx.at(.{ mu, lambda, nu, sigma }));
    };

    for (0..N) |i| for (0..N) |j| for (0..N) |nu| for (0..N) |lambda| for (0..N) |sigma| {
        addTo(g_ppxx.ptr(.{ i, j, nu, sigma }), C.at(lambda, j), g_pxxx.at(.{ i, lambda, nu, sigma }));
    };

    for (0..N) |i| for (0..N) |j| for (0..N) |a| for (0..N) |nu| for (0..N) |sigma| {
        addTo(g_pppx.ptr(.{ i, j, a, sigma }), C.at(nu, a), g_ppxx.at(.{ i, j, nu, sigma }));
    };

    for (0..N) |i| for (0..N) |j| for (0..N) |a| for (0..N) |sigma| for (0..N) |b| {
        addTo(g_pppp.ptr(.{ i, j, a, b }), C.at(sigma, b), g_pppx.at(.{ i, j, a, sigma }));
    };
}

pub fn ao2mo_oovv(comptime T: type, g_oovv: *Tensor(T, 4), g_xxxx: Tensor(T, 4), C: Matrix(T), nocc: usize, gpa: Allocator) !void {
    const N = g_xxxx.shape[0];

    std.debug.assert(g_xxxx.shape[1] == N);
    std.debug.assert(g_xxxx.shape[2] == N);
    std.debug.assert(g_xxxx.shape[3] == N);

    std.debug.assert(C.shape[0] == N);
    std.debug.assert(C.shape[1] == N);

    std.debug.assert(nocc <= N);

    const nvir = N - nocc;

    std.debug.assert(g_oovv.shape[0] == nocc);
    std.debug.assert(g_oovv.shape[1] == nocc);
    std.debug.assert(g_oovv.shape[2] == nvir);
    std.debug.assert(g_oovv.shape[3] == nvir);

    var g_oxxx = try Tensor(T, 4).initZero(.{ nocc, N, N, N }, gpa);
    defer g_oxxx.deinit(gpa);

    var g_ooxx = try Tensor(T, 4).initZero(.{ nocc, nocc, N, N }, gpa);
    defer g_ooxx.deinit(gpa);

    var g_oovx = try Tensor(T, 4).initZero(.{ nocc, nocc, nvir, N }, gpa);
    defer g_oovx.deinit(gpa);

    g_oovv.zero();

    const addTo = struct {
        pub fn addTo(ptr: *T, c: T, g: T) void {
            const val_c = Value(T).init(c);
            const val_g = Value(T).init(g);

            ptr.* = Value(T).init(ptr.*).add(val_c.mul(val_g)).val;
        }
    }.addTo;

    for (0..nocc) |i| for (0..N) |lambda| for (0..N) |nu| for (0..N) |mu| for (0..N) |sigma| {
        addTo(g_oxxx.ptr(.{ i, lambda, nu, sigma }), C.at(mu, i), g_xxxx.at(.{ mu, lambda, nu, sigma }));
    };

    for (0..nocc) |i| for (0..nocc) |j| for (0..N) |nu| for (0..N) |lambda| for (0..N) |sigma| {
        addTo(g_ooxx.ptr(.{ i, j, nu, sigma }), C.at(lambda, j), g_oxxx.at(.{ i, lambda, nu, sigma }));
    };

    for (0..nocc) |i| for (0..nocc) |j| for (nocc..N) |a| for (0..N) |nu| for (0..N) |sigma| {
        addTo(g_oovx.ptr(.{ i, j, a - nocc, sigma }), C.at(nu, a), g_ooxx.at(.{ i, j, nu, sigma }));
    };

    for (0..nocc) |i| for (0..nocc) |j| for (0..nvir) |a| for (0..N) |sigma| for (nocc..N) |b| {
        addTo(g_oovv.ptr(.{ i, j, a, b - nocc }), C.at(sigma, b), g_oovx.at(.{ i, j, a, sigma }));
    };
}

pub fn ao2mo_pp(comptime T: type, A_pp: *Matrix(T), A_xx: Matrix(T), C: Matrix(T)) void {
    const N = C.shape[0];

    std.debug.assert(A_pp.nrow() == N);
    std.debug.assert(A_pp.ncol() == N);
    std.debug.assert(A_xx.nrow() == N);
    std.debug.assert(A_xx.ncol() == N);

    for (0..N) |p| for (0..N) |q| {
        var sum: T = 0;

        for (0..N) |mu| for (0..N) |nu| {
            sum += C.at(mu, p) * C.at(nu, q) * A_xx.at(mu, nu);
        };

        A_pp.ptr(p, q).* = sum;
    };
}

pub fn mo2ao_xx(comptime T: type, A_xx: *Matrix(T), A_pp: Matrix(T), C: Matrix(T)) void {
    const N = C.shape[0];

    std.debug.assert(A_pp.nrow() == N);
    std.debug.assert(A_pp.ncol() == N);
    std.debug.assert(A_xx.nrow() == N);
    std.debug.assert(A_xx.ncol() == N);

    for (0..N) |mu| for (0..N) |p| {
        var sum: T = 0;

        for (0..N) |q| {
            sum += C.at(mu, q) * A_pp.at(q, p);
        }

        A_xx.ptr(mu, p).* = sum;
    };
}
