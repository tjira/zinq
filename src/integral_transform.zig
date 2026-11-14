//! Implements transformations of one- and two-electron integrals between atomic, spin, and molecular orbital bases.

const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;
const Value = @import("value.zig").Value;

const contract = @import("tensor_contraction.zig").contract;
const gemm = @import("linear_algebra.zig").gemm;

/// Transforms a one-electron operator from the spatial atomic orbital basis to the atomic spin-orbital basis.
pub fn ao2as_1e(comptime T: type, A_as: *Matrix(T), A_ao: Matrix(T)) void {
    const nbf = A_ao.shape[0];

    std.debug.assert(A_as.shape[0] == 2 * nbf);
    std.debug.assert(A_as.shape[1] == 2 * nbf);

    A_as.zero();

    for (0..nbf) |mu| for (0..nbf) |nu| {
        A_as.ptr(mu, nu).* = A_ao.at(mu, nu);

        A_as.ptr(mu + nbf, nu + nbf).* = A_ao.at(mu, nu);
    };
}

/// Transforms two-electron repulsion integrals from the spatial atomic orbital basis to the atomic spin-orbital basis.
pub fn ao2as_2e(comptime T: type, g_as: *Tensor(T, 4), g_ao: Tensor(T, 4)) void {
    const nbf = g_ao.shape[0];

    std.debug.assert(g_ao.shape[1] == nbf);
    std.debug.assert(g_ao.shape[2] == nbf);
    std.debug.assert(g_ao.shape[3] == nbf);

    std.debug.assert(g_as.shape[0] == 2 * nbf);
    std.debug.assert(g_as.shape[1] == 2 * nbf);
    std.debug.assert(g_as.shape[2] == 2 * nbf);
    std.debug.assert(g_as.shape[3] == 2 * nbf);

    g_as.zero();

    for (0..nbf) |mu| for (0..nbf) |nu| for (0..nbf) |lambda| for (0..nbf) |sigma| {
        g_as.ptr(.{ mu, nu, lambda, sigma }).* = g_ao.at(.{ mu, nu, lambda, sigma });

        g_as.ptr(.{ mu + nbf, nu, lambda + nbf, sigma }).* = g_ao.at(.{ mu, nu, lambda, sigma });
        g_as.ptr(.{ mu, nu + nbf, lambda, sigma + nbf }).* = g_ao.at(.{ mu, nu, lambda, sigma });

        g_as.ptr(.{ mu + nbf, nu + nbf, lambda + nbf, sigma + nbf }).* = g_ao.at(.{ mu, nu, lambda, sigma });
    };
}

/// Transforms a one-electron operator matrix from the atomic orbital basis to the molecular orbital basis.
pub fn ao2mo_1e(comptime T: type, A_mo: *Matrix(T), A_ao: Matrix(T), C: Matrix(T)) void {
    const N = C.shape[0];

    std.debug.assert(A_mo.nrow() == N);
    std.debug.assert(A_mo.ncol() == N);
    std.debug.assert(A_ao.nrow() == N);
    std.debug.assert(A_ao.ncol() == N);

    for (0..N) |p| for (0..N) |q| {
        var sum = Value(T).fromFloat(0);

        for (0..N) |mu| for (0..N) |nu| {
            const c_mu_p = Value(T).init(C.at(mu, p));
            const c_nu_q = Value(T).init(C.at(nu, q));

            const a_val = Value(T).init(A_ao.at(mu, nu));

            sum = sum.add(c_mu_p.mul(c_nu_q).mul(a_val));
        };

        A_mo.ptr(p, q).* = sum.val;
    };
}

/// Transforms two-electron repulsion integrals from the atomic orbital basis to the molecular orbital basis.
pub fn ao2mo_2e(comptime T: type, g_mo: *Tensor(T, 4), g_ao: Tensor(T, 4), C: Matrix(T), gpa: Allocator) !void {
    const N = g_ao.shape[0];

    std.debug.assert(g_ao.shape[1] == N);
    std.debug.assert(g_ao.shape[2] == N);
    std.debug.assert(g_ao.shape[3] == N);

    std.debug.assert(C.shape[0] == N);
    std.debug.assert(C.shape[1] == N);

    if (T == f64) {
        var g_1mo = try Tensor(T, 4).init(.{ N, N, N, N }, gpa);
        defer g_1mo.deinit(gpa);

        var g_2mo = try Tensor(T, 4).init(.{ N, N, N, N }, gpa);
        defer g_2mo.deinit(gpa);

        var g_3mo = try Tensor(T, 4).init(.{ N, N, N, N }, gpa);
        defer g_3mo.deinit(gpa);

        try contract(f64, "mi,mlns->ilns", &g_1mo, C, g_ao, 1, 0, gpa);

        try contract(f64, "ilns,lj->ijns", &g_2mo, g_1mo, C, 1, 0, gpa);
        try contract(f64, "ijns,na->ijas", &g_3mo, g_2mo, C, 1, 0, gpa);

        try contract(f64, "ijas,sb->ijab", g_mo, g_3mo, C, 1, 0, gpa);

        return;
    }

    var g_1mo = try Tensor(T, 4).initZero(.{ N, N, N, N }, gpa);
    defer g_1mo.deinit(gpa);

    var g_2mo = try Tensor(T, 4).initZero(.{ N, N, N, N }, gpa);
    defer g_2mo.deinit(gpa);

    var g_3mo = try Tensor(T, 4).initZero(.{ N, N, N, N }, gpa);
    defer g_3mo.deinit(gpa);

    g_mo.zero();

    const addTo = struct {
        pub fn addTo(ptr: *T, c: T, val: T) void {
            ptr.* = Value(T).init(ptr.*).add(Value(T).init(c).mul(Value(T).init(val))).val;
        }
    }.addTo;

    for (0..N) |p| for (0..N) |lambda| for (0..N) |nu| for (0..N) |mu| for (0..N) |sigma| {
        addTo(g_1mo.ptr(.{ p, lambda, nu, sigma }), C.at(mu, p), g_ao.at(.{ mu, lambda, nu, sigma }));
    };

    for (0..N) |p| for (0..N) |q| for (0..N) |nu| for (0..N) |lambda| for (0..N) |sigma| {
        addTo(g_2mo.ptr(.{ p, q, nu, sigma }), C.at(lambda, q), g_1mo.at(.{ p, lambda, nu, sigma }));
    };

    for (0..N) |p| for (0..N) |q| for (0..N) |r| for (0..N) |nu| for (0..N) |sigma| {
        addTo(g_3mo.ptr(.{ p, q, r, sigma }), C.at(nu, r), g_2mo.at(.{ p, q, nu, sigma }));
    };

    for (0..N) |p| for (0..N) |q| for (0..N) |r| for (0..N) |sigma| for (0..N) |s| {
        addTo(g_mo.ptr(.{ p, q, r, s }), C.at(sigma, s), g_3mo.at(.{ p, q, r, sigma }));
    };
}

/// Transforms atomic orbital integrals to a restricted molecular orbital subset of occupied and virtual orbitals.
pub fn ao2mo_oovv(comptime T: type, g_oovv: *Tensor(T, 4), g_ao: Tensor(T, 4), C: Matrix(T), nocc: usize, gpa: Allocator) !void {
    const N = g_ao.shape[0];

    std.debug.assert(g_ao.shape[1] == N);
    std.debug.assert(g_ao.shape[2] == N);
    std.debug.assert(g_ao.shape[3] == N);

    std.debug.assert(C.shape[0] == N);
    std.debug.assert(C.shape[1] == N);

    std.debug.assert(nocc <= N);

    const nvir = N - nocc;

    std.debug.assert(g_oovv.shape[0] == nocc);
    std.debug.assert(g_oovv.shape[1] == nocc);
    std.debug.assert(g_oovv.shape[2] == nvir);
    std.debug.assert(g_oovv.shape[3] == nvir);

    if (T == f64) {
        var C_occ = try Matrix(T).init(N, nocc, gpa);
        defer C_occ.deinit(gpa);

        var C_vir = try Matrix(T).init(N, nvir, gpa);
        defer C_vir.deinit(gpa);

        for (0..N) |mu| for (0..nocc) |i| {
            C_occ.ptr(mu, i).* = C.at(mu, i);
        };

        for (0..N) |mu| for (0..nvir) |a| {
            C_vir.ptr(mu, a).* = C.at(mu, nocc + a);
        };

        var g_1mo = try Tensor(T, 4).init(.{ nocc, N, N, N }, gpa);
        defer g_1mo.deinit(gpa);

        var g_2mo = try Tensor(T, 4).init(.{ nocc, nocc, N, N }, gpa);
        defer g_2mo.deinit(gpa);

        var g_3mo = try Tensor(T, 4).init(.{ nocc, nocc, nvir, N }, gpa);
        defer g_3mo.deinit(gpa);

        try contract(f64, "mi,mlns->ilns", &g_1mo, C_occ, g_ao, 1, 0, gpa);

        try contract(f64, "ilns,lj->ijns", &g_2mo, g_1mo, C_occ, 1, 0, gpa);
        try contract(f64, "ijns,na->ijas", &g_3mo, g_2mo, C_vir, 1, 0, gpa);

        try contract(f64, "ijas,sb->ijab", g_oovv, g_3mo, C_vir, 1, 0, gpa);

        return;
    }

    var g_1mo = try Tensor(T, 4).initZero(.{ nocc, N, N, N }, gpa);
    defer g_1mo.deinit(gpa);

    var g_2mo = try Tensor(T, 4).initZero(.{ nocc, nocc, N, N }, gpa);
    defer g_2mo.deinit(gpa);

    var g_3mo = try Tensor(T, 4).initZero(.{ nocc, nocc, nvir, N }, gpa);
    defer g_3mo.deinit(gpa);

    g_oovv.zero();

    const addTo = struct {
        pub fn addTo(ptr: *T, c: T, g: T) void {
            const val_c = Value(T).init(c);
            const val_g = Value(T).init(g);

            ptr.* = Value(T).init(ptr.*).add(val_c.mul(val_g)).val;
        }
    }.addTo;

    for (0..nocc) |i| for (0..N) |lambda| for (0..N) |nu| for (0..N) |mu| for (0..N) |sigma| {
        addTo(g_1mo.ptr(.{ i, lambda, nu, sigma }), C.at(mu, i), g_ao.at(.{ mu, lambda, nu, sigma }));
    };

    for (0..nocc) |i| for (0..nocc) |j| for (0..N) |nu| for (0..N) |lambda| for (0..N) |sigma| {
        addTo(g_2mo.ptr(.{ i, j, nu, sigma }), C.at(lambda, j), g_1mo.at(.{ i, lambda, nu, sigma }));
    };

    for (0..nocc) |i| for (0..nocc) |j| for (nocc..N) |a| for (0..N) |nu| for (0..N) |sigma| {
        addTo(g_3mo.ptr(.{ i, j, a - nocc, sigma }), C.at(nu, a), g_2mo.at(.{ i, j, nu, sigma }));
    };

    for (0..nocc) |i| for (0..nocc) |j| for (0..nvir) |a| for (0..N) |sigma| for (nocc..N) |b| {
        addTo(g_oovv.ptr(.{ i, j, a, b - nocc }), C.at(sigma, b), g_3mo.at(.{ i, j, a, sigma }));
    };
}

/// Transforms three-center density fitting tensors from the atomic orbital basis to occupied-virtual molecular orbitals.
pub fn ao2mo_ov_3c(comptime T: type, B_mo: *Tensor(T, 3), B_ao: Tensor(T, 3), C: Matrix(T), gpa: Allocator) !void {
    const nocc, const nvir = .{ B_mo.shape[0], B_mo.shape[1] };

    const nbf, const naux = .{ B_ao.shape[0], B_ao.shape[2] };

    const nblocks: usize = if (C.shape[0] != nbf) 2 else 1;

    var C_occ = try Matrix(T).init(nbf, nocc, gpa);
    defer C_occ.deinit(gpa);

    var C_vir = try Matrix(T).init(nbf, nvir, gpa);
    defer C_vir.deinit(gpa);

    var B1 = try Tensor(T, 3).initZero(.{ nocc, nbf, naux }, gpa);
    defer B1.deinit(gpa);

    const B_mat = B_ao.matrixView(nbf, nbf * naux);

    var B1_mat = B1.matrixView(nocc, nbf * naux);

    B_mo.zero();

    for (0..nblocks) |b| {
        const row_off = b * nbf;

        for (0..nbf) |mu| {
            for (0..nocc) |i| {
                C_occ.ptr(mu, i).* = C.at(row_off + mu, i);
            }

            for (0..nvir) |a| {
                C_vir.ptr(mu, a).* = C.at(row_off + mu, nocc + a);
            }
        }

        if (comptime T == f64) {
            gemm(T, &B1_mat, C_occ, B_mat, 1, 0, true, false);

            for (0..nocc) |i| {
                const B1_i = B1.slice(i).asMatrix();

                var B_mo_i = B_mo.slice(i).asMatrix();

                gemm(T, &B_mo_i, C_vir, B1_i, 1, 1, true, false);
            }
        }

        if (comptime T != f64) {
            for (0..nocc) |i| for (0..nbf) |nu| for (0..naux) |P| {
                var sum = Value(T).fromFloat(0);

                for (0..nbf) |mu| {
                    sum = sum.add(Value(T).init(C_occ.at(mu, i)).mul(Value(T).init(B_ao.at(.{ mu, nu, P }))));
                }

                B1.ptr(.{ i, nu, P }).* = sum.val;
            };

            for (0..nocc) |i| for (0..nvir) |a| for (0..naux) |P| {
                var sum = Value(T).init(B_mo.at(.{ i, a, P }));

                for (0..nbf) |nu| {
                    sum = sum.add(Value(T).init(C_vir.at(nu, a)).mul(Value(T).init(B1.at(.{ i, nu, P }))));
                }

                B_mo.ptr(.{ i, a, P }).* = sum.val;
            };
        }
    }
}

/// Transforms three-center density fitting tensors from the atomic orbital basis to the molecular orbital basis.
pub fn ao2mo_3c(comptime T: type, B_mo: *Tensor(T, 3), B_ao: Tensor(T, 3), C: Matrix(T), gpa: Allocator) !void {
    const nbf, const naux, const nmo = .{ B_ao.shape[0], B_ao.shape[2], C.shape[1] };

    var B1 = try Tensor(T, 3).initZero(.{ nmo, nbf, naux }, gpa);
    defer B1.deinit(gpa);

    if (comptime T == f64) {
        const B_mat = B_ao.matrixView(nbf, nbf * naux);

        var B1_mat = B1.matrixView(nmo, nbf * naux);

        gemm(T, &B1_mat, C, B_mat, 1, 0, true, false);

        for (0..nmo) |p| {
            const B1_p = B1.slice(p).asMatrix();

            var B_mo_p = B_mo.slice(p).asMatrix();

            gemm(T, &B_mo_p, C, B1_p, 1, 0, true, false);
        }

        return;
    }

    for (0..nmo) |p| for (0..nbf) |nu| for (0..naux) |P| {
        var sum = Value(T).fromFloat(0);

        for (0..nbf) |mu| {
            sum = sum.add(Value(T).init(C.at(mu, p)).mul(Value(T).init(B_ao.at(.{ mu, nu, P }))));
        }

        B1.ptr(.{ p, nu, P }).* = sum.val;
    };

    for (0..nmo) |p| for (0..nmo) |q| for (0..naux) |P| {
        var sum = Value(T).fromFloat(0);

        for (0..nbf) |nu| {
            sum = sum.add(Value(T).init(C.at(nu, q)).mul(Value(T).init(B1.at(.{ p, nu, P }))));
        }

        B_mo.ptr(.{ p, q, P }).* = sum.val;
    };
}

/// Transforms a one-electron operator matrix from the molecular orbital basis back to the atomic orbital basis.
pub fn mo2ao_1e(comptime T: type, A_ao: *Matrix(T), A_mo: Matrix(T), C: Matrix(T)) !void {
    const N = C.shape[0];

    std.debug.assert(A_mo.nrow() == N);
    std.debug.assert(A_mo.ncol() == N);
    std.debug.assert(A_ao.nrow() == N);
    std.debug.assert(A_ao.ncol() == N);

    if (T == f64) {
        return try contract(f64, "mq,qp->mp", A_ao, C, A_mo, 1, 0, null);
    }

    for (0..N) |mu| for (0..N) |p| {
        var sum = Value(T).fromFloat(0);

        for (0..N) |q| {
            sum = sum.add(Value(T).init(C.at(mu, q)).mul(Value(T).init(A_mo.at(q, p))));
        }

        A_ao.ptr(mu, p).* = sum.val;
    };
}

/// Expands spatial molecular orbital expansion coefficients into a molecular spin-orbital block matrix.
pub fn mo2ms_coef(comptime T: type, C_ms: *Matrix(T), C_mo: Matrix(T)) void {
    const nbf = C_mo.shape[0];

    std.debug.assert(C_ms.shape[0] == 2 * nbf);
    std.debug.assert(C_ms.shape[1] == 2 * nbf);

    C_ms.zero();

    for (0..nbf) |p| {
        const col_a = 2 * p + 0;
        const col_b = 2 * p + 1;

        for (0..nbf) |mu| {
            const indmu0, const indmu1 = .{ mu, mu + nbf };

            C_ms.ptr(indmu0, col_a).* = C_mo.at(mu, p);
            C_ms.ptr(indmu1, col_b).* = C_mo.at(mu, p);
        }
    }
}

/// Transforms core Hamiltonian and two-electron repulsion integrals to the molecular spin-orbital basis.
pub fn transformInts(comptime T: type, C: Matrix(T), H: Matrix(T), g: Tensor(T, 4), generalized: bool, gpa: Allocator) !struct { Matrix(T), Tensor(T, 4) } {
    const nsp = if (generalized) C.shape[0] else 2 * C.shape[0];

    var C_ms = if (generalized) C else try Matrix(T).initZero(nsp, nsp, gpa);
    defer if (!generalized) C_ms.deinit(gpa);

    if (!generalized) {
        mo2ms_coef(T, &C_ms, C);
    }

    var H_as = if (generalized) H else try Matrix(T).initZero(nsp, nsp, gpa);
    defer if (!generalized) H_as.deinit(gpa);

    if (!generalized) {
        ao2as_1e(T, &H_as, H);
    }

    var g_as = if (generalized) g else try Tensor(T, 4).initZero(.{ nsp, nsp, nsp, nsp }, gpa);
    defer if (!generalized) g_as.deinit(gpa);

    if (!generalized) {
        ao2as_2e(T, &g_as, g);
    }

    var h_ms = try Matrix(T).init(nsp, nsp, gpa);
    errdefer h_ms.deinit(gpa);

    ao2mo_1e(T, &h_ms, H_as, C_ms);

    var g_ms = try Tensor(T, 4).init(.{ nsp, nsp, nsp, nsp }, gpa);
    errdefer g_ms.deinit(gpa);

    try ao2mo_2e(T, &g_ms, g_as, C_ms, gpa);

    return .{ h_ms, g_ms };
}

/// Transforms core Hamiltonian and density fitting three-center tensors to the molecular spin-orbital basis.
pub fn transformIntsRi(comptime T: type, C: Matrix(T), H: Matrix(T), B_ao: Tensor(T, 3), generalized: bool, gpa: Allocator) !struct { Matrix(T), Tensor(T, 4) } {
    const naux = B_ao.shape[2];

    if (!generalized) {
        const nmo = C.shape[1];
        const nsp = 2 * nmo;

        var h_mo = try Matrix(T).init(nmo, nmo, gpa);
        defer h_mo.deinit(gpa);

        ao2mo_1e(T, &h_mo, H, C);

        var h_ms = try Matrix(T).initZero(nsp, nsp, gpa);
        errdefer h_ms.deinit(gpa);

        for (0..nmo) |p| for (0..nmo) |q| {
            h_ms.ptr(2 * p, 2 * q).* = h_mo.at(p, q);
            h_ms.ptr(2 * p + 1, 2 * q + 1).* = h_mo.at(p, q);
        };

        var B_mo = try Tensor(T, 3).initZero(.{ nmo, nmo, naux }, gpa);
        defer B_mo.deinit(gpa);

        try ao2mo_3c(T, &B_mo, B_ao, C, gpa);

        var g_ms = try Tensor(T, 4).initZero(.{ nsp, nsp, nsp, nsp }, gpa);
        errdefer g_ms.deinit(gpa);

        for (0..nmo) |p| for (0..nmo) |r| {
            for (0..nmo) |q| for (0..nmo) |s| {
                var sum = Value(T).fromFloat(0);

                for (0..naux) |P| {
                    const b_pr = Value(T).init(B_mo.at(.{ p, r, P }));
                    const b_qs = Value(T).init(B_mo.at(.{ q, s, P }));

                    sum = sum.add(b_pr.mul(b_qs));
                }

                g_ms.ptr(.{ 2 * p, 2 * q, 2 * r, 2 * s }).* = sum.val;
                g_ms.ptr(.{ 2 * p + 1, 2 * q + 1, 2 * r + 1, 2 * s + 1 }).* = sum.val;
                g_ms.ptr(.{ 2 * p, 2 * q + 1, 2 * r, 2 * s + 1 }).* = sum.val;
                g_ms.ptr(.{ 2 * p + 1, 2 * q, 2 * r + 1, 2 * s }).* = sum.val;
            };
        };

        return .{ h_ms, g_ms };
    }

    const nbf = B_ao.shape[0];
    const nsp = C.shape[1];

    var h_ms = try Matrix(T).init(nsp, nsp, gpa);
    errdefer h_ms.deinit(gpa);

    ao2mo_1e(T, &h_ms, H, C);

    var C_a = try Matrix(T).init(nbf, nsp, gpa);
    defer C_a.deinit(gpa);

    var C_b = try Matrix(T).init(nbf, nsp, gpa);
    defer C_b.deinit(gpa);

    for (0..nbf) |mu| for (0..nsp) |p| {
        C_a.ptr(mu, p).* = C.at(mu, p);
        C_b.ptr(mu, p).* = C.at(nbf + mu, p);
    };

    var B_a = try Tensor(T, 3).initZero(.{ nsp, nsp, naux }, gpa);
    defer B_a.deinit(gpa);

    try ao2mo_3c(T, &B_a, B_ao, C_a, gpa);

    var B_b = try Tensor(T, 3).initZero(.{ nsp, nsp, naux }, gpa);
    defer B_b.deinit(gpa);

    try ao2mo_3c(T, &B_b, B_ao, C_b, gpa);

    var g_ms = try Tensor(T, 4).initZero(.{ nsp, nsp, nsp, nsp }, gpa);
    errdefer g_ms.deinit(gpa);

    for (0..nsp) |p| for (0..nsp) |r| {
        for (0..nsp) |q| for (0..nsp) |s| {
            var sum = Value(T).fromFloat(0);

            for (0..naux) |P| {
                const b_pr = Value(T).init(B_a.at(.{ p, r, P })).add(Value(T).init(B_b.at(.{ p, r, P })));
                const b_qs = Value(T).init(B_a.at(.{ q, s, P })).add(Value(T).init(B_b.at(.{ q, s, P })));

                sum = sum.add(b_pr.mul(b_qs));
            }

            g_ms.ptr(.{ p, q, r, s }).* = sum.val;
        };
    };

    return .{ h_ms, g_ms };
}
