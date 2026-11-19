//! Implements transformations of one- and two-electron integrals between atomic, spin, and molecular orbital bases.

const std = @import("std");

const Allocator = std.mem.Allocator;

const DualMatrix = @import("dual.zig").DualMatrix;
const DualTensor = @import("dual.zig").DualTensor;
const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;

const contract = @import("tensor_contraction.zig").contract;
const gemm = @import("linear_algebra.zig").gemm;
const isDualTensor = @import("value.zig").isDualTensor;

/// Transforms a one-electron operator matrix from atomic orbital to molecular orbital basis using BLAS.
pub fn ao2mo_1e(comptime T: type, A_mo: anytype, A_ao: anytype, C: anytype, gpa: Allocator) !void {
    var W = try @TypeOf(A_ao).init(A_ao.nrow(), C.ncol(), gpa);
    defer W.deinit(gpa);

    gemm(T, &W, A_ao, C, 1, 0, false, false);

    gemm(T, A_mo, C, W, 1, 0, true, false);
}

/// Transforms two-electron repulsion integrals from atomic orbital to molecular orbital basis using BLAS.
pub fn ao2mo_2e(comptime T: type, g_mo: anytype, g_ao: anytype, C: anytype, gpa: Allocator) !void {
    const N = C.shape[0];

    std.debug.assert(g_ao.shape[1] == N);
    std.debug.assert(g_ao.shape[2] == N);
    std.debug.assert(g_ao.shape[3] == N);

    std.debug.assert(C.shape[0] == N);
    std.debug.assert(C.shape[1] == N);

    var g_1mo = try @TypeOf(g_ao).init(.{ N, N, N, N }, gpa);
    defer g_1mo.deinit(gpa);

    var g_2mo = try @TypeOf(g_ao).init(.{ N, N, N, N }, gpa);
    defer g_2mo.deinit(gpa);

    try contract(T, "mi,mlns->ilns", &g_1mo, C, g_ao, 1, 0, gpa);

    try contract(T, "ilns,lj->ijns", &g_2mo, g_1mo, C, 1, 0, gpa);
    try contract(T, "ijns,na->ijas", &g_1mo, g_2mo, C, 1, 0, gpa);

    try contract(T, "ijas,sb->ijab", g_mo, g_1mo, C, 1, 0, gpa);
}

/// Transforms atomic orbital integrals to a restricted molecular orbital subset of occupied and virtual orbitals.
pub fn ao2mo_oovv(comptime T: type, g_oovv: anytype, g_ao: anytype, C: anytype, nocc: usize, gpa: Allocator) !void {
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

    var C_occ = try @TypeOf(C).init(N, nocc, gpa);
    defer C_occ.deinit(gpa);

    var C_vir = try @TypeOf(C).init(N, nvir, gpa);
    defer C_vir.deinit(gpa);

    for (0..N) |mu| for (0..nocc) |i| {
        C_occ.set(mu, i, C.at(mu, i));
    };

    for (0..N) |mu| for (0..nvir) |a| {
        C_vir.set(mu, a, C.at(mu, nocc + a));
    };

    var g_1mo = try @TypeOf(g_ao).init(.{ nocc, N, N, N }, gpa);
    defer g_1mo.deinit(gpa);

    var g_2mo = try @TypeOf(g_ao).init(.{ nocc, nocc, N, N }, gpa);
    defer g_2mo.deinit(gpa);

    var g_3mo = try @TypeOf(g_ao).init(.{ nocc, nocc, nvir, N }, gpa);
    defer g_3mo.deinit(gpa);

    try contract(T, "mi,mlns->ilns", &g_1mo, C_occ, g_ao, 1, 0, gpa);

    try contract(T, "ilns,lj->ijns", &g_2mo, g_1mo, C_occ, 1, 0, gpa);
    try contract(T, "ijns,na->ijas", &g_3mo, g_2mo, C_vir, 1, 0, gpa);
    try contract(T, "ijas,sb->ijab", g_oovv, g_3mo, C_vir, 1, 0, gpa);
}

/// Transforms three-center density fitting tensors from the atomic orbital basis to occupied-virtual molecular orbitals.
pub fn ao2mo_ov_3c(comptime T: type, B_mo: anytype, B_ao: anytype, C: anytype, gpa: Allocator) !void {
    const nbf, const naux = .{ B_ao.shape[0], B_ao.shape[2] };

    const nblocks: usize = if (C.shape[0] != nbf) 2 else 1;

    const nocc, const nvir = .{ B_mo.shape[0], B_mo.shape[1] };

    var C_occ = try @TypeOf(C).init(nbf, nocc, gpa);
    defer C_occ.deinit(gpa);

    var C_vir = try @TypeOf(C).init(nbf, nvir, gpa);
    defer C_vir.deinit(gpa);

    var B1 = try @TypeOf(B_ao).initZero(.{ nocc, nbf, naux }, gpa);
    defer B1.deinit(gpa);

    const B_mat, var B1_mat = .{ B_ao.matrixView(nbf, nbf * naux), B1.matrixView(nocc, nbf * naux) };

    B_mo.zero();

    for (0..nblocks) |b| {
        const row_off = b * nbf;

        for (0..nbf) |mu| {
            for (0..nocc) |i| C_occ.set(mu, i, C.at(row_off + mu, i));

            for (0..nvir) |a| C_vir.set(mu, a, C.at(row_off + mu, nocc + a));
        }

        gemm(T, &B1_mat, C_occ, B_mat, 1, 0, true, false);

        for (0..nocc) |i| {
            var B_mo_i = B_mo.slice(i).asMatrix();

            gemm(T, &B_mo_i, C_vir, B1.slice(i).asMatrix(), 1, 1, true, false);
        }
    }
}

/// Transforms three-center density fitting tensors from the atomic orbital basis to the molecular orbital basis using BLAS.
pub fn ao2mo_3c(comptime T: type, B_mo: anytype, B_ao: anytype, C: anytype, gpa: Allocator) !void {
    const nbf, const naux, const nmo = .{ B_ao.shape[0], B_ao.shape[2], C.shape[1] };

    var B1 = try @TypeOf(B_ao).initZero(.{ nmo, nbf, naux }, gpa);
    defer B1.deinit(gpa);

    const B_mat, var B1_mat = .{ B_ao.matrixView(nbf, nbf * naux), B1.matrixView(nmo, nbf * naux) };

    gemm(T, &B1_mat, C, B_mat, 1, 0, true, false);

    for (0..nmo) |p| {
        var B_mo_p = B_mo.slice(p).asMatrix();

        gemm(T, &B_mo_p, C, B1.slice(p).asMatrix(), 1, 0, true, false);
    }
}

/// Reconstructs spin-orbital two-electron repulsion integrals from 3-center fitting tensors using BLAS GEMM.
pub fn b2ms(comptime T: type, g_ms: anytype, B_mo: anytype, gpa: Allocator) !void {
    const is_dual = comptime isDualTensor(@TypeOf(B_mo));

    const MatType = if (is_dual) DualMatrix(T) else Matrix(T);

    const n, const naux = .{ g_ms.shape[0] / 2, B_mo.shape[2] };

    var M = try MatType.init(n * n, naux, gpa);
    defer M.deinit(gpa);

    if (comptime is_dual) {
        copyOrbitalBlock(M.val.data, B_mo.val.data, n, B_mo.shape[1], naux);
        copyOrbitalBlock(M.der.data, B_mo.der.data, n, B_mo.shape[1], naux);
    }

    if (comptime !is_dual) {
        copyOrbitalBlock(M.data, B_mo.data, n, B_mo.shape[1], naux);
    }

    var G = try MatType.init(n * n, n * n, gpa);
    defer G.deinit(gpa);

    gemm(T, &G, M, M, 1, 0, false, true);

    scatterSpatialToSpin(g_ms, G, n);
}

/// Computes atomic orbital response matrix from molecular orbital rotations and coefficients using BLAS.
pub fn mo2ao_1e(comptime T: type, A_ao: *Matrix(T), A_mo: Matrix(T), C: Matrix(T)) !void {
    const N = C.shape[0];

    std.debug.assert(A_mo.nrow() == N);
    std.debug.assert(A_mo.ncol() == N);
    std.debug.assert(A_ao.nrow() == N);
    std.debug.assert(A_ao.ncol() == N);

    gemm(T, A_ao, C, A_mo, 1, 0, false, false);
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

            C_ms.set(indmu0, col_a, C_mo.at(mu, p));
            C_ms.set(indmu1, col_b, C_mo.at(mu, p));
        }
    }
}

/// Transforms core Hamiltonian and two-electron repulsion or fitting tensors to molecular spin-orbital basis.
pub fn transformInts(comptime T: type, C: anytype, H: anytype, two_e: anytype, generalized: bool, gpa: Allocator) !struct { @TypeOf(H), if (isDualTensor(@TypeOf(two_e))) DualTensor(T, 4) else Tensor(T, 4) } {
    if (two_e.shape.len == 3) {
        return transformIntsRi(T, C, H, two_e, generalized, gpa);
    }

    const MatType, const TenType = .{ @TypeOf(H), @TypeOf(two_e) };

    if (!generalized) {
        const nmo, const nsp = .{ C.shape[1], 2 * C.shape[1] };

        var h_mo = try MatType.init(nmo, nmo, gpa);
        defer h_mo.deinit(gpa);

        try ao2mo_1e(T, &h_mo, H, C, gpa);

        var h_ms = try MatType.initZero(nsp, nsp, gpa);
        errdefer h_ms.deinit(gpa);

        for (0..nmo) |p| for (0..nmo) |q| {
            h_ms.set(2 * p + 0, 2 * q + 0, h_mo.at(p, q));
            h_ms.set(2 * p + 1, 2 * q + 1, h_mo.at(p, q));
        };

        var g_mo = try TenType.init(.{ nmo, nmo, nmo, nmo }, gpa);
        defer g_mo.deinit(gpa);

        try ao2mo_2e(T, &g_mo, two_e, C, gpa);

        var g_ms = try TenType.initZero(.{ nsp, nsp, nsp, nsp }, gpa);
        errdefer g_ms.deinit(gpa);

        for (0..nmo) |p| for (0..nmo) |r| {
            for (0..nmo) |q| for (0..nmo) |s| {
                const val = g_mo.at(.{ p, q, r, s });

                g_ms.set(.{ 2 * p + 0, 2 * q + 0, 2 * r + 0, 2 * s + 0 }, val);
                g_ms.set(.{ 2 * p + 0, 2 * q + 1, 2 * r + 0, 2 * s + 1 }, val);
                g_ms.set(.{ 2 * p + 1, 2 * q + 0, 2 * r + 1, 2 * s + 0 }, val);
                g_ms.set(.{ 2 * p + 1, 2 * q + 1, 2 * r + 1, 2 * s + 1 }, val);
            };
        };

        return .{ h_ms, g_ms };
    }

    const nsp = C.shape[1];

    var h_ms = try MatType.init(nsp, nsp, gpa);
    errdefer h_ms.deinit(gpa);

    try ao2mo_1e(T, &h_ms, H, C, gpa);

    var g_ms = try TenType.init(.{ nsp, nsp, nsp, nsp }, gpa);
    errdefer g_ms.deinit(gpa);

    try ao2mo_2e(T, &g_ms, two_e, C, gpa);

    return .{ h_ms, g_ms };
}

/// Transforms core Hamiltonian and density fitting 3-center tensors to molecular spin-orbital basis using BLAS.
pub fn transformIntsRi(comptime T: type, C: anytype, H: anytype, B_ao: anytype, generalized: bool, gpa: Allocator) !struct { @TypeOf(H), if (isDualTensor(@TypeOf(B_ao))) DualTensor(T, 4) else Tensor(T, 4) } {
    const Ten4Type = if (comptime isDualTensor(@TypeOf(B_ao))) DualTensor(T, 4) else Tensor(T, 4);

    const MatType, const naux = .{ @TypeOf(H), B_ao.shape[2] };

    const Ten3Type = @TypeOf(B_ao);

    if (!generalized) {
        const nmo, const nsp = .{ C.shape[1], 2 * C.shape[1] };

        var h_mo = try MatType.init(nmo, nmo, gpa);
        defer h_mo.deinit(gpa);

        try ao2mo_1e(T, &h_mo, H, C, gpa);

        var h_ms = try MatType.initZero(nsp, nsp, gpa);
        errdefer h_ms.deinit(gpa);

        for (0..nmo) |p| for (0..nmo) |q| {
            h_ms.set(2 * p + 0, 2 * q + 0, h_mo.at(p, q));
            h_ms.set(2 * p + 1, 2 * q + 1, h_mo.at(p, q));
        };

        var B_mo = try Ten3Type.initZero(.{ nmo, nmo, naux }, gpa);
        defer B_mo.deinit(gpa);

        try ao2mo_3c(T, &B_mo, B_ao, C, gpa);

        var g_ms = try Ten4Type.initZero(.{ nsp, nsp, nsp, nsp }, gpa);
        errdefer g_ms.deinit(gpa);

        try b2ms(T, &g_ms, B_mo, gpa);

        return .{ h_ms, g_ms };
    }

    const nbf, const nsp = .{ B_ao.shape[0], C.shape[1] };

    var h_ms = try MatType.init(nsp, nsp, gpa);
    errdefer h_ms.deinit(gpa);

    try ao2mo_1e(T, &h_ms, H, C, gpa);

    var C_a = try MatType.init(nbf, nsp, gpa);
    defer C_a.deinit(gpa);

    var C_b = try MatType.init(nbf, nsp, gpa);
    defer C_b.deinit(gpa);

    for (0..nbf) |mu| for (0..nsp) |p| {
        C_a.set(mu, p, C.at(mu, p));

        C_b.set(mu, p, C.at(nbf + mu, p));
    };

    var B_a = try Ten3Type.initZero(.{ nsp, nsp, naux }, gpa);
    defer B_a.deinit(gpa);

    try ao2mo_3c(T, &B_a, B_ao, C_a, gpa);

    var B_b = try Ten3Type.initZero(.{ nsp, nsp, naux }, gpa);
    defer B_b.deinit(gpa);

    try ao2mo_3c(T, &B_b, B_ao, C_b, gpa);

    var g_ms = try Ten4Type.initZero(.{ nsp, nsp, nsp, nsp }, gpa);
    errdefer g_ms.deinit(gpa);

    var B_tot = try Ten3Type.init(.{ nsp, nsp, naux }, gpa);
    defer B_tot.deinit(gpa);

    if (comptime isDualTensor(Ten3Type)) {
        for (0..B_tot.val.data.len) |i| {
            B_tot.val.data[i] = B_a.val.data[i] + B_b.val.data[i];
            B_tot.der.data[i] = B_a.der.data[i] + B_b.der.data[i];
        }
    }

    if (comptime !isDualTensor(Ten3Type)) {
        for (0..B_tot.data.len) |i| {
            B_tot.data[i] = B_a.data[i] + B_b.data[i];
        }
    }

    const M = B_tot.matrixView(nsp * nsp, naux);

    var G = try MatType.init(nsp * nsp, nsp * nsp, gpa);
    defer G.deinit(gpa);

    gemm(T, &G, M, M, 1, 0, false, true);

    for (0..nsp) |p| for (0..nsp) |r| {
        const pr = p * nsp + r;

        for (0..nsp) |q| for (0..nsp) |s| {
            const qs = q * nsp + s;

            g_ms.set(.{ p, q, r, s }, G.at(pr, qs));
        };
    };

    return .{ h_ms, g_ms };
}

/// Extracts an n-by-n spatial orbital block from a 3-center tensor into a flattened matrix buffer.
fn copyOrbitalBlock(dst: anytype, src: anytype, n: usize, stride: usize, naux: usize) void {
    if (stride == n) {
        @memcpy(dst[0 .. n * n * naux], src[0 .. n * n * naux]);

        return;
    }

    for (0..n) |p| for (0..n) |r| {
        const s, const d = .{ (p * stride + r) * naux, (p * n + r) * naux };

        @memcpy(dst[d .. d + naux], src[s .. s + naux]);
    };
}

/// Scatters spatial two-electron integrals to molecular spin-orbital block components.
fn scatterSpatialToSpin(g_ms: anytype, G: anytype, n: usize) void {
    for (0..n) |p| for (0..n) |r| {
        const pr = p * n + r;

        for (0..n) |q| for (0..n) |s| {
            const val = G.at(pr, q * n + s);

            g_ms.set(.{ 2 * p + 0, 2 * q + 0, 2 * r + 0, 2 * s + 0 }, val);
            g_ms.set(.{ 2 * p + 0, 2 * q + 1, 2 * r + 0, 2 * s + 1 }, val);
            g_ms.set(.{ 2 * p + 1, 2 * q + 0, 2 * r + 1, 2 * s + 0 }, val);
            g_ms.set(.{ 2 * p + 1, 2 * q + 1, 2 * r + 1, 2 * s + 1 }, val);
        };
    };
}
