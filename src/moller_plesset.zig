const std = @import("std");

const Allocator = std.mem.Allocator;

const HartreeFockOptions = @import("hartree_fock.zig").Options;
const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

const ao2mo_oovv = @import("integral_transform.zig").ao2mo_oovv;
const ao2mo_pp_general = @import("integral_transform.zig").ao2mo_pp_general;
const ao2mo_pppp = @import("integral_transform.zig").ao2mo_pppp;
const gradientCoef = @import("hartree_fock.zig").gradientCoef;
const gradientOrben = @import("hartree_fock.zig").gradientOrben;
const hartree_fock_run = @import("hartree_fock.zig").run;
const printf = @import("read_write.zig").printf;

// OPTIONS =============================================================================================================

pub const Options = struct {
    hartree_fock: HartreeFockOptions,

    order: u32 = 2,
    gradient: bool = false,
};

// MOLLER-PLESSET FUNCTIONS=============================================================================================

pub fn mp2(comptime T: type, hfres: HartreeFockResult(T), generalized: bool, gpa: Allocator) !T {
    const nocc = if (generalized) hfres.ints.sys.nel else hfres.ints.sys.nel / 2;

    var g_oovv = try ao2mo_oovv(T, hfres.ints.g.?, hfres.C, nocc, gpa);
    defer g_oovv.deinit(gpa);

    var energy: T = 0;

    for (0..nocc) |i| for (0..nocc) |j| for (0..hfres.C.shape[1] - nocc) |a| for (0..hfres.C.shape[1] - nocc) |b| {
        const denom = hfres.e.at(i) + hfres.e.at(j) - hfres.e.at(a + nocc) - hfres.e.at(b + nocc);

        if (generalized) {
            const term = g_oovv.at(.{ i, j, a, b }) - g_oovv.at(.{ i, j, b, a });

            energy += (term * term) / denom;
        }

        if (!generalized) {
            const term1 = g_oovv.at(.{ i, j, a, b });
            const term2 = g_oovv.at(.{ i, j, b, a });

            energy += (term1 * (2 * term1 - term2)) / denom;
        }
    };

    return if (generalized) 0.25 * energy else energy;
}

pub fn gradientMp2(comptime T: type, hfres: HartreeFockResult(T), generalized: bool, gpa: Allocator) !Matrix(T) {
    const dg = hfres.ints.dg orelse return error.CoulombDerivativeMatrixNotCalculated;

    const S = hfres.ints.S orelse return error.OverlapMatrixNotCalculated;

    const nbf = hfres.C.shape[0];

    var G = try Matrix(T).initZero(hfres.ints.sys.atoms.len, 3, gpa);
    errdefer G.deinit(gpa);

    var dC = try gradientCoef(T, hfres.ints, hfres.C, hfres.P, hfres.e, generalized, gpa);
    defer dC.deinit(gpa);

    var de = try gradientOrben(T, hfres.ints, hfres.C, hfres.P, hfres.e, generalized, gpa);
    defer de.deinit(gpa);

    var g_mo = try ao2mo_pppp(T, hfres.ints.g.?, hfres.C, gpa);
    defer g_mo.deinit(gpa);

    const nocc = if (generalized) hfres.ints.sys.nel else hfres.ints.sys.nel / 2;

    for (0..3 * hfres.ints.sys.atoms.len) |x| {
        var U_x = try Matrix(T).initZero(nbf, nbf, gpa);
        defer U_x.deinit(gpa);

        const dC_x = Matrix(T){ .data = dC.data[x * nbf * nbf .. (x + 1) * nbf * nbf], .shape = .{ nbf, nbf } };

        ao2mo_pp_general(T, &U_x, S, hfres.C, dC_x);

        var dg_x = try Tensor(T, 4).init(.{ nbf, nbf, nbf, nbf }, gpa);
        defer dg_x.deinit(gpa);

        for (0..nbf) |mu| for (0..nbf) |lambda| for (0..nbf) |nu| for (0..nbf) |sigma| {
            dg_x.ptr(.{ mu, lambda, nu, sigma }).* = dg.at(.{ x, mu, lambda, nu, sigma });
        };

        var dg_x_mo = try ao2mo_pppp(T, dg_x, hfres.C, gpa);
        defer dg_x_mo.deinit(gpa);

        var dE_mp2: T = 0;

        for (0..nocc) |i| for (0..nocc) |j| for (0..nbf - nocc) |a| for (0..nbf - nocc) |b| {
            const a_mo = a + nocc;
            const b_mo = b + nocc;

            const denom = hfres.e.at(i) + hfres.e.at(j) - hfres.e.at(a_mo) - hfres.e.at(b_mo);

            var dg_mo_ijab: T = dg_x_mo.at(.{ i, j, a_mo, b_mo });
            var dg_mo_ijba: T = dg_x_mo.at(.{ i, j, b_mo, a_mo });

            const ddenom = de.at(x, i) + de.at(x, j) - de.at(x, a_mo) - de.at(x, b_mo);

            for (0..nbf) |t| {
                const U_x_ti = U_x.at(t, i);
                const U_x_tj = U_x.at(t, j);

                const U_x_ta = U_x.at(t, a_mo);
                const U_x_tb = U_x.at(t, b_mo);

                const g_mo_tjab = g_mo.at(.{ t, j, a_mo, b_mo });
                const g_mo_itab = g_mo.at(.{ i, t, a_mo, b_mo });
                const g_mo_tjba = g_mo.at(.{ t, j, b_mo, a_mo });
                const g_mo_itba = g_mo.at(.{ i, t, b_mo, a_mo });

                const g_mo_ijtb = g_mo.at(.{ i, j, t, b_mo });
                const g_mo_ijat = g_mo.at(.{ i, j, a_mo, t });
                const g_mo_ijta = g_mo.at(.{ i, j, t, a_mo });
                const g_mo_ijbt = g_mo.at(.{ i, j, b_mo, t });

                dg_mo_ijab += U_x_ti * g_mo_tjab + U_x_tj * g_mo_itab + U_x_ta * g_mo_ijtb + U_x_tb * g_mo_ijat;
                dg_mo_ijba += U_x_ti * g_mo_tjba + U_x_tj * g_mo_itba + U_x_tb * g_mo_ijta + U_x_ta * g_mo_ijbt;
            }

            if (generalized) {
                const term = g_mo.at(.{ i, j, a_mo, b_mo }) - g_mo.at(.{ i, j, b_mo, a_mo });

                dE_mp2 += (2 * term * (dg_mo_ijab - dg_mo_ijba)) / denom - (term * term * ddenom) / (denom * denom);
            }

            if (!generalized) {
                const term1 = g_mo.at(.{ i, j, a_mo, b_mo });
                const term2 = g_mo.at(.{ i, j, b_mo, a_mo });

                const val = term1 * (2 * term1 - term2);

                const dterm1 = dg_mo_ijab;
                const dterm2 = dg_mo_ijba;

                const dval = dterm1 * (2 * term1 - term2) + term1 * (2 * dterm1 - dterm2);

                dE_mp2 += dval / denom - (val * ddenom) / (denom * denom);
            }
        };

        G.ptr(x / 3, x % 3).* = if (generalized) 0.25 * dE_mp2 else dE_mp2;
    }

    return G;
}

// RESULT STRUCT =======================================================================================================

pub fn Result(comptime T: type) type {
    return struct {
        hartree_fock: HartreeFockResult(T),

        energy: []T,

        gradient: []Matrix(T),

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.hartree_fock.deinit(gpa);

            gpa.free(self.energy);

            for (self.gradient) |*G| {
                G.deinit(gpa);
            }

            gpa.free(self.gradient);
        }
    };
}

// RUN =================================================================================================================

pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator) !Result(T) {
    if (opt.order < 2) {
        return error.MollerPlessetOrderTooLow;
    }

    if (opt.order > 2) {
        return error.MollerPlessetOrderNotSupported;
    }

    var hf_opt = opt.hartree_fock;

    if (opt.gradient) {
        hf_opt.gradient = true;
    }

    var hfres = try hartree_fock_run(T, io, hf_opt, log, gpa);
    errdefer hfres.deinit(gpa);

    var energy = try gpa.alloc(T, 1);
    errdefer gpa.free(energy);

    const grad = try gpa.alloc(Matrix(T), if (opt.gradient) 1 else 0);
    errdefer if (opt.gradient) gpa.free(grad);

    const corr_energy = switch (opt.order) {
        2 => try mp2(T, hfres, opt.hartree_fock.generalized, gpa),
        else => unreachable,
    };

    energy[0] = hfres.energy[0] + corr_energy;

    if (log) try printf(io, "\nMP{d} CORRELATION ENERGY: {d:.14}\n", .{ opt.order, corr_energy });

    if (log) try printf(io, "\nFINAL MP{d} ENERGY: {d:.14}\n", .{ opt.order, energy[0] });

    if (opt.gradient) {
        var corr_grad = try gradientMp2(T, hfres, opt.hartree_fock.generalized, gpa);
        defer corr_grad.deinit(gpa);

        var total_grad = try Matrix(T).init(corr_grad.nrow(), corr_grad.ncol(), gpa);
        errdefer total_grad.deinit(gpa);

        for (0..total_grad.nrow()) |i| for (0..total_grad.ncol()) |j| {
            total_grad.ptr(i, j).* = hfres.gradient[0].at(i, j) + corr_grad.at(i, j);
        };

        grad[0] = total_grad;
    }

    if (log and opt.gradient) {
        try printf(io, "\nMP{d} NUCLEAR ENERGY GRADIENT\n", .{opt.order});

        for (0..grad[0].shape[0]) |i| for (0..grad[0].shape[1]) |j| {
            try printf(io, "{d:20.14}{s}", .{ grad[0].at(i, j), if (j == 2) "\n" else " " });
        };
    }

    return Result(T){ .hartree_fock = hfres, .energy = energy, .gradient = grad };
}
