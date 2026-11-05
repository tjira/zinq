//! Solves the Coupled-Perturbed Hartree-Fock equations to obtain the first-order orbital response to perturbations.

const std = @import("std");

const Allocator = std.mem.Allocator;

const Diis = @import("diis.zig").Diis;
const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

const ao2mo_1e = @import("integral_transform.zig").ao2mo_1e;
const mo2ao_1e = @import("integral_transform.zig").mo2ao_1e;
const printf = @import("read_write.zig").printf;

/// Computes the first-order derivative of molecular orbital coefficients and energies using CPHF.
pub fn orbitalResponse(comptime T: type, io: std.Io, hfres: HartreeFockResult(T), opt: anytype, log: bool, gpa: Allocator) !struct { Tensor(T, 3), Matrix(T) } {
    const nbf = hfres.C_mo.shape[0];

    var dC = try Tensor(T, 3).initZero(.{ 3 * hfres.ints.sys.atoms.len, nbf, nbf }, gpa);
    errdefer dC.deinit(gpa);

    var de = try Matrix(T).initZero(3 * hfres.ints.sys.atoms.len, nbf, gpa);
    errdefer de.deinit(gpa);

    try solveCPHF(T, io, &dC, &de, hfres, opt, log, gpa);

    return .{ dC, de };
}

/// Iteratively solves the CPHF linear equations to determine the first-order orbital response matrix.
pub fn solveCPHF(comptime T: type, io: std.Io, dC: *Tensor(T, 3), de: *Matrix(T), hfres: HartreeFockResult(T), opt: anytype, log: bool, gpa: Allocator) !void {
    const generalized = hfres.ints.sys.nbf != hfres.C_mo.shape[0];

    const dS = hfres.ints.dS orelse unreachable;
    const dH = hfres.ints.dH orelse unreachable;
    const dg = hfres.ints.dg orelse unreachable;

    const C, const P, const e = .{ hfres.C_mo, hfres.P_ao, hfres.e_mo };

    const nocc, const nbf = .{ if (generalized) hfres.ints.sys.nel else hfres.ints.sys.nel / 2, C.shape[0] };

    var F_x = try Matrix(T).init(nbf, nbf, gpa);
    defer F_x.deinit(gpa);

    var P_x = try Matrix(T).init(nbf, nbf, gpa);
    defer P_x.deinit(gpa);

    var U_x = try Matrix(T).init(nbf, nbf, gpa);
    defer U_x.deinit(gpa);

    var U_p = try Matrix(T).init(nbf, nbf, gpa);
    defer U_p.deinit(gpa);

    var V_x = try Matrix(T).init(nbf, nbf, gpa);
    defer V_x.deinit(gpa);

    var F_x_mo = try Matrix(T).init(nbf, nbf, gpa);
    defer F_x_mo.deinit(gpa);

    var S_x_mo = try Matrix(T).init(nbf, nbf, gpa);
    defer S_x_mo.deinit(gpa);

    var V_x_mo = try Matrix(T).init(nbf, nbf, gpa);
    defer V_x_mo.deinit(gpa);

    var dC_x = try Matrix(T).init(nbf, nbf, gpa);
    defer dC_x.deinit(gpa);

    const exch_factor: T = if (generalized) 1 else 0.5;
    const dens_factor: T = if (generalized) 1 else 2.0;

    for (0..3 * hfres.ints.sys.atoms.len) |c| {
        if (log) {
            try printf(io, "\nSOLVING CPHF FOR COORDINATE {d}\n", .{c + 1});
        }

        if (log) {
            try printf(io, "{s:4} {s:9} {s:9} {s:9}\n", .{ "ITER", "RMS(DU)", "RMS(DP)", "TIME" });
        }

        var diis_acc: ?Diis(T) = if (opt.diis) |d| if (d > 0) Diis(T).init(d) else null else null;
        defer if (diis_acc) |*d| d.deinit(gpa);

        F_x.zero();

        for (0..nbf) |mu| for (0..nbf) |nu| {
            F_x.ptr(mu, nu).* = dH.at(.{ c, mu, nu });
        };

        for (0..nbf) |lambda| for (0..nbf) |sigma| for (0..nbf) |mu| for (0..nbf) |nu| {
            const coul = dg.at(.{ c, lambda, mu, sigma, nu });
            const exch = dg.at(.{ c, lambda, sigma, mu, nu });

            F_x.ptr(mu, nu).* += P.at(lambda, sigma) * (coul - exch_factor * exch);
        };

        for (0..nbf) |mu| for (mu + 1..nbf) |nu| {
            const avg = (F_x.at(mu, nu) + F_x.at(nu, mu)) / 2;

            F_x.ptr(mu, nu).* = avg;
            F_x.ptr(nu, mu).* = avg;
        };

        const S_x = Matrix(T){ .data = dS.data[c * nbf * nbf .. (c + 1) * nbf * nbf], .shape = .{ nbf, nbf } };

        ao2mo_1e(T, &F_x_mo, F_x, C);
        ao2mo_1e(T, &S_x_mo, S_x, C);

        for (0..nbf) |p| for (0..nbf) |q| {
            U_x.ptr(p, q).* = -0.5 * S_x_mo.at(p, q);
        };

        for (nocc..nbf) |a| for (0..nocc) |i| {
            const diff = e.at(a) - e.at(i);

            if (@abs(diff) > 1e-12) {
                const u_ai = (e.at(i) * S_x_mo.at(a, i) - F_x_mo.at(a, i)) / diff;

                U_x.ptr(a, i).*, U_x.ptr(i, a).* = .{ u_ai, -S_x_mo.at(i, a) - u_ai };
            }
        };

        P_x.zero();

        for (0..opt.iterations) |i| {
            var timer = std.Io.Timestamp.now(io, .real);

            if (opt.diis) |d| if (d > 0) for (0..nbf) |p| for (0..nbf) |q| {
                U_p.ptr(p, q).* = U_x.at(p, q);
            };

            try mo2ao_1e(T, &dC_x, U_x, C);

            V_x.zero();

            var sum_sq_dp: T = 0;
            var sum_sq_du: T = 0;

            for (0..nbf) |mu| for (0..nbf) |nu| {
                var sum: T = 0;

                for (0..nocc) |k| {
                    sum += dens_factor * (dC_x.at(mu, k) * C.at(nu, k) + C.at(mu, k) * dC_x.at(nu, k));
                }

                sum_sq_dp += (sum - P_x.at(mu, nu)) * (sum - P_x.at(mu, nu));

                P_x.ptr(mu, nu).* = sum;
            };

            const rms_dp = @sqrt(sum_sq_dp / @as(T, @floatFromInt(nbf * nbf)));

            for (0..nbf) |mu| for (0..nbf) |nu| for (0..nbf) |lambda| for (0..nbf) |sigma| {
                const g_mlns = hfres.ints.g.?.at(.{ mu, lambda, nu, sigma });
                const g_mnls = hfres.ints.g.?.at(.{ mu, nu, lambda, sigma });

                V_x.ptr(lambda, sigma).* += P_x.at(mu, nu) * (g_mlns - exch_factor * g_mnls);
            };

            ao2mo_1e(T, &V_x_mo, V_x, C);

            for (nocc..nbf) |a| for (0..nocc) |k| {
                const prev, const delta_e = .{ U_x.at(a, k), e.at(a) - e.at(k) };

                const f = F_x_mo.at(a, k);
                const v = V_x_mo.at(a, k);
                const s = S_x_mo.at(a, k);

                U_x.ptr(a, k).* = if (@abs(delta_e) > 1e-12) (e.at(k) * s - f - v) / delta_e else -0.5 * s;

                sum_sq_du += (U_x.at(a, k) - prev) * (U_x.at(a, k) - prev);
            };

            const rms_du = @sqrt(sum_sq_du / @as(T, @floatFromInt((nbf - nocc) * nocc)));

            for (0..nocc) |k| for (nocc..nbf) |a| {
                U_x.ptr(k, a).* = -S_x_mo.at(k, a) - U_x.at(a, k);
            };

            if (diis_acc) |*d| {
                var e_diis = try Matrix(T).init(nbf, nbf, gpa);
                errdefer e_diis.deinit(gpa);

                for (0..nbf) |p| for (0..nbf) |q| {
                    e_diis.ptr(p, q).* = U_x.at(p, q) - U_p.at(p, q);
                };

                try d.update(&U_x, e_diis, false, gpa);

                for (0..nocc) |k| for (nocc..nbf) |a| {
                    U_x.ptr(k, a).* = -S_x_mo.at(k, a) - U_x.at(a, k);
                };
            }

            const elapsed = timer.untilNow(io, .real);

            if (log) {
                const fmt = "{d:4} {e:9.3} {e:9.3} {s:9} {s}\n";

                const time_str = try std.fmt.allocPrint(gpa, "{f}", .{elapsed});
                defer gpa.free(time_str);

                const diis_tag = if (diis_acc) |d| (if (d.err_hist.items.len >= 2) "DIIS" else "") else "";

                try printf(io, fmt, .{ i + 1, rms_du, rms_dp, time_str, diis_tag });
            }

            if (rms_du < opt.threshold and rms_dp < opt.threshold) {
                break;
            }
        } else return error.CphfDidNotConverge;

        for (0..nbf) |p| for (0..nbf) |q| {
            const delta_e = e.at(q) - e.at(p);

            const f = F_x_mo.at(p, q);
            const v = V_x_mo.at(p, q);
            const s = S_x_mo.at(p, q);

            U_x.ptr(p, q).* = if (@abs(delta_e) > 1e-12) (f + v - e.at(q) * s) / delta_e else -0.5 * s;
        };

        try mo2ao_1e(T, &dC_x, U_x, C);

        for (0..nbf) |mu| for (0..nbf) |p| {
            dC.ptr(.{ c, mu, p }).* = dC_x.at(mu, p);
        };

        for (0..nbf) |p| {
            de.ptr(c, p).* = F_x_mo.at(p, p) + V_x_mo.at(p, p) - S_x_mo.at(p, p) * e.at(p);
        }
    }
}