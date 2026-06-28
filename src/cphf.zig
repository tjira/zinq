const std = @import("std");

const Allocator = std.mem.Allocator;

const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

const ao2mo_pp = @import("integral_transform.zig").ao2mo_pp;
const diis = @import("hartree_fock.zig").diis;
const mo2ao_xx = @import("integral_transform.zig").mo2ao_xx;
const printf = @import("read_write.zig").printf;

pub fn orbitalResponse(comptime T: type, io: std.Io, hfres: HartreeFockResult(T), opt: anytype, log: bool, gpa: Allocator) !struct { Tensor(T, 3), Matrix(T) } {
    const nbf = hfres.C.shape[0];

    var dC = try Tensor(T, 3).initZero(.{ 3 * hfres.ints.sys.atoms.len, nbf, nbf }, gpa);
    errdefer dC.deinit(gpa);

    var de = try Matrix(T).initZero(3 * hfres.ints.sys.atoms.len, nbf, gpa);
    errdefer de.deinit(gpa);

    try solveCPHF(T, io, &dC, &de, hfres, opt, log, gpa);

    return .{ dC, de };
}

pub fn solveCPHF(comptime T: type, io: std.Io, dC: *Tensor(T, 3), de: *Matrix(T), hfres: HartreeFockResult(T), opt: anytype, log: bool, gpa: Allocator) !void {
    const generalized = hfres.ints.sys.nbf != hfres.C.shape[0];

    const dS = hfres.ints.dS orelse unreachable;
    const dH = hfres.ints.dH orelse unreachable;
    const dg = hfres.ints.dg orelse unreachable;

    const C = hfres.C;
    const P = hfres.P;
    const e = hfres.e;

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

    var F_x_MO = try Matrix(T).init(nbf, nbf, gpa);
    defer F_x_MO.deinit(gpa);

    var S_x_MO = try Matrix(T).init(nbf, nbf, gpa);
    defer S_x_MO.deinit(gpa);

    var V_x_MO = try Matrix(T).init(nbf, nbf, gpa);
    defer V_x_MO.deinit(gpa);

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

        var uxm_hist = std.ArrayList(Matrix(T)).empty;

        defer {
            for (0..uxm_hist.items.len) |i| uxm_hist.items[i].deinit(gpa);

            uxm_hist.deinit(gpa);
        }

        var err_hist = std.ArrayList(Matrix(T)).empty;

        defer {
            for (0..err_hist.items.len) |i| err_hist.items[i].deinit(gpa);

            err_hist.deinit(gpa);
        }

        F_x.zero();

        for (0..nbf) |i| for (0..nbf) |j| {
            F_x.ptr(i, j).* = dH.at(.{ c, i, j });
        };

        for (0..nbf) |i| for (0..nbf) |j| for (0..nbf) |k| for (0..nbf) |l| {
            F_x.ptr(k, l).* += P.at(i, j) * (dg.at(.{ c, i, k, j, l }) - exch_factor * dg.at(.{ c, i, j, k, l }));
        };

        for (0..nbf) |i| for (i + 1..nbf) |j| {
            const avg = (F_x.at(i, j) + F_x.at(j, i)) / 2;

            F_x.ptr(i, j).* = avg;
            F_x.ptr(j, i).* = avg;
        };

        const S_x = Matrix(T){ .data = dS.data[c * nbf * nbf .. (c + 1) * nbf * nbf], .shape = .{ nbf, nbf } };

        ao2mo_pp(T, &F_x_MO, F_x, C);
        ao2mo_pp(T, &S_x_MO, S_x, C);

        for (0..nbf) |i| for (0..nbf) |j| {
            U_x.ptr(i, j).* = -0.5 * S_x_MO.at(i, j);
        };

        for (nocc..nbf) |a| for (0..nocc) |i| {
            const diff = e.at(a) - e.at(i);

            if (@abs(diff) > 1e-12) {
                const u_ai = (e.at(i) * S_x_MO.at(a, i) - F_x_MO.at(a, i)) / diff;

                U_x.ptr(a, i).*, U_x.ptr(i, a).* = .{ u_ai, -S_x_MO.at(i, a) - u_ai };
            }
        };

        P_x.zero();

        for (0..opt.iterations) |i| {
            var timer = std.Io.Timestamp.now(io, .real);

            if (opt.diis != null and opt.diis.? > 0) for (0..nbf) |j| for (0..nbf) |k| {
                U_p.ptr(j, k).* = U_x.at(j, k);
            };

            mo2ao_xx(T, &dC_x, U_x, C);

            V_x.zero();

            var sum_sq_dp: T = 0;
            var sum_sq_du: T = 0;

            for (0..nbf) |mu| for (0..nbf) |nu| {
                var sum: T = 0;

                for (0..nocc) |j| {
                    sum += dens_factor * (dC_x.at(mu, j) * C.at(nu, j) + C.at(mu, j) * dC_x.at(nu, j));
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

            ao2mo_pp(T, &V_x_MO, V_x, C);

            for (nocc..nbf) |a| for (0..nocc) |j| {
                const prev, const delta_e = .{ U_x.at(a, j), e.at(a) - e.at(j) };

                const f = F_x_MO.at(a, j);
                const v = V_x_MO.at(a, j);
                const s = S_x_MO.at(a, j);

                U_x.ptr(a, j).* = if (@abs(delta_e) > 1e-12) (e.at(j) * s - f - v) / delta_e else -0.5 * s;

                sum_sq_du += (U_x.at(a, j) - prev) * (U_x.at(a, j) - prev);
            };

            const rms_du = @sqrt(sum_sq_du / @as(T, @floatFromInt((nbf - nocc) * nocc)));

            for (0..nocc) |j| for (nocc..nbf) |a| {
                U_x.ptr(j, a).* = -S_x_MO.at(j, a) - U_x.at(a, j);
            };

            if (opt.diis != null and opt.diis.? > 0) {
                try uxm_hist.ensureUnusedCapacity(gpa, 1);
                try err_hist.ensureUnusedCapacity(gpa, 1);

                var u_diis_opt: ?Matrix(T) = try Matrix(T).init(nbf, nbf, gpa);
                errdefer if (u_diis_opt) |*m| m.deinit(gpa);

                var e_diis_opt: ?Matrix(T) = try Matrix(T).init(nbf, nbf, gpa);
                errdefer if (e_diis_opt) |*m| m.deinit(gpa);

                for (0..nbf) |j| for (0..nbf) |k| {
                    u_diis_opt.?.ptr(j, k).* = U_x.at(j, k);
                    e_diis_opt.?.ptr(j, k).* = U_x.at(j, k) - U_p.at(j, k);
                };

                if (uxm_hist.items.len >= opt.diis.?) {
                    var old_u = uxm_hist.orderedRemove(0);
                    var old_e = err_hist.orderedRemove(0);

                    old_u.deinit(gpa);
                    old_e.deinit(gpa);
                }

                uxm_hist.appendAssumeCapacity(u_diis_opt.?);
                err_hist.appendAssumeCapacity(e_diis_opt.?);

                diis(T, uxm_hist.items, err_hist.items, &U_x, false, gpa) catch {
                    var u_retry = try Matrix(T).init(nbf, nbf, gpa);
                    errdefer u_retry.deinit(gpa);

                    var e_retry = try Matrix(T).init(nbf, nbf, gpa);
                    errdefer e_retry.deinit(gpa);

                    @memcpy(u_retry.data, U_x.data);

                    for (0..nbf) |j| for (0..nbf) |k| {
                        e_retry.ptr(j, k).* = err_hist.items[err_hist.items.len - 1].at(j, k);
                    };

                    for (0..uxm_hist.items.len) |j| uxm_hist.items[j].deinit(gpa);
                    for (0..err_hist.items.len) |j| err_hist.items[j].deinit(gpa);

                    uxm_hist.clearRetainingCapacity();
                    err_hist.clearRetainingCapacity();

                    try uxm_hist.ensureUnusedCapacity(gpa, 1);
                    try err_hist.ensureUnusedCapacity(gpa, 1);

                    uxm_hist.appendAssumeCapacity(u_retry);
                    err_hist.appendAssumeCapacity(e_retry);
                };

                for (0..nocc) |j| for (nocc..nbf) |a| {
                    U_x.ptr(j, a).* = -S_x_MO.at(j, a) - U_x.at(a, j);
                };
            }

            const elapsed = timer.untilNow(io, .real);

            if (log) {
                const fmt = "{d:4} {e:9.3} {e:9.3} {s:9} {s}\n";

                const time_str = try std.fmt.allocPrint(gpa, "{f}", .{elapsed});
                defer gpa.free(time_str);

                try printf(io, fmt, .{ i + 1, rms_du, rms_dp, time_str, if (err_hist.items.len >= 2) "DIIS" else "" });
            }

            if (rms_du < opt.threshold and rms_dp < opt.threshold) {
                break;
            }
        } else return error.CphfDidNotConverge;

        for (0..nbf) |i| for (0..nbf) |j| {
            const delta_e = e.at(j) - e.at(i);

            const f = F_x_MO.at(i, j);
            const v = V_x_MO.at(i, j);
            const s = S_x_MO.at(i, j);

            U_x.ptr(i, j).* = if (@abs(delta_e) > 1e-12) (f + v - e.at(j) * s) / delta_e else -0.5 * s;
        };

        mo2ao_xx(T, &dC_x, U_x, C);

        for (0..nbf) |i| for (0..nbf) |q| {
            dC.ptr(.{ c, i, q }).* = dC_x.at(i, q);
        };

        for (0..nbf) |q| {
            de.ptr(c, q).* = F_x_MO.at(q, q) + V_x_MO.at(q, q) - S_x_MO.at(q, q) * e.at(q);
        }
    }
}
