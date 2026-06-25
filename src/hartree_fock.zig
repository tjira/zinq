const std = @import("std");

const Allocator = std.mem.Allocator;

const DftPotential = @import("density_functional_theory.zig").DftPotential;
const Integrals = @import("molecular_integrals.zig").Integrals;
const Matrix = @import("tensor.zig").Matrix;
const MolecularIntegralsOptions = @import("molecular_integrals.zig").Options;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

const ao2mo_pp = @import("integral_transform.zig").ao2mo_pp;
const geigh = @import("linear_algebra.zig").geigh;
const luFactorize = @import("linear_algebra.zig").luFactorize;
const luSolve = @import("linear_algebra.zig").luSolve;
const mo2ao_xx = @import("integral_transform.zig").mo2ao_xx;
const molecular_integrals_run = @import("molecular_integrals.zig").run;
const printf = @import("read_write.zig").printf;
const writeMatrix = @import("read_write.zig").writeMatrix;

// OPTIONS =============================================================================================================

const Write = struct {
    coefficients: ?[]const u8 = null,
    density: ?[]const u8 = null,
    fock: ?[]const u8 = null,
};

pub const Options = struct {
    system: []const u8,
    basis: []const u8,

    write: Write = .{},
    diis: ?u32 = 8,
    generalized: bool = false,
    iterations: u32 = 100,
    threshold: f64 = 1e-8,

    gradient: ?union(enum) {
        analytic: struct {},
    } = null,

    response: ?struct {
        iterations: u32 = 100,
        threshold: f64 = 1e-8,
    } = null,

    dft: ?struct {
        exchange: ?[]const u8 = null,
        correlation: ?[]const u8 = null,
        exchange_correlation: ?[]const u8 = null,

        grid: struct {
            radial: usize = 50,
            angular: usize = 302,
        } = .{},
    } = null,
};
// HARTREE-FOCK FUNCTIONS ==============================================================================================

pub fn getFock(comptime T: type, F: *Matrix(T), ints: Integrals(T), P: Matrix(T), generalized: bool, dft: ?*DftPotential(T), gpa: Allocator) !void {
    std.debug.assert(F.shape[0] == ints.H.?.shape[0]);
    std.debug.assert(F.shape[1] == ints.H.?.shape[1]);
    std.debug.assert(P.shape[0] == ints.H.?.shape[0]);
    std.debug.assert(P.shape[1] == ints.H.?.shape[1]);

    std.debug.assert(ints.g.?.shape[0] == ints.H.?.shape[0]);
    std.debug.assert(ints.g.?.shape[1] == ints.H.?.shape[0]);
    std.debug.assert(ints.g.?.shape[2] == ints.H.?.shape[0]);
    std.debug.assert(ints.g.?.shape[3] == ints.H.?.shape[0]);

    for (0..F.shape[0]) |i| for (0..F.shape[1]) |j| {
        F.ptr(i, j).* = ints.H.?.at(i, j);
    };

    if (dft) |pot| {
        for (0..ints.g.?.shape[0]) |j| for (0..ints.g.?.shape[1]) |k| for (0..ints.g.?.shape[2]) |l| for (0..ints.g.?.shape[3]) |m| {
            F.ptr(l, m).* += P.at(j, k) * ints.g.?.at(.{ j, l, k, m });
        };

        if (pot.exx_coef > 0) {
            const factor: T = if (generalized) 1 else 0.5;

            for (0..ints.g.?.shape[0]) |i| for (0..ints.g.?.shape[1]) |j| for (0..ints.g.?.shape[2]) |k| for (0..ints.g.?.shape[3]) |l| {
                F.ptr(k, l).* -= factor * pot.exx_coef * P.at(i, j) * ints.g.?.at(.{ i, j, k, l });
            };
        }

        try pot.evaluate(ints.sys, P, gpa);

        for (0..F.shape[0]) |j| for (0..F.shape[1]) |k| {
            F.ptr(j, k).* += pot.Vxc.at(j, k);
        };
    }

    if (dft == null) {
        const exch_factor: T = if (generalized) 1.0 else 0.5;

        for (0..ints.g.?.shape[0]) |i| for (0..ints.g.?.shape[1]) |j| for (0..ints.g.?.shape[2]) |k| for (0..ints.g.?.shape[3]) |l| {
            F.ptr(k, l).* += P.at(i, j) * (ints.g.?.at(.{ i, k, j, l }) - exch_factor * ints.g.?.at(.{ i, j, k, l }));
        };
    }

    for (0..F.shape[0]) |i| for (i + 1..F.shape[1]) |j| {
        const avg = (F.at(i, j) + F.at(j, i)) / 2;

        F.ptr(i, j).* = avg;
        F.ptr(j, i).* = avg;
    };
}

pub fn getDensity(comptime T: type, P: *Matrix(T), C: Matrix(T), nocc: usize, generalized: bool) T {
    std.debug.assert(C.shape[0] == P.shape[0]);
    std.debug.assert(C.shape[1] == P.shape[1]);

    const factor: T, var sum_sq_dp: T = .{ if (generalized) 1 else 2, 0 };

    for (0..P.shape[0]) |i| for (0..P.shape[1]) |j| {
        var sum: T = 0;

        for (0..nocc) |k| {
            sum += factor * C.at(i, k) * C.at(j, k);
        }

        sum_sq_dp += (sum - P.at(i, j)) * (sum - P.at(i, j));

        P.ptr(i, j).* = sum;
    };

    return @sqrt(sum_sq_dp / @as(T, @floatFromInt(P.shape[0] * P.shape[1])));
}

pub fn getEnergy(comptime T: type, ints: Integrals(T), F: Matrix(T), P: Matrix(T), dft: ?*DftPotential(T)) T {
    var energy: T = if (dft) |pot| pot.Exc else 0;

    if (dft) |pot| for (0..P.shape[0]) |j| for (0..P.shape[1]) |k| {
        energy += 0.5 * P.at(j, k) * (ints.H.?.at(j, k) + F.at(j, k) - pot.Vxc.at(j, k));
    };

    if (dft == null) for (0..P.shape[0]) |i| for (0..P.shape[1]) |j| {
        energy += 0.5 * P.at(i, j) * (ints.H.?.at(i, j) + F.at(i, j));
    };

    return energy;
}

// DIIS FUNCTIONS ======================================================================================================

pub fn getError(comptime T: type, err: *Matrix(T), F: Matrix(T), P: Matrix(T), S: Matrix(T), gpa: Allocator) !void {
    const nbf = F.shape[0];

    std.debug.assert(F.shape[1] == nbf);
    std.debug.assert(P.shape[0] == nbf);
    std.debug.assert(P.shape[1] == nbf);
    std.debug.assert(S.shape[0] == nbf);
    std.debug.assert(S.shape[1] == nbf);

    std.debug.assert(err.shape[0] == nbf);
    std.debug.assert(err.shape[1] == nbf);

    var FP = try Matrix(T).init(nbf, nbf, gpa);
    defer FP.deinit(gpa);

    var FPS = try Matrix(T).init(nbf, nbf, gpa);
    defer FPS.deinit(gpa);

    for (0..nbf) |i| for (0..nbf) |j| {
        var sum: T = 0;

        for (0..nbf) |k| {
            sum += F.at(i, k) * P.at(k, j);
        }

        FP.ptr(i, j).* = sum;
    };

    for (0..nbf) |i| for (0..nbf) |j| {
        var sum: T = 0;

        for (0..nbf) |k| {
            sum += FP.at(i, k) * S.at(k, j);
        }

        FPS.ptr(i, j).* = sum;
    };

    for (0..nbf) |i| for (0..nbf) |j| {
        err.ptr(i, j).* = FPS.at(i, j) - FPS.at(j, i);
    };
}

pub fn diis(comptime T: type, fck_hist: []const Matrix(T), err_hist: []const Matrix(T), F: *Matrix(T), gpa: Allocator) !void {
    if (fck_hist.len < 2) return;

    var B = try Matrix(T).initZero(fck_hist.len + 1, fck_hist.len + 1, gpa);
    defer B.deinit(gpa);

    var b = try Matrix(T).initZero(fck_hist.len + 1, 1, gpa);
    defer b.deinit(gpa);

    for (0..fck_hist.len) |i| for (i..fck_hist.len) |j| {
        var sum: T = 0;

        const ei = err_hist[i];
        const ej = err_hist[j];

        for (0..ei.data.len) |k| {
            sum += ei.data[k] * ej.data[k];
        }

        B.ptr(i, j).* = sum;
        B.ptr(j, i).* = sum;
    };

    for (0..fck_hist.len) |i| {
        B.ptr(i, fck_hist.len).* = -1.0;
        B.ptr(fck_hist.len, i).* = -1.0;
    }

    for (0..fck_hist.len) |i| {
        b.ptr(i, 0).* = 0.0;
    }

    b.ptr(fck_hist.len, 0).* = -1.0;

    const ipiv = try gpa.alloc(i32, fck_hist.len + 1);
    defer gpa.free(ipiv);

    try luFactorize(T, &B, ipiv);

    var c = try Matrix(T).init(fck_hist.len + 1, 1, gpa);
    defer c.deinit(gpa);

    try luSolve(T, &c, B, ipiv, b);

    F.zero();

    for (0..fck_hist.len) |i| for (0..F.shape[0]) |j| for (0..F.shape[1]) |k| {
        F.ptr(j, k).* += c.at(i, 0) * fck_hist[i].at(j, k);
    };

    for (0..F.shape[0]) |i| for (i + 1..F.shape[1]) |j| {
        const avg = (F.at(i, j) + F.at(j, i)) / 2;

        F.ptr(i, j).* = avg;
        F.ptr(j, i).* = avg;
    };
}

// GRADIENT FUNCTIONS ==================================================================================================

pub fn gradient(comptime T: type, ints: Integrals(T), C: Matrix(T), P: Matrix(T), e: Vector(T), generalized: bool, gpa: Allocator) !Matrix(T) {
    const dS = ints.dS orelse unreachable;
    const dH = ints.dH orelse unreachable;
    const dg = ints.dg orelse unreachable;

    const nocc = if (generalized) ints.sys.nel else ints.sys.nel / 2;

    var G = try Matrix(T).initZero(ints.sys.atoms.len, 3, gpa);
    errdefer G.deinit(gpa);

    const exch_factor: T = if (generalized) 1 else 0.5;

    var W = try Matrix(T).init(P.shape[0], P.shape[0], gpa);
    defer W.deinit(gpa);

    const factor: T = if (generalized) 1 else 2;

    for (0..W.nrow()) |i| for (0..W.ncol()) |j| {
        var sum: T = 0;

        for (0..nocc) |k| {
            sum += factor * e.at(k) * C.at(i, k) * C.at(j, k);
        }

        W.ptr(i, j).* = sum;
    };

    for (0..ints.sys.atoms.len) |i| for (0..3) |j| {
        var n_G: T = 0;
        var h_G: T = 0;
        var s_G: T = 0;
        var g_G: T = 0;

        for (0..ints.sys.atoms.len) |k| {
            if (k == i) continue;

            const Zi = @as(T, @floatFromInt(ints.sys.atoms[i]));
            const Zj = @as(T, @floatFromInt(ints.sys.atoms[k]));

            const dx = ints.sys.coors[3 * i + 0] - ints.sys.coors[3 * k + 0];
            const dy = ints.sys.coors[3 * i + 1] - ints.sys.coors[3 * k + 1];
            const dz = ints.sys.coors[3 * i + 2] - ints.sys.coors[3 * k + 2];

            const dist = std.math.sqrt(dx * dx + dy * dy + dz * dz);

            n_G -= Zi * Zj * (ints.sys.coors[3 * i + j] - ints.sys.coors[3 * k + j]) / (dist * dist * dist);
        }

        for (0..dS.shape[1]) |p| for (0..dS.shape[2]) |q| {
            const dh_val = dH.at(.{ 3 * i + j, p, q });
            const ds_val = dS.at(.{ 3 * i + j, p, q });

            h_G += P.at(p, q) * dh_val;
            s_G -= W.at(p, q) * ds_val;
        };

        for (0..dg.shape[1]) |p| for (0..dg.shape[2]) |q| for (0..dg.shape[3]) |r| for (0..dg.shape[4]) |s| {
            const dg1 = dg.at(.{ 3 * i + j, p, r, q, s });
            const dg2 = dg.at(.{ 3 * i + j, p, q, r, s });

            g_G += 0.5 * P.at(p, q) * P.at(r, s) * (dg1 - exch_factor * dg2);
        };

        G.ptr(i, j).* = h_G + g_G + s_G + n_G;
    };

    return G;
}

pub fn orbitalResponse(comptime T: type, io: std.Io, hfres: Result(T), opt: anytype, log: bool, gpa: Allocator) !struct { Tensor(T, 3), Matrix(T) } {
    const generalized = hfres.ints.sys.nbf != hfres.C.shape[0];

    const dS = hfres.ints.dS orelse unreachable;
    const dH = hfres.ints.dH orelse unreachable;
    const dg = hfres.ints.dg orelse unreachable;

    const C = hfres.C;
    const P = hfres.P;
    const e = hfres.e;

    const nocc, const nbf = .{ if (generalized) hfres.ints.sys.nel else hfres.ints.sys.nel / 2, C.shape[0] };

    var dC = try Tensor(T, 3).initZero(.{ 3 * hfres.ints.sys.atoms.len, nbf, nbf }, gpa);
    errdefer dC.deinit(gpa);

    var de = try Matrix(T).initZero(3 * hfres.ints.sys.atoms.len, nbf, gpa);
    errdefer de.deinit(gpa);

    var F_x = try Matrix(T).init(nbf, nbf, gpa);
    defer F_x.deinit(gpa);

    var P_x = try Matrix(T).init(nbf, nbf, gpa);
    defer P_x.deinit(gpa);

    var U_x = try Matrix(T).init(nbf, nbf, gpa);
    defer U_x.deinit(gpa);

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

                U_x.ptr(a, j).* = if (@abs(delta_e) > 1e-12)
                    (e.at(j) * S_x_MO.at(a, j) - F_x_MO.at(a, j) - V_x_MO.at(a, j)) / delta_e
                else
                    -0.5 * S_x_MO.at(a, j);

                sum_sq_du += (U_x.at(a, j) - prev) * (U_x.at(a, j) - prev);
            };

            const rms_du = @sqrt(sum_sq_du / @as(T, @floatFromInt((nbf - nocc) * nocc)));

            for (0..nocc) |j| for (nocc..nbf) |a| {
                U_x.ptr(j, a).* = -S_x_MO.at(j, a) - U_x.at(a, j);
            };

            const elapsed = timer.untilNow(io, .real);

            if (log) {
                const fmt = "{d:4} {e:9.3} {e:9.3} {s}\n";

                const time_str = try std.fmt.allocPrint(gpa, "{f}", .{elapsed});
                defer gpa.free(time_str);

                try printf(io, fmt, .{ i + 1, rms_du, rms_dp, time_str });
            }

            if (rms_du < opt.threshold and rms_dp < opt.threshold) {
                break;
            }
        } else return error.CphfDidNotConverge;

        for (0..nbf) |i| for (0..nbf) |j| {
            const delta_e = e.at(j) - e.at(i);

            U_x.ptr(i, j).* = if (@abs(delta_e) > 1e-12)
                (F_x_MO.at(i, j) + V_x_MO.at(i, j) - e.at(j) * S_x_MO.at(i, j)) / delta_e
            else
                -0.5 * S_x_MO.at(i, j);
        };

        mo2ao_xx(T, &dC_x, U_x, C);

        for (0..nbf) |mu| for (0..nbf) |q| {
            dC.ptr(.{ c, mu, q }).* = dC_x.at(mu, q);
        };

        for (0..nbf) |q| {
            de.ptr(c, q).* = F_x_MO.at(q, q) + V_x_MO.at(q, q) - S_x_MO.at(q, q) * e.at(q);
        }
    }

    return .{ dC, de };
}

// RESULT STRUCT =======================================================================================================

pub fn Result(comptime T: type) type {
    return struct {
        ints: Integrals(T),

        C: Matrix(T),
        P: Matrix(T),
        F: Matrix(T),
        e: Vector(T),

        energy: []T,

        gradient: []Matrix(T) = &.{},

        dC: ?Tensor(T, 3) = null,

        de: ?Matrix(T) = null,

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.ints.deinit(gpa);

            self.C.deinit(gpa);
            self.P.deinit(gpa);
            self.F.deinit(gpa);
            self.e.deinit(gpa);

            gpa.free(self.energy);

            for (0..self.gradient.len) |i| {
                self.gradient[i].deinit(gpa);
            }

            gpa.free(self.gradient);

            if (self.dC) |*dC| dC.deinit(gpa);

            if (self.de) |*de| de.deinit(gpa);
        }
    };
}

// RUN =================================================================================================================

pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator) !Result(T) {
    if (opt.dft != null and opt.gradient != null) {
        return error.DftGradientNotSupported;
    }

    const molopts = MolecularIntegralsOptions{
        .system = opt.system,
        .basis = opt.basis,
        .spin = opt.generalized,
        .calculate = .{
            .kinetic_d1 = opt.gradient != null,
            .overlap_d1 = opt.gradient != null,
            .coulomb_d1 = opt.gradient != null,
            .nuclear_d1 = opt.gradient != null,
            .hmatrix_d1 = opt.gradient != null,
        },
    };

    var ints = try molecular_integrals_run(T, io, molopts, log, gpa);
    errdefer ints.deinit(gpa);

    var energy = try gpa.alloc(T, 1);
    errdefer gpa.free(energy);

    var grad = try gpa.alloc(Matrix(T), if (opt.gradient) |_| 1 else 0);
    errdefer if (opt.gradient) |_| gpa.free(grad);

    const VN = try ints.sys.nrep();

    const nocc = if (opt.generalized) ints.sys.nel else blk: {
        if (ints.sys.nel % 2 != 0) {
            return error.OnlyClosedShellSupported;
        }

        break :blk ints.sys.nel / 2;
    };

    const nbf = if (opt.generalized) 2 * ints.sys.nbf else ints.sys.nbf;

    if (log) {
        try printf(io, "\nNUMBER OF BASIS FUNCTIONS: {d}, NUMBER OF OCCUPIED ORBITALS: {d}\n", .{ nbf, nocc });
    }

    var B = try Matrix(T).init(nbf, nbf, gpa);
    defer B.deinit(gpa);

    var C = try Matrix(T).init(nbf, nbf, gpa);
    errdefer C.deinit(gpa);

    var P = try Matrix(T).initZero(nbf, nbf, gpa);
    errdefer P.deinit(gpa);

    var F = try Matrix(T).init(nbf, nbf, gpa);
    errdefer F.deinit(gpa);

    var e = try Vector(T).init(nbf, gpa);
    errdefer e.deinit(gpa);

    var e_old: T = VN;
    var e_new: T = VN;

    var dft: ?DftPotential(T) = null;

    if (opt.dft) |dft_opt| {
        const n_rad, const n_leb = .{ dft_opt.grid.radial, dft_opt.grid.angular };

        const funcs = .{ dft_opt.exchange, dft_opt.correlation, dft_opt.exchange_correlation };

        dft = try DftPotential(T).init(ints.sys, funcs, n_rad, n_leb, opt.generalized, gpa);
    }

    defer if (dft) |*pot| {
        pot.deinit(gpa);
    };

    {
        @memcpy(B.data, ints.S.?.data);

        try geigh(T, &e, &C, ints.H.?, &B);

        _ = getDensity(T, &P, C, nocc, opt.generalized);
    }

    if (opt.iterations > 0 and log) {
        const fmt = "\nSELF CONSISTENT FIELD\n{s:4} {s:20} {s:9} {s:9} {s:9}\n";

        try printf(io, fmt, .{ "ITER", "TOTAL ENERGY", "|DE|", "RMS(DP)", "TIME" });
    }

    var fck_hist = std.ArrayList(Matrix(T)).empty;

    defer {
        for (0..fck_hist.items.len) |i| fck_hist.items[i].deinit(gpa);

        fck_hist.deinit(gpa);
    }

    var err_hist = std.ArrayList(Matrix(T)).empty;

    defer {
        for (0..err_hist.items.len) |i| err_hist.items[i].deinit(gpa);

        err_hist.deinit(gpa);
    }

    for (0..opt.iterations) |i| {
        var timer = std.Io.Timestamp.now(io, .real);

        if (dft) |*pot| {
            try getFock(T, &F, ints, P, opt.generalized, pot, gpa);

            e_new = getEnergy(T, ints, F, P, pot) + VN;
        }

        if (dft == null) {
            try getFock(T, &F, ints, P, opt.generalized, null, gpa);

            e_new = getEnergy(T, ints, F, P, null) + VN;
        }

        if (opt.diis != null and opt.diis.? > 0) {
            try fck_hist.ensureUnusedCapacity(gpa, 1);
            try err_hist.ensureUnusedCapacity(gpa, 1);

            var f_diis = try Matrix(T).init(nbf, nbf, gpa);
            errdefer f_diis.deinit(gpa);

            var e_diis = try Matrix(T).init(nbf, nbf, gpa);
            errdefer e_diis.deinit(gpa);

            for (0..nbf) |j| for (0..nbf) |k| {
                f_diis.ptr(j, k).* = F.at(j, k);
            };

            try getError(T, &e_diis, F, P, ints.S.?, gpa);

            if (fck_hist.items.len >= opt.diis.?) {
                var old_f = fck_hist.orderedRemove(0);
                var old_e = err_hist.orderedRemove(0);

                old_f.deinit(gpa);
                old_e.deinit(gpa);
            }

            fck_hist.appendAssumeCapacity(f_diis);
            err_hist.appendAssumeCapacity(e_diis);

            diis(T, fck_hist.items, err_hist.items, &F, gpa) catch {
                for (0..fck_hist.items.len) |j| fck_hist.items[j].deinit(gpa);
                for (0..err_hist.items.len) |j| err_hist.items[j].deinit(gpa);

                fck_hist.clearRetainingCapacity();
                err_hist.clearRetainingCapacity();

                try fck_hist.ensureUnusedCapacity(gpa, 1);
                try err_hist.ensureUnusedCapacity(gpa, 1);

                var f_retry = try Matrix(T).init(nbf, nbf, gpa);
                errdefer f_retry.deinit(gpa);

                var e_retry = try Matrix(T).init(nbf, nbf, gpa);
                errdefer e_retry.deinit(gpa);

                for (0..nbf) |j| for (0..nbf) |k| {
                    f_retry.ptr(j, k).* = F.at(j, k);
                };

                try getError(T, &e_retry, F, P, ints.S.?, gpa);

                fck_hist.appendAssumeCapacity(f_retry);
                err_hist.appendAssumeCapacity(e_retry);
            };
        }

        @memcpy(B.data, ints.S.?.data);

        try geigh(T, &e, &C, F, &B);

        const p_rm = getDensity(T, &P, C, nocc, opt.generalized);

        const delta_energy = @abs(e_new - e_old);

        const elapsed = timer.untilNow(io, .real);

        if (log) {
            const fmt = "{d:4} {d:20.14} {e:9.3} {e:9.3} {s:9} {s}\n";

            const time_str = try std.fmt.allocPrint(gpa, "{f}", .{elapsed});
            defer gpa.free(time_str);

            const de = if (i == 0) 0 else delta_energy;

            try printf(io, fmt, .{ i + 1, e_new, de, p_rm, time_str, if (err_hist.items.len >= 2) "DIIS" else "" });
        }

        e_old = e_new;

        if (i > 0 and delta_energy < opt.threshold and p_rm < opt.threshold) {
            break;
        }
    } else return error.ScfDidNotConverge;

    energy[0] = e_new;

    if (log) {
        try printf(io, "\nNUCLEAR REPULSION ENERGY: {d:.14} Eh\n", .{VN});
    }

    if (log) if (dft) |*pot| {
        const names = try pot.getFunctionalNames(gpa);
        defer gpa.free(names);

        try printf(io, "\nFINAL DFT ENERGY ({s}): {d:.14} Eh\n", .{ names, energy[0] });
    };

    if (log and dft == null) {
        try printf(io, "\nFINAL HARTREE-FOCK ENERGY: {d:.14} Eh\n", .{energy[0]});
    }

    if (opt.write.coefficients) |fname| {
        try writeMatrix(T, io, fname, C);
    }

    if (opt.write.density) |fname| {
        try writeMatrix(T, io, fname, P);
    }

    if (opt.write.fock) |fname| {
        try writeMatrix(T, io, fname, F);
    }

    if (opt.gradient) |_| {
        grad[0] = try gradient(T, ints, C, P, e, opt.generalized, gpa);
    }

    errdefer {
        if (opt.gradient) |_| grad[0].deinit(gpa);
    }

    if (log) for (0..grad.len) |i| {
        try std.Io.File.stdout().writeStreamingAll(io, "\nHARTREE-FOCK NUCLEAR ENERGY GRADIENT\n");

        for (0..grad[i].shape[0]) |j| for (0..grad[i].shape[1]) |k| {
            try printf(io, "{d:20.14}{s}", .{ grad[i].at(j, k), if (k == 2) "\n" else " " });
        };
    };

    var result: Result(T) = .{ .ints = ints, .P = P, .C = C, .F = F, .e = e, .energy = energy, .gradient = grad };

    if (opt.response) |response| {
        result.dC, result.de = try orbitalResponse(T, io, result, response, log, gpa);
    }

    return result;
}
