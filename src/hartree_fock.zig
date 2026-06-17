const std = @import("std");

const Allocator = std.mem.Allocator;

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
    gradient: bool = false,
    iterations: usize = 100,
    threshold: f64 = 1e-8,
};

// HARTREE-FOCK FUNCTIONS ==============================================================================================

pub fn getFock(comptime T: type, F: *Matrix(T), H: Matrix(T), P: Matrix(T), g: Tensor(T, 4), generalized: bool) !void {
    std.debug.assert(F.shape[0] == H.shape[0]);
    std.debug.assert(F.shape[1] == H.shape[1]);
    std.debug.assert(P.shape[0] == H.shape[0]);
    std.debug.assert(P.shape[1] == H.shape[1]);
    std.debug.assert(g.shape[0] == H.shape[0]);
    std.debug.assert(g.shape[1] == H.shape[0]);
    std.debug.assert(g.shape[2] == H.shape[0]);
    std.debug.assert(g.shape[3] == H.shape[0]);

    for (0..F.shape[0]) |i| for (0..F.shape[1]) |j| {
        F.ptr(i, j).* = H.at(i, j);
    };

    const exch_factor: T = if (generalized) 1.0 else 0.5;

    for (0..g.shape[0]) |i| for (0..g.shape[1]) |j| for (0..g.shape[2]) |k| for (0..g.shape[3]) |l| {
        F.ptr(k, l).* += P.at(i, j) * (g.at(.{ i, k, j, l }) - exch_factor * g.at(.{ i, j, k, l }));
    };

    for (0..F.shape[0]) |i| for (i + 1..F.shape[1]) |j| {
        const avg = (F.at(i, j) + F.at(j, i)) / 2;

        F.ptr(i, j).* = avg;
        F.ptr(j, i).* = avg;
    };
}

pub fn getDensity(comptime T: type, P: *Matrix(T), C: Matrix(T), nocc: usize, generalized: bool) void {
    std.debug.assert(C.shape[0] == P.shape[0]);
    std.debug.assert(C.shape[1] == P.shape[1]);

    const factor: T = if (generalized) 1 else 2;

    for (0..P.shape[0]) |i| for (0..P.shape[1]) |j| {
        var sum: T = 0;

        for (0..nocc) |k| {
            sum += factor * C.at(i, k) * C.at(j, k);
        }

        P.ptr(i, j).* = sum;
    };
}

pub fn getDensityAlloc(comptime T: type, C: Matrix(T), nocc: usize, generalized: bool, gpa: Allocator) !Matrix(T) {
    var P = try Matrix(T).initZero(C.shape[0], C.shape[1], gpa);

    getDensity(T, &P, C, nocc, generalized);

    return P;
}

pub fn getEnergy(comptime T: type, H: Matrix(T), F: Matrix(T), P: Matrix(T)) T {
    var energy: T = 0;

    for (0..P.shape[0]) |i| for (0..P.shape[1]) |j| {
        energy += 0.5 * P.at(i, j) * (H.at(i, j) + F.at(i, j));
    };

    return energy;
}

pub fn getDensRms(comptime T: type, P_old: Matrix(T), P_new: Matrix(T)) T {
    std.debug.assert(P_old.shape[0] == P_new.shape[0]);
    std.debug.assert(P_old.shape[1] == P_new.shape[1]);

    var sum_sq_diff: T = 0;

    for (0..P_old.shape[0]) |i| for (0..P_old.shape[1]) |j| {
        const diff = P_new.at(i, j) - P_old.at(i, j);

        sum_sq_diff += diff * diff;
    };

    return @sqrt(sum_sq_diff / @as(T, @floatFromInt(P_old.shape[0] * P_old.shape[1])));
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

    for (0..nbf) |r| for (0..nbf) |c| {
        var sum: T = 0;

        for (0..nbf) |k| {
            sum += F.at(r, k) * P.at(k, c);
        }

        FP.ptr(r, c).* = sum;
    };

    for (0..nbf) |r| for (0..nbf) |c| {
        var sum: T = 0;

        for (0..nbf) |k| {
            sum += FP.at(r, k) * S.at(k, c);
        }

        FPS.ptr(r, c).* = sum;
    };

    for (0..nbf) |r| for (0..nbf) |c| {
        err.ptr(r, c).* = FPS.at(r, c) - FPS.at(c, r);
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

        for (ei.data, ej.data) |val_i, val_j| {
            sum += val_i * val_j;
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

    for (0..fck_hist.len) |i| for (0..F.shape[0]) |r| for (0..F.shape[1]) |l| {
        F.ptr(r, l).* += c.at(i, 0) * fck_hist[i].at(r, l);
    };
}

// GRADIENT FUNCTIONS ==================================================================================================

pub fn gradient(comptime T: type, ints: Integrals(T), C: Matrix(T), P: Matrix(T), e: Vector(T), generalized: bool, gpa: Allocator) !Matrix(T) {
    const dS = ints.dS orelse unreachable;
    const dK = ints.dK orelse unreachable;
    const dV = ints.dV orelse unreachable;
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

    for (0..ints.sys.atoms.len) |A| for (0..3) |c| {
        var n_G: T = 0;
        var h_G: T = 0;
        var s_G: T = 0;
        var g_G: T = 0;

        for (0..ints.sys.atoms.len) |B| {
            if (B == A) continue;

            const Zi = @as(T, @floatFromInt(ints.sys.atoms[A]));
            const Zj = @as(T, @floatFromInt(ints.sys.atoms[B]));

            const dx = ints.sys.coors[3 * A + 0] - ints.sys.coors[3 * B + 0];
            const dy = ints.sys.coors[3 * A + 1] - ints.sys.coors[3 * B + 1];
            const dz = ints.sys.coors[3 * A + 2] - ints.sys.coors[3 * B + 2];

            const dist = std.math.sqrt(dx * dx + dy * dy + dz * dz);

            n_G -= Zi * Zj * (ints.sys.coors[3 * A + c] - ints.sys.coors[3 * B + c]) / (dist * dist * dist);
        }

        for (0..dS.shape[1]) |i| for (0..dS.shape[2]) |j| {
            const dk_val = dK.at(.{ 3 * A + c, i, j });
            const dv_val = dV.at(.{ 3 * A + c, i, j });
            const ds_val = dS.at(.{ 3 * A + c, i, j });

            h_G += P.at(i, j) * (dk_val + dv_val);

            s_G -= W.at(i, j) * ds_val;
        };

        for (0..dg.shape[1]) |i| for (0..dg.shape[2]) |j| for (0..dg.shape[3]) |k| for (0..dg.shape[4]) |l| {
            const dg1 = dg.at(.{ 3 * A + c, i, k, j, l });
            const dg2 = dg.at(.{ 3 * A + c, i, j, k, l });

            g_G += 0.5 * P.at(i, j) * P.at(k, l) * (dg1 - exch_factor * dg2);
        };

        G.ptr(A, c).* = h_G + g_G + s_G + n_G;
    };

    return G;
}

pub fn gradientCoef(comptime T: type, ints: Integrals(T), C: Matrix(T), P: Matrix(T), e: Vector(T), generalized: bool, gpa: Allocator) !Tensor(T, 3) {
    const dS = ints.dS orelse unreachable;
    const dK = ints.dK orelse unreachable;
    const dV = ints.dV orelse unreachable;
    const dg = ints.dg orelse unreachable;

    const nbf = C.shape[0];

    var dC = try Tensor(T, 3).initZero(.{ 3 * ints.sys.atoms.len, nbf, nbf }, gpa);
    errdefer dC.deinit(gpa);

    const exch_factor: T = if (generalized) 1 else 0.5;

    var F_x = try Matrix(T).init(nbf, nbf, gpa);
    defer F_x.deinit(gpa);

    var F_x_MO = try Matrix(T).init(nbf, nbf, gpa);
    defer F_x_MO.deinit(gpa);

    var S_x_MO = try Matrix(T).init(nbf, nbf, gpa);
    defer S_x_MO.deinit(gpa);

    var U_x = try Matrix(T).init(nbf, nbf, gpa);
    defer U_x.deinit(gpa);

    for (0..3 * ints.sys.atoms.len) |x| {
        F_x.zero();
        U_x.zero();

        for (0..nbf) |k| for (0..nbf) |l| {
            F_x.ptr(k, l).* = dK.at(.{ x, k, l }) + dV.at(.{ x, k, l });
        };

        for (0..nbf) |i| for (0..nbf) |j| for (0..nbf) |k| for (0..nbf) |l| {
            F_x.ptr(k, l).* += P.at(i, j) * (dg.at(.{ x, i, k, j, l }) - exch_factor * dg.at(.{ x, i, j, k, l }));
        };

        for (0..nbf) |k| for (k + 1..nbf) |l| {
            const avg = (F_x.at(k, l) + F_x.at(l, k)) / 2;

            F_x.ptr(k, l).* = avg;
            F_x.ptr(l, k).* = avg;
        };

        const S_x = Matrix(T){ .data = dS.data[x * nbf * nbf .. (x + 1) * nbf * nbf], .shape = .{ nbf, nbf } };

        ao2mo_pp(T, &F_x_MO, F_x, C);
        ao2mo_pp(T, &S_x_MO, S_x, C);

        for (0..nbf) |p| for (0..nbf) |q| {
            const diff = e.at(q) - e.at(p);

            U_x.ptr(p, q).* = if (@abs(diff) > 1e-12) (F_x_MO.at(p, q) - e.at(q) * S_x_MO.at(p, q)) / diff else -0.5 * S_x_MO.at(p, q);
        };

        var dC_x = Matrix(T){ .data = dC.data[x * nbf * nbf .. (x + 1) * nbf * nbf], .shape = .{ nbf, nbf } };

        mo2ao_xx(T, &dC_x, U_x, C);
    }

    return dC;
}

pub fn gradientOrben(comptime T: type, ints: Integrals(T), C: Matrix(T), P: Matrix(T), e: Vector(T), generalized: bool, gpa: Allocator) !Matrix(T) {
    const dS = ints.dS orelse unreachable;
    const dK = ints.dK orelse unreachable;
    const dV = ints.dV orelse unreachable;
    const dg = ints.dg orelse unreachable;

    const nbf = C.shape[0];

    var de = try Matrix(T).initZero(3 * ints.sys.atoms.len, nbf, gpa);
    errdefer de.deinit(gpa);

    const exch_factor: T = if (generalized) 1.0 else 0.5;

    var F_x = try Matrix(T).init(nbf, nbf, gpa);
    defer F_x.deinit(gpa);

    var F_x_MO = try Matrix(T).init(nbf, nbf, gpa);
    defer F_x_MO.deinit(gpa);

    var S_x_MO = try Matrix(T).init(nbf, nbf, gpa);
    defer S_x_MO.deinit(gpa);

    for (0..3 * ints.sys.atoms.len) |x| {
        F_x.zero();

        for (0..nbf) |k| for (0..nbf) |l| {
            F_x.ptr(k, l).* = dK.at(.{ x, k, l }) + dV.at(.{ x, k, l });
        };

        for (0..nbf) |i| for (0..nbf) |j| for (0..nbf) |k| for (0..nbf) |l| {
            F_x.ptr(k, l).* += P.at(i, j) * (dg.at(.{ x, i, k, j, l }) - exch_factor * dg.at(.{ x, i, j, k, l }));
        };

        for (0..nbf) |k| for (k + 1..nbf) |l| {
            const avg = (F_x.at(k, l) + F_x.at(l, k)) / 2;

            F_x.ptr(k, l).* = avg;
            F_x.ptr(l, k).* = avg;
        };

        const S_x = Matrix(T){ .data = dS.data[x * nbf * nbf .. (x + 1) * nbf * nbf], .shape = .{ nbf, nbf } };

        ao2mo_pp(T, &F_x_MO, F_x, C);
        ao2mo_pp(T, &S_x_MO, S_x, C);

        for (0..nbf) |p| {
            de.ptr(x, p).* = F_x_MO.at(p, p) - e.at(p) * S_x_MO.at(p, p);
        }
    }

    return de;
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

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.ints.deinit(gpa);

            self.C.deinit(gpa);
            self.P.deinit(gpa);
            self.F.deinit(gpa);
            self.e.deinit(gpa);

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
    const molopts = MolecularIntegralsOptions{
        .system = opt.system,
        .basis = opt.basis,
        .spin = opt.generalized,
        .calculate = .{
            .kinetic_d1 = opt.gradient,
            .overlap_d1 = opt.gradient,
            .coulomb_d1 = opt.gradient,
            .nuclear_d1 = opt.gradient,
        },
    };

    var ints = try molecular_integrals_run(T, io, molopts, log, gpa);
    errdefer ints.deinit(gpa);

    var energy = try gpa.alloc(T, 1);
    errdefer gpa.free(energy);

    var grad = try gpa.alloc(Matrix(T), if (opt.gradient) 1 else 0);
    errdefer if (opt.gradient) gpa.free(grad);

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

    var H = try Matrix(T).init(nbf, nbf, gpa);
    defer H.deinit(gpa);

    var B = try Matrix(T).init(nbf, nbf, gpa);
    defer B.deinit(gpa);

    var C = try Matrix(T).init(nbf, nbf, gpa);
    errdefer C.deinit(gpa);

    var F = try Matrix(T).init(nbf, nbf, gpa);
    errdefer F.deinit(gpa);

    var e = try Vector(T).init(nbf, gpa);
    errdefer e.deinit(gpa);

    for (0..nbf) |i| for (0..nbf) |j| {
        H.ptr(i, j).* = ints.K.?.at(i, j) + ints.V.?.at(i, j);
    };

    if (log) {
        try printf(io, "\nNUCLEAR REPULSION ENERGY: {d:.14} Eh\n", .{VN});
    }

    var P_old = try Matrix(T).initZero(nbf, nbf, gpa);
    defer P_old.deinit(gpa);

    var P_new = try Matrix(T).initZero(nbf, nbf, gpa);
    errdefer P_new.deinit(gpa);

    var e_old: T = VN;
    var e_new: T = VN;

    {
        @memcpy(B.data, ints.S.?.data);

        try geigh(T, &e, &C, H, &B);

        getDensity(T, &P_old, C, nocc, opt.generalized);
    }

    if (opt.iterations > 0 and log) {
        const fmt = "\nSELF CONSISTENT FIELD\n{s:4} {s:20} {s:9} {s:9} {s:9}\n";

        try printf(io, fmt, .{ "ITER", "TOTAL ENERGY", "|DE|", "RMS(DP)", "TIME" });
    }

    var fck_hist = std.ArrayList(Matrix(T)).empty;

    defer {
        for (fck_hist.items) |*mat| mat.deinit(gpa);

        fck_hist.deinit(gpa);
    }

    var err_hist = std.ArrayList(Matrix(T)).empty;

    defer {
        for (err_hist.items) |*mat| mat.deinit(gpa);

        err_hist.deinit(gpa);
    }

    for (0..opt.iterations) |i| {
        var timer = std.Io.Timestamp.now(io, .real);

        try getFock(T, &F, H, P_old, ints.g.?, opt.generalized);

        e_new = getEnergy(T, H, F, P_old) + VN;

        if (opt.diis != null and opt.diis.? > 0) {
            var f_diis = try Matrix(T).init(nbf, nbf, gpa);
            errdefer f_diis.deinit(gpa);

            var e_diis = try Matrix(T).init(nbf, nbf, gpa);
            errdefer e_diis.deinit(gpa);

            for (0..nbf) |j| for (0..nbf) |k| {
                f_diis.ptr(j, k).* = F.at(j, k);
            };

            try getError(T, &e_diis, F, P_old, ints.S.?, gpa);

            if (fck_hist.items.len >= opt.diis.?) {
                var old_f = fck_hist.orderedRemove(0);
                var old_e = err_hist.orderedRemove(0);

                old_f.deinit(gpa);
                old_e.deinit(gpa);
            }

            try fck_hist.append(gpa, f_diis);
            try err_hist.append(gpa, e_diis);

            diis(T, fck_hist.items, err_hist.items, &F, gpa) catch {
                for (fck_hist.items) |*mat| mat.deinit(gpa);
                for (err_hist.items) |*mat| mat.deinit(gpa);

                fck_hist.clearRetainingCapacity();
                err_hist.clearRetainingCapacity();

                var f_retry = try Matrix(T).init(nbf, nbf, gpa);
                var e_retry = try Matrix(T).init(nbf, nbf, gpa);

                for (0..nbf) |j| for (0..nbf) |k| {
                    f_retry.ptr(j, k).* = F.at(j, k);
                };

                try getError(T, &e_retry, F, P_old, ints.S.?, gpa);

                try fck_hist.append(gpa, f_retry);
                try err_hist.append(gpa, e_retry);
            };
        }

        @memcpy(B.data, ints.S.?.data);

        try geigh(T, &e, &C, F, &B);

        getDensity(T, &P_new, C, nocc, opt.generalized);

        const p_rm = getDensRms(T, P_old, P_new);
        const delta_energy = @abs(e_new - e_old);

        const elapsed = timer.untilNow(io, .real);

        if (log) {
            const fmt = "{d:4} {d:20.14} {e:9.3} {e:9.3} {s:9} {s}\n";

            const time_str = try std.fmt.allocPrint(gpa, "{f}", .{elapsed});
            defer gpa.free(time_str);

            const de = if (i == 0) 0 else delta_energy;

            try printf(io, fmt, .{ i + 1, e_new, de, p_rm, time_str, if (err_hist.items.len >= 2) "DIIS" else "" });
        }

        @memcpy(P_old.data, P_new.data);

        e_old = e_new;

        if (i > 0 and delta_energy < opt.threshold and p_rm < opt.threshold) {
            break;
        }

        if (i == opt.iterations - 1) {
            return error.ScfDidNotConverge;
        }
    }

    energy[0] = e_new;

    if (log) {
        try printf(io, "\nFINAL HARTREE-FOCK ENERGY: {d:.14} Eh\n", .{energy[0]});
    }

    if (opt.write.coefficients) |fname| {
        try writeMatrix(T, io, fname, C);
    }

    if (opt.write.density) |fname| {
        try writeMatrix(T, io, fname, P_new);
    }

    if (opt.write.fock) |fname| {
        try writeMatrix(T, io, fname, F);
    }

    if (opt.gradient) {
        grad[0] = try gradient(T, ints, C, P_new, e, opt.generalized, gpa);
    }

    errdefer {
        if (opt.gradient) grad[0].deinit(gpa);
    }

    if (log) for (grad) |G| {
        try std.Io.File.stdout().writeStreamingAll(io, "\nHARTREE-FOCK NUCLEAR ENERGY GRADIENT\n");

        for (0..G.shape[0]) |i| for (0..G.shape[1]) |j| {
            try printf(io, "{d:20.14}{s}", .{ G.at(i, j), if (j == 2) "\n" else " " });
        };
    };

    return Result(T){ .ints = ints, .P = P_new, .C = C, .F = F, .e = e, .energy = energy, .gradient = grad };
}
