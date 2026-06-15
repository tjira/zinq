const std = @import("std");

const Allocator = std.mem.Allocator;

const Integrals = @import("molecular_integrals.zig").Integrals;
const Matrix = @import("tensor.zig").Matrix;
const MolecularIntegralsOptions = @import("molecular_integrals.zig").Options;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

const molecular_integrals_run = @import("molecular_integrals.zig").run;
const geigh = @import("linear_algebra.zig").geigh;
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
    iterations: usize = 100,
    generalized: bool = false,
    threshold: f64 = 1e-8,
    gradient: bool = false,
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

pub fn gradient(comptime T: type, ints: Integrals(T), C: Matrix(T), P: Matrix(T), e: Vector(T), generalized: bool, gpa: Allocator) !Matrix(T) {
    const dS = ints.dS orelse return error.OverlapDerivativeMatrixNotCalculated;
    const dK = ints.dK orelse return error.KineticDerivativeMatrixNotCalculated;
    const dV = ints.dV orelse return error.NuclearDerivativeMatrixNotCalculated;
    const dg = ints.dg orelse return error.CoulombDerivativeMatrixNotCalculated;

    const nocc = if (generalized) ints.sys.nel else ints.sys.nel / 2;

    var G = try Matrix(T).initZero(ints.sys.atoms.len, 3, gpa);
    errdefer G.deinit(gpa);

    var W = try Matrix(T).initZero(P.shape[0], P.shape[0], gpa);
    defer W.deinit(gpa);

    const factor: T = if (generalized) 1 else 2;

    for (0..W.nrow()) |i| for (0..W.ncol()) |j| {
        var sum: T = 0;

        for (0..nocc) |k| {
            sum += factor * e.at(k) * C.at(i, k) * C.at(j, k);
        }

        W.ptr(i, j).* = sum;
    };

    const exch_factor: T = if (generalized) 1 else 0.5;

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

    const VN = try ints.sys.nrep();

    const nocc = if (opt.generalized) ints.sys.nel else blk: {
        if (ints.sys.nel % 2 != 0) {
            return error.OnlyClosedShellSupported;
        }

        break :blk ints.sys.nel / 2;
    };

    const nbf = if (opt.generalized) 2 * ints.sys.nbf else ints.sys.nbf;

    var energy = try gpa.alloc(T, 1);
    errdefer gpa.free(energy);

    var grad = try gpa.alloc(Matrix(T), if (opt.gradient) 1 else 0);
    errdefer if (opt.gradient) gpa.free(grad);

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
        const fmt = "\nSELF CONSISTENT FIELD\n{s:4} {s:20} {s:9} {s:9} {s:4}\n";

        try printf(io, fmt, .{ "ITER", "TOTAL ENERGY", "|DE|", "RMS(DP)", "TIME" });
    }

    for (0..opt.iterations) |i| {
        var timer = std.Io.Timestamp.now(io, .real);

        try getFock(T, &F, H, P_old, ints.g.?, opt.generalized);

        e_new = getEnergy(T, H, F, P_old) + VN;

        @memcpy(B.data, ints.S.?.data);

        try geigh(T, &e, &C, F, &B);

        getDensity(T, &P_new, C, nocc, opt.generalized);

        const p_rm = getDensRms(T, P_old, P_new);
        const delta_energy = @abs(e_new - e_old);

        if (log) {
            const fmt = "{d:4} {d:20.14} {e:9.3} {e:9.3} {f}\n";

            try printf(io, fmt, .{ i + 1, e_new, if (i > 0) delta_energy else 0, p_rm, timer.untilNow(io, .real) });
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

    if (log) for (grad) |G| {
        try std.Io.File.stdout().writeStreamingAll(io, "\nNUCLEAR ENERGY GRADIENT\n");

        for (0..G.shape[0]) |i| for (0..G.shape[1]) |j| {
            try printf(io, "{d:20.14}{s}", .{ G.at(i, j), if (j == 2) "\n" else " " });
        };
    };

    return Result(T){ .ints = ints, .P = P_new, .C = C, .F = F, .e = e, .energy = energy, .gradient = grad };
}
