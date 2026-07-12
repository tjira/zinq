const std = @import("std");

const Allocator = std.mem.Allocator;

const DftPotential = @import("density_functional_theory.zig").DftPotential;
const Integrals = @import("molecular_integrals.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const MolecularIntegralsOptions = @import("molecular_integrals.zig").Options;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

const ao2mo_pp = @import("integral_transform.zig").ao2mo_pp;
const bfgs = @import("molecular_optimization.zig").bfgs;
const calculateHarmonicFrequencies = @import("frequency_analysis.zig").calculateHarmonicFrequencies;
const calculateNumericalGradient = @import("nuclear_derivative.zig").calculateNumericalGradient;
const calculateNumericalHessian = @import("nuclear_derivative.zig").calculateNumericalHessian;
const exportIfBuiltin = @import("molecular_integrals.zig").exportIfBuiltin;
const geigh = @import("linear_algebra.zig").geigh;
const getSymbol = @import("constant.zig").getSymbol;
const luFactorize = @import("linear_algebra.zig").luFactorize;
const luSolve = @import("linear_algebra.zig").luSolve;
const mo2ao_xx = @import("integral_transform.zig").mo2ao_xx;
const molecular_integrals_run = @import("molecular_integrals.zig").run;
const molecular_integrals_runFromSystem = @import("molecular_integrals.zig").runFromSystem;
const mulliken = @import("population_analysis.zig").mulliken;
const orbitalResponse = @import("cphf.zig").orbitalResponse;
const printHarmonicFrequencies = @import("frequency_analysis.zig").printHarmonicFrequencies;
const printMullikenCharges = @import("population_analysis.zig").printMullikenCharges;
const printf = @import("read_write.zig").printf;
const steepestDescent = @import("molecular_optimization.zig").steepestDescent;
const writeMatrix = @import("read_write.zig").writeMatrix;

const AN2SM = @import("constant.zig").AN2SM;
const AU2CM = @import("constant.zig").AU2CM;

pub const GradientOptions = union(enum) {
    analytic: struct {},
    numeric: struct {
        step: f64 = 1e-5,
    },
};

pub const Options = struct {
    system: []const u8,
    basis: []const u8,

    write: Write = .{},

    diis: ?u32 = 8,
    generalized: bool = false,
    charge: i32 = 0,
    multiplicity: u32 = 1,
    iterations: u32 = 100,
    threshold: f64 = 1e-8,
    mulliken: bool = false,

    dft: ?struct {
        exchange: ?[]const u8 = null,
        correlation: ?[]const u8 = null,
        exchange_correlation: ?[]const u8 = null,

        grid: struct {
            radial: usize = 50,
            angular: usize = 302,
        } = .{},
    } = null,

    gradient: ?GradientOptions = null,

    optimize: ?union(enum) {
        steepest_descent: struct {
            gradient: GradientOptions = .analytic,
            threshold: f64 = 1e-4,
            iterations: u32 = 100,
            step: f64 = 1e-1,
        },
        bfgs: struct {
            gradient: GradientOptions = .analytic,
            threshold: f64 = 1e-4,
            iterations: u32 = 100,
            step: f64 = 1,
        },
    } = null,

    hessian: ?union(enum) {
        numeric: struct {
            step: f64 = 1e-5,
        },
    } = null,

    response: ?struct {
        iterations: u32 = 100,
        threshold: f64 = 1e-8,
        diis: ?u32 = 8,
    } = null,
};

pub fn Result(comptime T: type) type {
    return struct {
        ints: Integrals(T),

        C: Matrix(T),
        P: Matrix(T),
        F: Matrix(T),
        e: Vector(T),

        energy: []T,

        grad: []Matrix(T) = &.{},
        hess: []Matrix(T) = &.{},

        dC: ?Tensor(T, 3) = null,

        de: ?Matrix(T) = null,

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.ints.deinit(gpa);

            self.C.deinit(gpa);
            self.P.deinit(gpa);
            self.F.deinit(gpa);
            self.e.deinit(gpa);

            gpa.free(self.energy);

            for (0..self.grad.len) |i| {
                self.grad[i].deinit(gpa);
            }

            gpa.free(self.grad);

            for (0..self.hess.len) |i| {
                self.hess[i].deinit(gpa);
            }

            gpa.free(self.hess);

            if (self.dC) |*dC| dC.deinit(gpa);
            if (self.de) |*de| de.deinit(gpa);
        }
    };
}

fn ScfWorkspace(comptime T: type) type {
    return struct {
        P: *Matrix(T),
        F: *Matrix(T),
        C: *Matrix(T),
        e: *Vector(T),
    };
}

const Write = struct {
    coefficients: ?[]const u8 = null,
    density: ?[]const u8 = null,
    fock: ?[]const u8 = null,
};

pub fn diis(comptime T: type, fck_hist: []const Matrix(T), err_hist: []const Matrix(T), F: *Matrix(T), symmetric: bool, gpa: Allocator) !void {
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

    if (symmetric) for (0..F.shape[0]) |i| for (i + 1..F.shape[1]) |j| {
        const avg = (F.at(i, j) + F.at(j, i)) / 2;

        F.ptr(i, j).* = avg;
        F.ptr(j, i).* = avg;
    };
}

pub fn gradient(comptime T: type, ints: Integrals(T), C: Matrix(T), P: Matrix(T), e: Vector(T), generalized: bool, dft: ?*DftPotential(T), gpa: Allocator) !Matrix(T) {
    const dS = ints.dS orelse unreachable;
    const dH = ints.dH orelse unreachable;
    const dg = ints.dg orelse unreachable;

    const nocc = if (generalized) ints.sys.nel else ints.sys.nel / 2;

    var G = try nuclearRepulsionGradient(T, ints.sys, gpa);
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

    const exx_val = if (dft) |d| d.exx_coef else 1;

    for (0..ints.sys.atoms.len) |i| for (0..3) |j| {
        var h_G: T = 0;
        var s_G: T = 0;
        var g_G: T = 0;

        for (0..dS.shape[1]) |p| for (0..dS.shape[2]) |q| {
            const dh_val = dH.at(.{ 3 * i + j, p, q });
            const ds_val = dS.at(.{ 3 * i + j, p, q });

            h_G += P.at(p, q) * dh_val;
            s_G -= W.at(p, q) * ds_val;
        };

        for (0..dg.shape[1]) |p| for (0..dg.shape[2]) |q| for (0..dg.shape[3]) |r| for (0..dg.shape[4]) |s| {
            const dg1 = dg.at(.{ 3 * i + j, p, r, q, s });
            const dg2 = dg.at(.{ 3 * i + j, p, q, r, s });

            g_G += 0.5 * P.at(p, q) * P.at(r, s) * (dg1 - exx_val * exch_factor * dg2);
        };

        G.ptr(i, j).* += h_G + g_G + s_G;
    };

    return G;
}

pub fn nuclearRepulsionGradient(comptime T: type, sys: MolecularSystem(T), gpa: Allocator) !Matrix(T) {
    var dVN = try Matrix(T).initZero(sys.atoms.len, 3, gpa);
    errdefer dVN.deinit(gpa);

    for (0..sys.atoms.len) |i| {
        const Zi = @as(T, @floatFromInt(sys.atoms[i]));

        const xi = sys.coors[3 * i + 0];
        const yi = sys.coors[3 * i + 1];
        const zi = sys.coors[3 * i + 2];

        for (0..sys.atoms.len) |k| {
            if (k == i) continue;

            const Zk = @as(T, @floatFromInt(sys.atoms[k]));

            const xk = sys.coors[3 * k + 0];
            const yk = sys.coors[3 * k + 1];
            const zk = sys.coors[3 * k + 2];

            const dx = xi - xk;
            const dy = yi - yk;
            const dz = zi - zk;

            const dist = std.math.sqrt(dx * dx + dy * dy + dz * dz);

            dVN.ptr(i, 0).* -= Zi * Zk * dx / (dist * dist * dist);
            dVN.ptr(i, 1).* -= Zi * Zk * dy / (dist * dist * dist);
            dVN.ptr(i, 2).* -= Zi * Zk * dz / (dist * dist * dist);
        }
    }

    return dVN;
}

pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator) !Result(T) {
    try checkInvalidInput(opt);

    const basis_path = try exportIfBuiltin(io, opt.basis, gpa);

    defer if (std.mem.startsWith(u8, opt.basis, "builtin:")) {
        std.Io.Dir.cwd().deleteFile(io, basis_path) catch {};
    };

    var sys = try MolecularSystem(T).init(opt.system, basis_path, opt.charge, opt.multiplicity, gpa);
    defer sys.deinit(gpa);

    if (std.mem.startsWith(u8, opt.basis, "builtin:")) {
        try std.Io.Dir.cwd().deleteFile(io, basis_path);
    }

    return try runFromSystem(T, io, opt, &sys, null, log, gpa);
}

pub fn runFromSystem(comptime T: type, io: std.Io, opt: Options, sys: *MolecularSystem(T), Pg: ?Matrix(T), log: bool, gpa: Allocator) !Result(T) {
    try checkInvalidInput(opt);

    var final_Pg: ?Matrix(T) = null;
    defer if (final_Pg) |*p| p.deinit(gpa);

    if (opt.optimize) |o| {
        switch (o) {
            .steepest_descent => final_Pg = try steepestDescent(T, io, runFromSystem, opt, sys, Pg, log, gpa),
            .bfgs => final_Pg = try bfgs(T, io, runFromSystem, opt, sys, Pg, log, gpa),
        }
    }

    const molopts = MolecularIntegralsOptions{
        .system = opt.system,
        .basis = opt.basis,
        .spin = opt.generalized,
        .charge = opt.charge,
        .multiplicity = opt.multiplicity,
        .calculate = .{
            .kinetic_d1 = opt.gradient != null and opt.gradient.? == .analytic,
            .overlap_d1 = opt.gradient != null and opt.gradient.? == .analytic,
            .coulomb_d1 = opt.gradient != null and opt.gradient.? == .analytic,
            .nuclear_d1 = opt.gradient != null and opt.gradient.? == .analytic,
            .hmatrix_d1 = opt.gradient != null and opt.gradient.? == .analytic,
        },
    };

    var ints = try molecular_integrals_runFromSystem(T, io, molopts, sys.*, log, gpa);
    errdefer ints.deinit(gpa);

    var energy = try gpa.alloc(T, 1);
    errdefer gpa.free(energy);

    var grad = try gpa.alloc(Matrix(T), if (opt.gradient) |_| 1 else 0);
    errdefer if (opt.gradient) |_| gpa.free(grad);

    const VN = try sys.nrep();

    const nocc = if (opt.generalized) sys.nel else blk: {
        if (sys.nel % 2 != 0) {
            return error.OnlyClosedShellSupported;
        }

        break :blk sys.nel / 2;
    };

    const nbf = if (opt.generalized) 2 * sys.nbf else sys.nbf;

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

    var dft: ?DftPotential(T) = null;

    if (opt.dft) |dft_opt| {
        const n_rad, const n_leb = .{ dft_opt.grid.radial, dft_opt.grid.angular };

        const funcs = .{ dft_opt.exchange, dft_opt.correlation, dft_opt.exchange_correlation };

        dft = try DftPotential(T).init(sys.*, funcs, n_rad, n_leb, opt.generalized, gpa);
    }

    defer if (dft) |*pot| {
        pot.deinit(gpa);
    };

    const guess_p = final_Pg orelse Pg;

    if (guess_p) |guess| {
        @memcpy(P.data, guess.data);
    }

    if (guess_p == null) {
        @memcpy(B.data, ints.S.?.data);

        try geigh(T, &e, &C, ints.H.?, &B);

        _ = getDensity(T, &P, C, nocc, opt.generalized);
    }

    const ws: ScfWorkspace(T) = .{ .P = &P, .F = &F, .C = &C, .e = &e };

    energy[0] = try scf(T, io, opt, ints, ws, if (dft) |*d| d else null, log, gpa);

    if (log and opt.mulliken) {
        var charges = try mulliken(T, sys.*, P, ints.S.?, gpa);
        defer charges.deinit(gpa);

        const method_str = if (dft) |_| "DFT" else "HARTREE-FOCK";

        try printMullikenCharges(T, io, sys.*, charges, method_str);
    }

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

    try exportMatrices(T, io, opt.write, C, P, F);

    if (opt.gradient) |gradopt| switch (gradopt) {
        .analytic => grad[0] = try gradient(T, ints, C, P, e, opt.generalized, if (dft) |*d| d else null, gpa),
        .numeric => grad[0] = try calculateNumericalGradient(T, io, runFromSystem, opt, sys, log, gpa),
    };

    errdefer {
        if (opt.gradient) |_| grad[0].deinit(gpa);
    }

    if (log) for (0..grad.len) |i| {
        const grad_type_str = if (opt.gradient.? == .analytic) "ANALYTICAL" else "NUMERICAL";

        const method_str = if (dft) |_| "DFT" else "HARTREE-FOCK";

        try printf(io, "\n{s} {s} NUCLEAR ENERGY GRADIENT (Eh/a0)\n", .{ method_str, grad_type_str });

        for (0..grad[i].shape[0]) |j| for (0..grad[i].shape[1]) |k| {
            try printf(io, "{d:20.14}{s}", .{ grad[i].at(j, k), if (k == 2) "\n" else " " });
        };
    };

    const hess = try handleHessianAndFrequencies(T, io, opt, runFromSystem, sys, log, gpa);

    errdefer {
        if (opt.hessian) |_| hess[0].deinit(gpa);

        gpa.free(hess);
    }

    var result: Result(T) = .{ .ints = ints, .P = P, .C = C, .F = F, .e = e, .energy = energy, .grad = grad, .hess = hess };

    if (opt.response) |response| {
        result.dC, result.de = try orbitalResponse(T, io, result, response, log, gpa);
    }

    return result;
}

fn checkInvalidInput(opt: Options) !void {
    if (opt.multiplicity == 0) {
        std.log.err("MULTIPLICITY MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.multiplicity > 1 and !opt.generalized) {
        std.log.err("OPEN SHELL IS NOT SUPPORTED UNLESS GENERALIZED HARTREE-FOCK IS ENABLED", .{});

        return error.OnlyClosedShellSupported;
    }

    if (opt.system.len == 0) {
        std.log.err("MOLECULAR SYSTEM XYZ PATH IS EMPTY", .{});

        return error.InvalidInput;
    }

    if (opt.basis.len == 0) {
        std.log.err("BASIS SET G94 PATH IS EMPTY", .{});

        return error.InvalidInput;
    }

    if (opt.iterations == 0) {
        std.log.err("SCF ITERATIONS MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.threshold <= 0) {
        std.log.err("SCF CONVERGENCE THRESHOLD MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    if (opt.dft) |d| {
        if (opt.gradient != null and opt.gradient.? == .analytic) {
            std.log.err("ANALYTIC GRADIENTS ARE NOT SUPPORTED FOR DFT CALCULATIONS", .{});

            return error.InvalidInput;
        }

        if (d.grid.radial == 0) {
            std.log.err("DFT RADIAL GRID POINTS MUST BE GREATER THAN 0", .{});

            return error.InvalidInput;
        }

        const valid_angular_points = switch (d.grid.angular) {
            50, 74, 110, 302, 590, 974 => true,

            else => false,
        };

        if (!valid_angular_points) {
            std.log.err("DFT ANGULAR GRID SIZE IS NOT SUPPORTED. USE 50, 74, 110, 302, 590, OR 974", .{});

            return error.InvalidInput;
        }

        if (d.exchange_correlation != null and (d.exchange != null or d.correlation != null)) {
            std.log.err("CANNOT SPECIFY BOTH EXCHANGE_CORRELATION AND INDIVIDUAL FUNCTIONALS", .{});

            return error.InvalidInput;
        }

        if (d.exchange_correlation == null and d.exchange == null and d.correlation == null) {
            std.log.err("MUST SPECIFY EITHER EXCHANGE_CORRELATION OR EXCHANGE AND/OR CORRELATION FUNCTIONALS", .{});

            return error.InvalidInput;
        }
    }

    if (opt.response) |r| {
        if (r.iterations == 0) {
            std.log.err("CPHF RESPONSE ITERATIONS MUST BE GREATER THAN 0", .{});

            return error.InvalidInput;
        }

        if (r.threshold <= 0) {
            std.log.err("CPHF RESPONSE THRESHOLD MUST BE GREATER THAN 0", .{});

            return error.InvalidInput;
        }
    }
}

fn getDensity(comptime T: type, P: *Matrix(T), C: Matrix(T), nocc: usize, generalized: bool) T {
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

fn getEnergy(comptime T: type, ints: Integrals(T), F: Matrix(T), P: Matrix(T), dft: ?*DftPotential(T)) T {
    var energy: T = if (dft) |pot| pot.Exc else 0;

    if (dft) |pot| for (0..P.shape[0]) |j| for (0..P.shape[1]) |k| {
        energy += 0.5 * P.at(j, k) * (ints.H.?.at(j, k) + F.at(j, k) - pot.Vxc.at(j, k));
    };

    if (dft == null) for (0..P.shape[0]) |i| for (0..P.shape[1]) |j| {
        energy += 0.5 * P.at(i, j) * (ints.H.?.at(i, j) + F.at(i, j));
    };

    return energy;
}

fn getError(comptime T: type, err: *Matrix(T), F: Matrix(T), P: Matrix(T), S: Matrix(T), gpa: Allocator) !void {
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

fn getFock(comptime T: type, F: *Matrix(T), ints: Integrals(T), P: Matrix(T), generalized: bool, dft: ?*DftPotential(T), gpa: Allocator) !void {
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

fn scf(comptime T: type, io: std.Io, opt: Options, ints: Integrals(T), ws: ScfWorkspace(T), dft: ?*DftPotential(T), log: bool, gpa: Allocator) !T {
    const VN = try ints.sys.nrep();

    const nocc = if (opt.generalized) ints.sys.nel else blk: {
        if (ints.sys.nel % 2 != 0) {
            return error.OnlyClosedShellSupported;
        }

        break :blk ints.sys.nel / 2;
    };

    const nbf = if (opt.generalized) 2 * ints.sys.nbf else ints.sys.nbf;

    var e_old: T = VN;
    var e_new: T = VN;

    var B = try Matrix(T).init(nbf, nbf, gpa);
    defer B.deinit(gpa);

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

        if (dft) |pot| {
            try getFock(T, ws.F, ints, ws.P.*, opt.generalized, pot, gpa);

            e_new = getEnergy(T, ints, ws.F.*, ws.P.*, pot) + VN;
        }

        if (dft == null) {
            try getFock(T, ws.F, ints, ws.P.*, opt.generalized, null, gpa);

            e_new = getEnergy(T, ints, ws.F.*, ws.P.*, null) + VN;
        }

        if (opt.diis != null and opt.diis.? > 0) {
            try fck_hist.ensureUnusedCapacity(gpa, 1);
            try err_hist.ensureUnusedCapacity(gpa, 1);

            var f_diis = try Matrix(T).init(nbf, nbf, gpa);
            errdefer f_diis.deinit(gpa);

            var e_diis = try Matrix(T).init(nbf, nbf, gpa);
            errdefer e_diis.deinit(gpa);

            for (0..nbf) |j| for (0..nbf) |k| {
                f_diis.ptr(j, k).* = ws.F.at(j, k);
            };

            try getError(T, &e_diis, ws.F.*, ws.P.*, ints.S.?, gpa);

            if (fck_hist.items.len >= opt.diis.?) {
                var old_f = fck_hist.orderedRemove(0);
                var old_e = err_hist.orderedRemove(0);

                old_f.deinit(gpa);
                old_e.deinit(gpa);
            }

            fck_hist.appendAssumeCapacity(f_diis);
            err_hist.appendAssumeCapacity(e_diis);

            diis(T, fck_hist.items, err_hist.items, ws.F, true, gpa) catch {
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
                    f_retry.ptr(j, k).* = ws.F.at(j, k);
                };

                try getError(T, &e_retry, ws.F.*, ws.P.*, ints.S.?, gpa);

                fck_hist.appendAssumeCapacity(f_retry);
                err_hist.appendAssumeCapacity(e_retry);
            };
        }

        @memcpy(B.data, ints.S.?.data);

        try geigh(T, ws.e, ws.C, ws.F.*, &B);

        const delta_energy, const p_rm = .{ @abs(e_new - e_old), getDensity(T, ws.P, ws.C.*, nocc, opt.generalized) };

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

    return e_new;
}

fn handleHessianAndFrequencies(comptime T: type, io: std.Io, opt: Options, runFn: anytype, sys: *MolecularSystem(T), log: bool, gpa: Allocator) ![]Matrix(T) {
    var hess = try gpa.alloc(Matrix(T), if (opt.hessian) |_| 1 else 0);
    errdefer if (opt.hessian) |_| gpa.free(hess);

    if (opt.hessian) |hessopt| switch (hessopt) {
        .numeric => hess[0] = try calculateNumericalHessian(T, io, runFn, opt, sys, log, gpa),
    };

    errdefer if (opt.hessian) |_| hess[0].deinit(gpa);

    if (log and opt.hessian != null) {
        var freqs = try calculateHarmonicFrequencies(T, hess[0], sys.atoms, gpa);
        defer freqs.deinit(gpa);

        const method_str = if (opt.dft != null) "DFT NUMERIC" else "HARTREE-FOCK NUMERIC";

        try printHarmonicFrequencies(T, io, freqs, method_str);
    }

    return hess;
}

fn exportMatrices(comptime T: type, io: std.Io, write: Write, C: Matrix(T), P: Matrix(T), F: Matrix(T)) !void {
    if (write.coefficients) |fname| {
        try writeMatrix(T, io, fname, C);
    }

    if (write.density) |fname| {
        try writeMatrix(T, io, fname, P);
    }

    if (write.fock) |fname| {
        try writeMatrix(T, io, fname, F);
    }
}
