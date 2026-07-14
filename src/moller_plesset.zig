//! Computes Møller-Plesset correlation energy corrections to a Hartree-Fock reference wavefunction.

const std = @import("std");

const Allocator = std.mem.Allocator;

const HartreeFockOptions = @import("hartree_fock.zig").Options;
const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const ScalarDual = @import("dual.zig").ScalarDual;
const Tensor = @import("tensor.zig").Tensor;
const Value = @import("value.zig").Value;
const Vector = @import("tensor.zig").Vector;

const ao2mo_pppp = @import("integral_transform.zig").ao2mo_pppp;
const ao2so_coef = @import("integral_transform.zig").ao2so_coef;
const ao2so_pppp = @import("integral_transform.zig").ao2so_pppp;
const bfgs = @import("molecular_optimization.zig").bfgs;
const calculateHarmonicFrequencies = @import("frequency_analysis.zig").calculateHarmonicFrequencies;
const calculateNumericalGradient = @import("nuclear_derivative.zig").calculateNumericalGradient;
const calculateNumericalHessian = @import("nuclear_derivative.zig").calculateNumericalHessian;
const exportIfBuiltin = @import("molecular_integrals.zig").exportIfBuiltin;
const generateDets = @import("configuration_interaction.zig").generateDets;
const hartree_fock_run = @import("hartree_fock.zig").run;
const hartree_fock_runFromSystem = @import("hartree_fock.zig").runFromSystem;
const printHarmonicFrequencies = @import("frequency_analysis.zig").printHarmonicFrequencies;
const printf = @import("read_write.zig").printf;
const writeMatrix = @import("read_write.zig").writeMatrix;
const writeXyzFile = @import("read_write.zig").writeXyzFile;
const slater = @import("configuration_interaction.zig").slater;
const steepestDescent = @import("molecular_optimization.zig").steepestDescent;

const AU2CM = @import("constant.zig").AU2CM;

/// Options for computing Moller-Plesset energy gradients analytically or numerically.
pub const GradientOptions = union(enum) {
    analytic: struct {},
    numeric: struct {
        step: f64 = 1e-5,
    },
};

/// Parameters governing the Møller-Plesset perturbation theory calculation and its derivatives.
pub const Options = struct {
    hartree_fock: HartreeFockOptions,

    order: u32 = 2,

    write: Write = .{},

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
};

/// Holds Moller-Plesset calculation outputs: HF reference result, perturbation energies, and derivatives.
pub fn Result(comptime T: type) type {
    return struct {
        hartree_fock: HartreeFockResult(T),

        energy: []T,

        grad: []Matrix(T),
        hess: []Matrix(T),

        /// Frees all allocated resources associated with the Møller-Plesset result.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.hartree_fock.deinit(gpa);

            gpa.free(self.energy);

            for (0..self.grad.len) |i| {
                self.grad[i].deinit(gpa);
            }

            gpa.free(self.grad);

            for (0..self.hess.len) |i| {
                self.hess[i].deinit(gpa);
            }

            gpa.free(self.hess);
        }
    };
}

/// Destination paths for outputting Moller-Plesset geometries, gradients, and Hessians.
const Write = struct {
    geometry: ?[]const u8 = null,
    gradient: ?[]const u8 = null,
    hessian: ?[]const u8 = null,
};

/// Executes a Møller-Plesset calculation on a molecular system specified by file paths.
pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator) !Result(T) {
    try checkInvalidInput(opt);

    const basis_path = try exportIfBuiltin(io, opt.hartree_fock.basis, gpa);

    defer if (std.mem.startsWith(u8, opt.hartree_fock.basis, "builtin:")) {
        std.Io.Dir.cwd().deleteFile(io, basis_path) catch {};
    };

    var sys = try MolecularSystem(T).init(opt.hartree_fock.system, basis_path, opt.hartree_fock.charge, opt.hartree_fock.multiplicity, gpa);
    defer sys.deinit(gpa);

    if (std.mem.startsWith(u8, opt.hartree_fock.basis, "builtin:")) {
        try std.Io.Dir.cwd().deleteFile(io, basis_path);
    }

    return try runFromSystem(T, io, opt, &sys, null, log, gpa);
}

/// Runs a Møller-Plesset calculation starting from an initialized MolecularSystem structure.
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

    try checkInvalidInput(opt);

    const generalized, var hf_opt = .{ opt.hartree_fock.generalized, opt.hartree_fock };

    if (opt.gradient != null and opt.gradient.? == .analytic) {
        hf_opt.gradient = .{ .analytic = .{} };

        if (hf_opt.response == null) {
            hf_opt.response = .{};
        }
    }

    var hfres = try hartree_fock_runFromSystem(T, io, hf_opt, sys, final_Pg orelse Pg, log, gpa);
    errdefer hfres.deinit(gpa);

    var energy = try gpa.alloc(T, 1);
    errdefer gpa.free(energy);

    const grad = try gpa.alloc(Matrix(T), if (opt.gradient) |_| 1 else 0);
    errdefer gpa.free(grad);

    energy[0] = hfres.energy[0];

    if (opt.gradient) |gradopt| if (gradopt == .analytic) {
        grad[0] = try hfres.grad[0].clone(gpa);
    };

    errdefer if (opt.gradient != null and opt.gradient.? == .analytic) {
        grad[0].deinit(gpa);
    };

    const nocc = if (generalized) hfres.ints.sys.nel else hfres.ints.sys.nel / 2;

    const energies = try mp(T, opt.order, hfres.ints.g.?, hfres.C, hfres.e, nocc, generalized, gpa);
    defer gpa.free(energies);

    const grads = if (opt.gradient != null and opt.gradient.? == .analytic) try gradient(T, opt.order, hfres, gpa) else null;

    defer if (grads != null) {
        for (0..(opt.order + 1)) |i| {
            grads.?[i].deinit(gpa);
        }

        gpa.free(grads.?);
    };

    for (2..(opt.order + 1)) |i| {
        energy[0] += energies[i];

        if (log) {
            try printf(io, "\nMP{d} TOTAL ENERGY: {d:.14} Eh\n", .{ i, energy[0] });
        }

        if (grads != null) {
            for (0..grad[0].nrow()) |j| for (0..grad[0].ncol()) |k| {
                grad[0].ptr(j, k).* += grads.?[i].at(j, k);
            };

            if (log) {
                try printf(io, "\nMP{d} ANALYTICAL NUCLEAR ENERGY GRADIENT (Eh/a0)\n", .{i});

                for (0..grad[0].nrow()) |j| for (0..grad[0].ncol()) |k| {
                    try printf(io, "{d:20.14}{s}", .{ grad[0].at(j, k), if (k == 2) "\n" else " " });
                };
            }
        }
    }

    if (opt.gradient) |gradopt| if (gradopt == .numeric) {
        grad[0] = try calculateNumericalGradient(T, io, runFromSystem, opt, sys, log, gpa);
        errdefer grad[0].deinit(gpa);

        if (log) {
            try printf(io, "\nMP{d} NUMERICAL NUCLEAR ENERGY GRADIENT (Eh/a0)\n", .{opt.order});

            for (0..grad[0].nrow()) |j| for (0..grad[0].ncol()) |k| {
                try printf(io, "{d:20.14}{s}", .{ grad[0].at(j, k), if (k == 2) "\n" else " " });
            };
        }
    };

    errdefer if (opt.gradient != null and opt.gradient.? == .numeric) {
        grad[0].deinit(gpa);
    };

    const hess = try handleHessianAndFrequencies(T, io, opt, runFromSystem, sys, log, gpa);

    errdefer {
        if (opt.hessian) |_| hess[0].deinit(gpa);

        gpa.free(hess);
    }

    if (opt.write.gradient) |fname| if (grad.len > 0) {
        try writeMatrix(T, io, fname, grad[0]);
    };

    if (opt.write.hessian) |fname| if (hess.len > 0) {
        try writeMatrix(T, io, fname, hess[0]);
    };

    if (opt.write.geometry) |fname| {
        try writeXyzFile(T, io, fname, sys.atoms, sys.coors);
    }

    return Result(T){ .hartree_fock = hfres, .energy = energy, .grad = grad, .hess = hess };
}

/// Validates Moller-Plesset inputs for method requirements and minimum perturbation order.
fn checkInvalidInput(opt: Options) !void {
    if (opt.write.gradient != null and opt.gradient == null) {
        std.log.err("GRADIENT WRITE REQUESTED BUT GRADIENT IS NOT CALCULATED", .{});

        return error.InvalidInput;
    }

    if (opt.write.hessian != null and opt.hessian == null) {
        std.log.err("HESSIAN WRITE REQUESTED BUT HESSIAN IS NOT CALCULATED", .{});

        return error.InvalidInput;
    }

    if (opt.order < 2) {
        std.log.err("MØLLER-PLESSET PERTURBATION ORDER MUST BE AT LEAST 2 (MP2)", .{});

        return error.InvalidInput;
    }
}

/// Computes the nuclear gradient of the Møller-Plesset energy corrections using dual number differentiation.
fn gradient(comptime T: type, order: usize, hfres: HartreeFockResult(T), gpa: Allocator) ![]Matrix(T) {
    const nbf, const generalized = .{ hfres.C.shape[0], hfres.ints.sys.nbf != hfres.C.shape[0] };

    const nocc = if (generalized) hfres.ints.sys.nel else hfres.ints.sys.nel / 2;

    const grads = try gpa.alloc(Matrix(T), order + 1);
    errdefer gpa.free(grads);

    for (0..order + 1) |k| {
        grads[k] = Matrix(T).initZero(hfres.ints.sys.atoms.len, 3, gpa) catch |err| {
            for (0..k) |i| {
                grads[i].deinit(gpa);
            }

            return err;
        };
    }

    for (0..3 * hfres.ints.sys.atoms.len) |i| {
        var C_d = try Matrix(ScalarDual(T)).init(nbf, nbf, gpa);
        defer C_d.deinit(gpa);

        for (0..nbf) |mu| for (0..nbf) |p| {
            C_d.ptr(mu, p).* = ScalarDual(T).init(hfres.C.at(mu, p), hfres.dC.?.at(.{ i, mu, p }));
        };

        var e_d = try Vector(ScalarDual(T)).init(nbf, gpa);
        defer e_d.deinit(gpa);

        for (0..nbf) |j| {
            e_d.ptr(j).* = ScalarDual(T).init(hfres.e.at(j), hfres.de.?.at(i, j));
        }

        var g_d = try Tensor(ScalarDual(T), 4).init(hfres.ints.g.?.shape, gpa);
        defer g_d.deinit(gpa);

        for (0..nbf) |mu| for (0..nbf) |lambda| for (0..nbf) |nu| for (0..nbf) |sigma| {
            const dg_x = hfres.ints.dg.?.at(.{ i, mu, lambda, nu, sigma });

            g_d.ptr(.{ mu, lambda, nu, sigma }).* = ScalarDual(T).init(hfres.ints.g.?.at(.{ mu, lambda, nu, sigma }), dg_x);
        };

        const corr_energies_d = try mp(ScalarDual(T), order, g_d, C_d, e_d, nocc, generalized, gpa);
        defer gpa.free(corr_energies_d);

        for (2..order + 1) |k| {
            grads[k].ptr(i / 3, i % 3).* = corr_energies_d[k].der;
        }
    }

    return grads;
}

/// Computes Møller-Plesset energy corrections up to the specified perturbation order.
fn mp(comptime T: type, order: usize, g: Tensor(T, 4), C: Matrix(T), e: Vector(T), nocc: usize, generalized: bool, gpa: Allocator) ![]T {
    const nsp, const nel = if (generalized) .{ C.shape[0], nocc } else .{ 2 * C.shape[0], 2 * nocc };

    var C_SO = if (generalized) C else try Matrix(T).initZero(nsp, nsp, gpa);
    defer if (!generalized) C_SO.deinit(gpa);

    if (!generalized) {
        ao2so_coef(T, &C_SO, C);
    }

    var g_SO = if (generalized) g else try Tensor(T, 4).initZero(.{ nsp, nsp, nsp, nsp }, gpa);
    defer if (!generalized) g_SO.deinit(gpa);

    if (!generalized) {
        ao2so_pppp(T, &g_SO, g);
    }

    var g_MS = try Tensor(T, 4).init(.{ nsp, nsp, nsp, nsp }, gpa);
    defer g_MS.deinit(gpa);

    try ao2mo_pppp(T, &g_MS, g_SO, C_SO, gpa);

    var H_MS = try Matrix(T).initZero(nsp, nsp, gpa);
    defer H_MS.deinit(gpa);

    var e_SO = try Vector(T).init(nsp, gpa);
    defer e_SO.deinit(gpa);

    if (generalized) for (0..nsp) |p| {
        e_SO.ptr(p).* = e.at(p);
    };

    if (!generalized) for (0..C.shape[0]) |i| {
        e_SO.ptr(2 * i + 0).* = e.at(i);
        e_SO.ptr(2 * i + 1).* = e.at(i);
    };

    for (0..nsp) |p| for (0..nsp) |q| {
        var sum = Value(T).fromFloat(0);

        for (0..nel) |i| {
            const term1 = Value(T).init(g_MS.at(.{ p, i, q, i }));
            const term2 = Value(T).init(g_MS.at(.{ p, i, i, q }));

            sum = sum.add(term1.sub(term2));
        }

        H_MS.ptr(p, q).* = (if (p == q) Value(T).init(e_SO.at(p)) else Value(T).fromFloat(0)).sub(sum).val;
    };

    var excitations = try gpa.alloc(u32, @min(order, @min(nel, nsp - nel)));
    defer gpa.free(excitations);

    for (0..excitations.len) |i| {
        excitations[i] = @intCast(i + 1);
    }

    var dets = try generateDets(nel, nsp, excitations, gpa);

    defer {
        for (0..dets.items.len) |i| gpa.free(dets.items[i]);

        dets.deinit(gpa);
    }

    var H_CI = try Matrix(T).init(dets.items.len, dets.items.len, gpa);
    defer H_CI.deinit(gpa);

    for (0..dets.items.len) |i| for (i..dets.items.len) |j| {
        const val = slater(T, dets.items[i], dets.items[j], H_MS, g_MS);

        H_CI.ptr(i, j).* = val;
        H_CI.ptr(j, i).* = val;
    };

    var E0 = try gpa.alloc(Value(T), dets.items.len);
    defer gpa.free(E0);

    for (0..dets.items.len) |i| {
        var sum = Value(T).fromFloat(0);

        for (dets.items[i]) |p| {
            sum = sum.add(Value(T).init(e_SO.at(p)));
        }

        E0[i] = sum;
    }

    var C_coeff = try gpa.alloc(Value(T), (order + 1) * dets.items.len);
    defer gpa.free(C_coeff);

    @memset(C_coeff, Value(T).fromFloat(0));

    var E = try gpa.alloc(Value(T), order + 1);
    defer gpa.free(E);

    E[0], C_coeff[0] = .{ E0[0], Value(T).fromFloat(1) };

    for (1..order + 1) |k| {
        if (k == 1) {
            E[1] = Value(T).init(H_CI.at(0, 0)).sub(E0[0]);
        }

        if (k > 1) {
            var sum = Value(T).fromFloat(0);

            for (1..dets.items.len) |j| {
                sum = sum.add(Value(T).init(H_CI.at(0, j)).mul(C_coeff[(k - 1) * dets.items.len + j]));
            }

            E[k] = sum;
        }

        if (k == order) break;

        for (1..dets.items.len) |i| {
            var term1 = Value(T).fromFloat(0);
            var term3 = Value(T).fromFloat(0);

            for (0..dets.items.len) |j| {
                term1 = term1.add(Value(T).init(H_CI.at(i, j)).mul(C_coeff[(k - 1) * dets.items.len + j]));
            }

            const term2 = E0[i].mul(C_coeff[(k - 1) * dets.items.len + i]);

            for (1..k + 1) |j| {
                term3 = term3.add(E[j].mul(C_coeff[(k - j) * dets.items.len + i]));
            }

            C_coeff[k * dets.items.len + i] = term1.sub(term2).sub(term3).div(E0[0].sub(E0[i]));
        }
    }

    const energies = try gpa.alloc(T, order + 1);
    errdefer gpa.free(energies);

    for (0..order + 1) |k| {
        energies[k] = E[k].val;
    }

    return energies;
}

/// Computes the nuclear Hessian of the MP energy numerically and performs harmonic frequency analysis.
fn handleHessianAndFrequencies(comptime T: type, io: std.Io, opt: Options, runFn: anytype, sys: *MolecularSystem(T), log: bool, gpa: Allocator) ![]Matrix(T) {
    var hess = try gpa.alloc(Matrix(T), if (opt.hessian) |_| 1 else 0);
    errdefer if (opt.hessian) |_| gpa.free(hess);

    if (opt.hessian) |hessopt| switch (hessopt) {
        .numeric => hess[0] = try calculateNumericalHessian(T, io, runFn, opt, sys, log, gpa),
    };

    errdefer if (opt.hessian) |_| hess[0].deinit(gpa);

    if (log and opt.hessian != null) {
        var freqs = try calculateHarmonicFrequencies(T, hess[0], sys.*, gpa);
        defer freqs.deinit(gpa);

        const method_str = try std.fmt.allocPrint(gpa, "MP{d} NUMERIC", .{opt.order});
        defer gpa.free(method_str);

        try printHarmonicFrequencies(T, io, freqs, method_str);
    }

    return hess;
}
