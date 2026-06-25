const std = @import("std");

const Allocator = std.mem.Allocator;

const HartreeFockOptions = @import("hartree_fock.zig").Options;
const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

const ao2mo_oovv = @import("integral_transform.zig").ao2mo_oovv;
const ao2mo_pppp = @import("integral_transform.zig").ao2mo_pppp;
const gradientResponse = @import("hartree_fock.zig").gradientResponse;
const hartree_fock_run = @import("hartree_fock.zig").run;
const printf = @import("read_write.zig").printf;

const ScalarDual = @import("dual.zig").ScalarDual;
const Value = @import("value.zig").Value;
const Integrals = @import("molecular_integrals.zig").Integrals;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;

// OPTIONS =============================================================================================================

pub const Options = struct {
    hartree_fock: HartreeFockOptions,

    order: u32 = 2,

    gradient: ?union(enum) {
        analytic: struct {},
    } = null,
};

// MOLLER-PLESSET FUNCTIONS ============================================================================================

pub fn mp(comptime T: type, order: usize, g: Tensor(T, 4), C: Matrix(T), e: Vector(T), nocc: usize, generalized: bool, gpa: Allocator) !T {
    return switch (order) {
        2 => try mp2(T, g, C, e, nocc, generalized, gpa),
        else => error.MollerPlessetInvalidOrder,
    };
}

pub fn mp2(comptime T: type, g: Tensor(T, 4), C: Matrix(T), e: Vector(T), nocc: usize, generalized: bool, gpa: Allocator) !T {
    var g_oovv = try Tensor(T, 4).initZero(.{ nocc, nocc, C.shape[1] - nocc, C.shape[1] - nocc }, gpa);
    defer g_oovv.deinit(gpa);

    try ao2mo_oovv(T, &g_oovv, g, C, nocc, gpa);

    var energy = Value(T).fromFloat(0);

    for (0..nocc) |i| for (0..nocc) |j| for (0..C.shape[1] - nocc) |a| for (0..C.shape[1] - nocc) |b| {
        const e_i = Value(T).init(e.at(i));
        const e_j = Value(T).init(e.at(j));

        const e_a = Value(T).init(e.at(a + nocc));
        const e_b = Value(T).init(e.at(b + nocc));

        const denom = e_i.add(e_j).sub(e_a).sub(e_b);

        const term1 = Value(T).init(g_oovv.at(.{ i, j, a, b }));
        const term2 = Value(T).init(g_oovv.at(.{ i, j, b, a }));

        if (generalized) {
            const term = term1.sub(term2);

            energy = energy.add(term.mul(term).div(denom));
        }

        if (!generalized) {
            const term = term1.muls(2).sub(term2);

            energy = energy.add(term1.mul(term).div(denom));
        }
    };

    return (if (generalized) energy.muls(0.25) else energy).val;
}

pub fn gradient(comptime T: type, order: usize, hfres: HartreeFockResult(T), gpa: Allocator) !Matrix(T) {
    const nbf, const generalized = .{ hfres.C.shape[0], hfres.ints.sys.nbf != hfres.C.shape[0] };

    var G = try Matrix(T).initZero(hfres.ints.sys.atoms.len, 3, gpa);
    errdefer G.deinit(gpa);

    const nocc = if (generalized) hfres.ints.sys.nel else hfres.ints.sys.nel / 2;

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

        const corr_energy_d = try mp(ScalarDual(T), order, g_d, C_d, e_d, nocc, generalized, gpa);

        G.ptr(i / 3, i % 3).* = corr_energy_d.der;
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

            for (0..self.gradient.len) |i| {
                self.gradient[i].deinit(gpa);
            }

            gpa.free(self.gradient);
        }
    };
}

// RUN =================================================================================================================

pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator) !Result(T) {
    const generalized = opt.hartree_fock.generalized;

    var hf_opt = opt.hartree_fock;

    if (opt.gradient != null) {
        hf_opt.gradient = .{ .analytic = .{} };

        if (hf_opt.response == null) {
            hf_opt.response = .{};
        }
    }

    var hfres = try hartree_fock_run(T, io, hf_opt, log, gpa);
    errdefer hfres.deinit(gpa);

    var energy = try gpa.alloc(T, 1);
    errdefer gpa.free(energy);

    const grad = try gpa.alloc(Matrix(T), if (opt.gradient) |_| 1 else 0);
    errdefer if (opt.gradient) |_| gpa.free(grad);

    energy[0] = hfres.energy[0];

    if (opt.gradient) |_| {
        grad[0] = try hfres.gradient[0].clone(gpa);
    }

    errdefer {
        if (opt.gradient) |_| grad[0].deinit(gpa);
    }

    const nocc = if (generalized) hfres.ints.sys.nel else hfres.ints.sys.nel / 2;

    for (2..(opt.order + 1)) |i| {
        const step_energy = try mp(T, i, hfres.ints.g.?, hfres.C, hfres.e, nocc, generalized, gpa);

        energy[0] += step_energy;

        if (log) {
            try printf(io, "\nMP{d} TOTAL ENERGY: {d:.14} Eh\n", .{ i, energy[0] });
        }

        if (opt.gradient) |_| {
            var step_grad = try gradient(T, i, hfres, gpa);
            defer step_grad.deinit(gpa);

            for (0..grad[0].nrow()) |j| for (0..grad[0].ncol()) |k| {
                grad[0].ptr(j, k).* += step_grad.at(j, k);
            };

            if (log) {
                try printf(io, "\nMP{d} NUCLEAR ENERGY GRADIENT\n", .{i});

                for (0..grad[0].nrow()) |j| for (0..grad[0].ncol()) |k| {
                    try printf(io, "{d:20.14}{s}", .{ grad[0].at(j, k), if (k == 2) "\n" else " " });
                };
            }
        }
    }

    return Result(T){ .hartree_fock = hfres, .energy = energy, .gradient = grad };
}
