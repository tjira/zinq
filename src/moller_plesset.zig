const std = @import("std");

const Allocator = std.mem.Allocator;

const HartreeFockOptions = @import("hartree_fock.zig").Options;
const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

const hartree_fock_run = @import("hartree_fock.zig").run;
const printf = @import("read_write.zig").printf;
const ao2mo_oovv = @import("integral_transform.zig").ao2mo_oovv;

// OPTIONS =============================================================================================================

pub const Options = struct {
    hartree_fock: HartreeFockOptions,

    order: u32 = 2,
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

// RESULT STRUCT =======================================================================================================

pub fn Result(comptime T: type) type {
    return struct {
        hartree_fock: HartreeFockResult(T),

        correlation_energy: T,

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.hartree_fock.deinit(gpa);
        }
    };
}

// RUN =================================================================================================================

pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator) !Result(T) {
    if (opt.order < 2) {
        @panic("MOLLER-PLESSET ORDER LOWER THAN 2 IS RIDICULOUS");
    }

    if (opt.order > 2) {
        @panic("MOLLER-PLESSET ORDER HIGHER THAN 2 IS NOT SUPPORTED");
    }

    var hfres = try hartree_fock_run(T, io, opt.hartree_fock, log, gpa);
    errdefer hfres.deinit(gpa);

    const corr_energy = switch(opt.order) {
        2 => try mp2(T, hfres, opt.hartree_fock.generalized, gpa),
        else => unreachable,
    };

    if (log) try printf(io, "\nMP{d} CORRELATION ENERGY: {d:.14}\n", .{ opt.order, corr_energy });

    if (log) try printf(io, "\nFINAL MP{d} ENERGY: {d:.14}\n", .{ opt.order, hfres.energy + corr_energy });

    return Result(T){ .hartree_fock = hfres, .correlation_energy = corr_energy };
}

