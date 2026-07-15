//! Computes potential energy surfaces on a multidimensional grid for plotting and interpolation.

const std = @import("std");

const Allocator = std.mem.Allocator;

const Grid = @import("wavepacket.zig").Grid;
const Matrix = @import("tensor.zig").Matrix;
const Potential = @import("potential.zig").Potential;
const PotentialOptions = @import("potential.zig").Options;

const eighBatch = @import("linear_algebra.zig").eighBatch;
const printf = @import("read_write.zig").printf;
const writeMatrixHjoin = @import("read_write.zig").writeMatrixHjoin;

/// Configuration options for evaluating and exporting potential energy matrices on a grid.
pub const Options = struct {
    potential: PotentialOptions,

    adiabatic: bool = false,
    time: f64 = 0,

    grid: struct {
        bounds: []const [2]f64,
        npoint: u32,
    },

    write: Write = .{},
};

/// Result containing the evaluated potential energy matrix.
pub fn Result(comptime T: type) type {
    return struct {
        U: Matrix(T),

        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.U.deinit(gpa);
        }
    };
}

/// File paths for saving potential energy grid values.
const Write = struct {
    potential: ?[]const u8 = null,
};

/// Evaluates the potential on a grid, diagonalizes if adiabatic, and writes the output.
pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator) !Result(T) {
    var timer = std.Io.Timestamp.now(io, .real);

    var pot = try Potential(T).init(io, opt.potential, gpa);
    defer pot.deinit(gpa);

    var grid = try Grid(T).init(opt.grid.bounds, opt.grid.npoint, gpa);
    defer grid.deinit(gpa);

    const nstate, const nrow = .{ pot.nstate(), grid.r.nrow() };

    var U = try Matrix(T).initZero(nrow, nstate * nstate, gpa);
    errdefer U.deinit(gpa);

    var W = if (opt.adiabatic) try Matrix(T).init(nrow, nstate, gpa) else null;
    defer if (opt.adiabatic) W.?.deinit(gpa);

    var C = if (opt.adiabatic) try Matrix(T).init(nrow, nstate * nstate, gpa) else null;
    defer if (opt.adiabatic) C.?.deinit(gpa);

    if (log) {
        try printf(io, "\nMATRIX ALLOCATION: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    pot.evalBatch(T, &U, grid.r, opt.time);

    if (log) {
        try printf(io, "COMPUTE POTENTIAL: {f}\n", .{timer.untilNow(io, .real)});
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.adiabatic) {
        try eighBatch(T, &W.?, &C.?, U);

        U.zero();

        for (0..nrow) |i| for (0..nstate) |j| {
            U.ptr(i, j * nstate + j).* = W.?.at(i, j);
        };

        if (log) {
            try printf(io, "SOLVE EIGENVALUES: {f}\n", .{timer.untilNow(io, .real)});
        }
    }

    timer = std.Io.Timestamp.now(io, .real);

    if (opt.write.potential) |path| {
        try writeMatrixHjoin(T, io, path, grid.r, null, U, null);

        if (log) {
            try printf(io, "POTENTIAL WRITTEN: {f}\n", .{timer.untilNow(io, .real)});
        }
    }

    return Result(T){ .U = U };
}
