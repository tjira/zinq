const std = @import("std");

const Allocator = std.mem.Allocator;

// zig fmt: off
const Matrix           = @import("tensor.zig"   ).   Matrix;
const Potential        = @import("potential.zig").Potential;
const PotentialOptions = @import("potential.zig").  Options;
const Vector           = @import("tensor.zig"   )   .Vector;
// zig fmt: on

// zig fmt: off
const eighBatch         = @import("openblas.zig"  )        .eighBatch;
const printf            = @import("read_write.zig")           .printf;
const writeMatrixHjoin  = @import("read_write.zig") .writeMatrixHjoin;
const writeMatrixLspace = @import("read_write.zig").writeMatrixLspace;
// zig fmt: on

// OPTIONS =============================================================================================================

// zig fmt: off
const InitialConditions = struct {
    position: []const f64,
    momentum: []const f64,
    gamma:    []const f64,

    state: u32 = 0,
    seed:  u32 = 1,
};
// zig fmt: on

// zig fmt: off
const Write = struct {
    kinetic_energy:   ?[]const u8 = null,
    momentum:         ?[]const u8 = null,
    population:       ?[]const u8 = null,
    position:         ?[]const u8 = null,
    potential_energy: ?[]const u8 = null,
    total_energy:     ?[]const u8 = null,
};
// zig fmt: on

// zig fmt: off
pub const Options = struct {
    initial_conditions: InitialConditions,
    potential:           PotentialOptions,
    time_step:                        f64,
    iterations:                       u32,
    mass:                             f64,

    write:        Write     =   .{},
    adiabatic:    bool      =  true,
    log_interval: u32       =     1,
};
// zig fmt: on

// RUN =================================================================================================================

pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator, arena: Allocator) !void {
    if (log) try std.Io.File.stdout().writeStreamingAll(io, "\nCLASSICAL DYNAMICS INIT: ");

    var timer = std.Io.Timestamp.now(io, .real);

    if (log) try printf(io, "{f}\n", .{timer.untilNow(io, .real)});

    _ = T;
    _ = opt;
    _ = gpa;
    _ = arena;
}
