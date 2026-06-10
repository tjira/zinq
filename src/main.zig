const std = @import("std");

const Allocator = std.mem.Allocator;

pub const ClassicalDynamicsOptions = @import("classical_dynamics.zig").Options;
pub const classical_dynamics_run = @import("classical_dynamics.zig").run;
pub const QuantumDynamicsOptions = @import("quantum_dynamics.zig").Options;
pub const quantum_dynamics_run = @import("quantum_dynamics.zig").run;
pub const SurfaceHopping = @import("surface_hopping.zig").SurfaceHopping;
pub const SurfaceHoppingOptions = @import("surface_hopping.zig").Options;

pub const PotentialOptions = @import("potential.zig").Options;

const printf = @import("read_write.zig").printf;

// OPTIONS =============================================================================================================

const Options = struct {
    zinq: []union(enum) {
        classical_dynamics: ClassicalDynamicsOptions,
        quantum_dynamics: QuantumDynamicsOptions,
    },
};

// INIT FUNCTIONS ======================================================================================================

fn parse(comptime T: type, io: std.Io, fname: []const u8, arena: Allocator) !std.json.Parsed(T) {
    const fcontent = try std.Io.Dir.cwd().readFileAlloc(io, fname, arena, .unlimited);

    return try std.json.parseFromSlice(T, arena, fcontent, .{});
}

fn run(comptime T: type, io: std.Io, fname: []const u8, gpa: Allocator, arena: Allocator) !void {
    const inputs = (try parse(Options, io, fname, arena)).value.zinq;

    for (inputs, 0..) |e, i| {
        try printf(io, "/#{d}\n", .{i + 1});

        switch (e) {
            .classical_dynamics => |field| _ = try classical_dynamics_run(T, io, field, true, gpa, arena),
            .quantum_dynamics => |field| _ = try quantum_dynamics_run(T, io, field, true, gpa, arena),
        }
    }
}

fn targets(args: []const []const u8) []const []const u8 {
    return if (args.len > 1) args[1..] else &[_][]const u8{"input.json"};
}

pub fn main(init: std.process.Init) !void {
    var timer = std.Io.Timestamp.now(init.io, .real);

    try std.Io.File.stdout().writeStreamingAll(init.io, "ZINQ\n");

    for (targets(try init.minimal.args.toSlice(init.arena.allocator()))) |e| {
        try printf(init.io, "\nRUNNING TARGET: {s}", .{e});

        try run(f64, init.io, e, init.gpa, init.arena.allocator());
    }

    try printf(init.io, "\nTOTAL EXECUTION TIME: {f}\n", .{timer.untilNow(init.io, .real)});
}
