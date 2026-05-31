const std = @import("std");

const Allocator = std.mem.Allocator;

const QuantumDynamicsOptions = @import("quantum_dynamics.zig").Options;

const printf = @import("read_write.zig").printf;
const quantum_dynamics_run = @import("quantum_dynamics.zig").run;

// OPTIONS =============================================================================================================

const Options = struct {
    zinq: []union(enum) {
        quantum_dynamics: QuantumDynamicsOptions,
    },
};

// INIT FUNCTIONS ======================================================================================================

fn parse(comptime T: type, io: std.Io, fname: []const u8, arena: Allocator) !std.json.Parsed(T) {
    const fcontent = try std.Io.Dir.cwd().readFileAlloc(io, fname, arena, .unlimited);

    return try std.json.parseFromSlice(T, arena, fcontent, .{});
}

fn run(comptime T: type, io: std.Io, fname: []const u8, gpa: Allocator, arena: Allocator) !void {
    for ((try parse(Options, io, fname, arena)).value.zinq) |opt| switch (opt) {
        .quantum_dynamics => |field| try quantum_dynamics_run(T, io, field, gpa, arena),
    };
}

fn targets(args: []const []const u8) []const []const u8 {
    return if (args.len > 1) args[1..] else &[_][]const u8{"input.json"};
}

pub fn main(init: std.process.Init) !void {
    var timer = std.Io.Timestamp.now(init.io, .real);

    try std.Io.File.stdout().writeStreamingAll(init.io, "ZINQ\n");

    for (targets(try init.minimal.args.toSlice(init.arena.allocator()))) |target| {
        try run(f64, init.io, target, init.gpa, init.arena.allocator());
    }

    try printf(init.io, "\nTOTAL EXECUTION TIME: {f}\n", .{timer.untilNow(init.io, .real)});
}
