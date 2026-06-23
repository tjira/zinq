const builtin = @import("builtin");
const config = @import("config");
const std = @import("std");

const Allocator = std.mem.Allocator;

pub const ClassicalDynamicsOptions = @import("classical_dynamics.zig").Options;
pub const classical_dynamics_run = @import("classical_dynamics.zig").run;
pub const ConfigurationInteractionOptions = @import("configuration_interaction.zig").Options;
pub const configuration_interaction_run = @import("configuration_interaction.zig").run;
pub const HartreeFockOptions = @import("hartree_fock.zig").Options;
pub const hartree_fock_run = @import("hartree_fock.zig").run;
pub const MolecularIntegralsOptions = @import("molecular_integrals.zig").Options;
pub const molecular_integrals_run = @import("molecular_integrals.zig").run;
pub const MollerPlessetOptions = @import("moller_plesset.zig").Options;
pub const moller_plesset_run = @import("moller_plesset.zig").run;
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
        configuration_interaction: ConfigurationInteractionOptions,
        hartree_fock: HartreeFockOptions,
        molecular_integrals: MolecularIntegralsOptions,
        moller_plesset: MollerPlessetOptions,
        quantum_dynamics: QuantumDynamicsOptions,
    },
};

// INIT FUNCTIONS ======================================================================================================

fn parse(comptime T: type, io: std.Io, fname: []const u8, arena: Allocator) !?std.json.Parsed(T) {
    const fcontent = std.Io.Dir.cwd().readFileAlloc(io, fname, arena, .unlimited) catch |err| {
        if (err == error.FileNotFound) {
            try printf(io, "\nINPUT FILE '{s}' NOT FOUND\n", .{fname});

            return null;
        }

        return err;
    };

    return try std.json.parseFromSlice(T, arena, fcontent, .{});
}

fn run(comptime T: type, io: std.Io, fname: []const u8, gpa: Allocator, arena: Allocator) !void {
    const parsed = try parse(Options, io, fname, arena) orelse return;

    for (0..parsed.value.zinq.len) |i| {
        try printf(io, "\nRUNNING TARGET: {s}/#{d}\n", .{ fname, i + 1 });

        switch (parsed.value.zinq[i]) {
            inline else => |field, tag| {
                var result = try @field(@This(), @tagName(tag) ++ "_run")(T, io, field, true, gpa);
                defer result.deinit(gpa);
            },
        }
    }
}

fn targets(args: []const []const u8) []const []const u8 {
    return if (args.len > 1) args[1..] else &[_][]const u8{"input.json"};
}

pub fn main(init: std.process.Init) !void {
    var timer = std.Io.Timestamp.now(init.io, .real);

    const v_major = builtin.zig_version.major;
    const v_minor = builtin.zig_version.minor;
    const v_patch = builtin.zig_version.patch;

    try printf(init.io, "ZIG: v{d}.{d}.{d}, ZINQ: {s}\n", .{ v_major, v_minor, v_patch, config.version });

    for (targets(try init.minimal.args.toSlice(init.arena.allocator()))) |e| {
        try run(f64, init.io, e, init.gpa, init.arena.allocator());
    }

    try printf(init.io, "\nTOTAL EXECUTION TIME: {f}\n", .{timer.untilNow(init.io, .real)});
}
