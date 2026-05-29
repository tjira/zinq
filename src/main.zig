const std = @import("std");

const quantum_dynamics = @import("quantum_dynamics.zig");
const writer = @import("writer.zig");

const printf = writer.printf;

pub fn parse(comptime T: type, io: std.Io, fname: []const u8, gpa: std.mem.Allocator) !std.json.Parsed(T) {
    const fcontent = try std.Io.Dir.cwd().readFileAlloc(io, fname, gpa, .unlimited);
    defer gpa.free(fcontent);

    return try std.json.parseFromSlice(T, gpa, fcontent, .{});
}

pub fn run(comptime T: type, io: std.Io, fname: []const u8, gpa: std.mem.Allocator) !void {
    const options = try parse(quantum_dynamics.Options, io, fname, gpa);
    defer options.deinit();

    try quantum_dynamics.run(T, io, options.value, gpa);
}

pub fn main(init: std.process.Init) !void {
    const gpa = init.gpa;
    const io = init.io;

    var timer = std.Io.Timestamp.now(io, .real);

    try std.Io.File.stdout().writeStreamingAll(io, "ZINQ\n");

    const args = try init.minimal.args.toSlice(init.arena.allocator());

    if (args.len == 1) {
        try run(f64, io, "input.json", gpa);
    }

    for (args[1..]) |arg| {
        try run(f64, io, arg, gpa);
    }

    try printf(io, "\nTOTAL EXECUTION TIME: {f}\n", .{timer.untilNow(io, .real)});
}
