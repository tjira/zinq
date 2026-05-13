const std = @import("std");
const zinq = @import("zinq");

const RealMatrix = zinq.real_matrix.RealMatrix;

const exportRealMatrix = zinq.device_write.exportRealMatrix;
const readRealMatrix = zinq.device_read.readRealMatrix;
const print = zinq.device_write.print;

pub fn help(io: std.Io) !void {
    try print(io, "USAGE: zinq-mtrans [-i INPUT] [-o OUTPUT] [-h]\n", .{});
}

pub fn parse(io: std.Io, input: *[]const u8, result: *[]const u8, args: []const []const u8, h: *bool) !void {
    var i: usize = 0;

    while (i < args.len) : (i += 1) {
        const arg = args[i];

        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            h.* = true;
            try help(io);
            return;
        }

        if (std.mem.eql(u8, arg, "-i") or std.mem.eql(u8, arg, "--input")) {
            i += 1;
            input.* = if (i < args.len) args[i] else return error.InvalidArgument;
        }
        if (std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--output")) {
            i += 1;
            result.* = if (i < args.len) args[i] else return error.InvalidArgument;
        }
    }
}

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    const io = init.io;

    var input: []const u8 = "A.mat";
    var result: []const u8 = "B.mat";

    var timer_total = std.Io.Timestamp.now(io, .real);
    var h = false;

    const args = try init.minimal.args.toSlice(init.arena.allocator());

    try parse(io, &input, &result, args[1..], &h);

    if (h) return;

    try print(io, "MATRIX TRANSPOSE - INPUT: {s}", .{input});

    if (result.len > 0) try print(io, ", OUTPUT: {s}", .{result});
    try print(io, "\n", .{});

    {
        var A = try readRealMatrix(f64, io, input, allocator);
        defer A.deinit(allocator);

        var timer_mtrans = std.Io.Timestamp.now(io, .real);
        try print(io, "\nTRANSPOSING THE MATRIX: ", .{});

        try A.transpose();

        try print(io, "{f}\n", .{timer_mtrans.untilNow(io, .real)});

        if (result.len > 0) try exportRealMatrix(f64, io, result, A);
    }

    try print(io, "\nTOTAL EXECUTION TIME: {f}\n", .{timer_total.untilNow(io, .real)});
}
