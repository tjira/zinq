const std = @import("std");
const zinq = @import("zinq");

const RealMatrix = zinq.real_matrix.RealMatrix;

const exportRealMatrix = zinq.device_write.exportRealMatrix;
const readRealMatrix = zinq.device_read.readRealMatrix;
const print = zinq.device_write.print;

pub fn help() !void {
    try print("USAGE: zinq-mtrans [-i INPUT] [-o OUTPUT] [-h]\n", .{});
}

pub fn parse(input: *[]const u8, result: *[]const u8, allocator: std.mem.Allocator, h: *bool) !std.process.ArgIterator {
    var argc: usize = 0;
    var argv = try std.process.argsWithAllocator(allocator);
    _ = argv.next();

    while (argv.next()) |arg| {
        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            h.* = true;
            try help();
            return argv;
        }

        if (std.mem.eql(u8, arg, "-i") or std.mem.eql(u8, arg, "--input")) {
            input.* = argv.next() orelse return error.InvalidArgument;
            argc += 1;
        }
        if (std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--output")) {
            result.* = argv.next() orelse return error.InvalidArgument;
            argc += 1;
        }

        argc += 1;
    }

    return argv;
}

pub fn main() !void {
    var input: []const u8 = "A.mat";
    var result: []const u8 = "B.mat";

    var timer_total = try std.time.Timer.start();
    var h = false;

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();

    defer {
        if (gpa.deinit() == .leak) std.log.err("MEMORY LEAK DETECTED IN THE ALLOCATOR\n", .{});
    }

    var argv = try parse(&input, &result, allocator, &h);
    defer argv.deinit();
    if (h) return;

    try print("MATRIX TRANSPOSE - INPUT: {s}", .{input});

    if (result.len > 0) try print(", OUTPUT: {s}", .{result});
    try print("\n", .{});

    {
        var A = try readRealMatrix(f64, input, allocator);
        defer A.deinit(allocator);

        var timer_mtrans = try std.time.Timer.start();
        try print("\nTRANSPOSING THE MATRIX: ", .{});

        try A.transpose();

        try print("{D}\n", .{timer_mtrans.read()});

        if (result.len > 0) try exportRealMatrix(f64, result, A);
    }

    try print("\nTOTAL EXECUTION TIME: {D}\n", .{timer_total.read()});
}
