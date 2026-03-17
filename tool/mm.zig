const std = @import("std");
const zinq = @import("zinq");

const RealMatrix = zinq.real_matrix.RealMatrix;

const exportRealMatrix = zinq.device_write.exportRealMatrix;
const readRealMatrix = zinq.device_read.readRealMatrix;
const mm = zinq.matrix_multiplication.mm;
const print = zinq.device_write.print;

pub fn help() !void {
    try print("USAGE: zinq-mm [-i FIRST_MATRIX SECOND_MATRIX] [-o RESULT_MATRIX] [-h]\n", .{});
}

pub fn parse(first: *[]const u8, second: *[]const u8, result: *[]const u8, allocator: std.mem.Allocator, h: *bool) !std.process.ArgIterator {
    var argc: usize = 0; var argv = try std.process.argsWithAllocator(allocator); _ = argv.next();

    while (argv.next()) |arg| {

        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {h.* = true; try help(); return argv;}

        if (std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--output")) {result.* = argv.next() orelse return error.InvalidArgument; argc += 1;}

        if (std.mem.eql(u8, arg, "-i") or std.mem.eql(u8, arg, "--input")) {
            first.*  = argv.next() orelse return error.InvalidArgument; argc += 1;
            second.* = argv.next() orelse return error.InvalidArgument; argc += 1;
        }

        argc += 1;
    }

    return argv;
}

pub fn main() !void {
    var first: []const u8 = "A.mat"; var second: []const u8 = "B.mat"; var result: []const u8 = "C.mat";

    var timer_total = try std.time.Timer.start(); var h = false;

    var gpa = std.heap.GeneralPurposeAllocator(.{}){}; const allocator = gpa.allocator();

    defer {
        if (gpa.deinit() == .leak) std.log.err("MEMORY LEAK DETECTED IN THE ALLOCATOR\n", .{});
    }

    var argv = try parse(&first, &second, &result, allocator, &h); defer argv.deinit(); if (h) return;

    try print("MATRIX MULTIPLICATION - INPUTS: {s}/{s}", .{first, second});

    if (result.len > 0) try print(", OUTPUT: {s}", .{result}); try print("\n", .{});

    {
        var A = try readRealMatrix(f64, first,  allocator); defer A.deinit(allocator);
        var B = try readRealMatrix(f64, second, allocator); defer B.deinit(allocator);

        var C = try RealMatrix(f64).init(A.rows, B.cols, allocator); defer C.deinit(allocator);

        var timer_mm = try std.time.Timer.start(); try print("\nMULTIPLYING THE MATRICES: ", .{});

        try mm(f64, &C, A, false, B, false);

        try print("{D}\n", .{timer_mm.read()});

        if (result.len > 0) try exportRealMatrix(f64, result, C);
    }

    try print("\nTOTAL EXECUTION TIME: {D}\n", .{timer_total.read()});
}
