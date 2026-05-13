const std = @import("std");
const zinq = @import("zinq");

const RealMatrix = zinq.real_matrix.RealMatrix;

const exportRealMatrix = zinq.device_write.exportRealMatrix;
const readRealMatrix = zinq.device_read.readRealMatrix;
const mm = zinq.matrix_multiplication.mm;
const print = zinq.device_write.print;

pub fn help(io: std.Io) !void {
    try print(io, "USAGE: zinq-mm [-i FIRST_MATRIX SECOND_MATRIX] [-o RESULT_MATRIX] [-h]\n", .{});
}

pub fn parse(io: std.Io, first: *[]const u8, second: *[]const u8, result: *[]const u8, args: []const []const u8, h: *bool) !void {
    var i: usize = 0;

    while (i < args.len) : (i += 1) {
        const arg = args[i];

        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            h.* = true;
            try help(io);
            return;
        }

        if (std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--output")) {
            i += 1;
            result.* = if (i < args.len) args[i] else return error.InvalidArgument;
        }

        if (std.mem.eql(u8, arg, "-i") or std.mem.eql(u8, arg, "--input")) {
            i += 1;
            first.* = if (i < args.len) args[i] else return error.InvalidArgument;
            i += 1;
            second.* = if (i < args.len) args[i] else return error.InvalidArgument;
        }
    }
}

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    const io = init.io;

    var first: []const u8 = "A.mat";
    var second: []const u8 = "B.mat";
    var result: []const u8 = "C.mat";

    var timer_total = std.Io.Timestamp.now(io, .real);
    var h = false;

    const args = try init.minimal.args.toSlice(init.arena.allocator());

    try parse(io, &first, &second, &result, args[1..], &h);

    if (h) return;

    try print(io, "MATRIX MULTIPLICATION - INPUTS: {s}/{s}", .{ first, second });

    if (result.len > 0) try print(io, ", OUTPUT: {s}", .{ result });
    try print(io, "\n", .{});

    {
        var A = try readRealMatrix(f64, io, first, allocator);
        defer A.deinit(allocator);

        var B = try readRealMatrix(f64, io, second, allocator);
        defer B.deinit(allocator);

        var C = try RealMatrix(f64).init(A.rows, B.cols, allocator);
        defer C.deinit(allocator);

        var timer_mm = std.Io.Timestamp.now(io, .real);
        try print(io, "\nMULTIPLYING THE MATRICES: ", .{});

        try mm(f64, &C, A, false, B, false);

        try print(io, "{f}\n", .{timer_mm.untilNow(io, .real)});

        if (result.len > 0) try exportRealMatrix(f64, io, result, C);
    }

    try print(io, "\nTOTAL EXECUTION TIME: {f}\n", .{timer_total.untilNow(io, .real)});
}
