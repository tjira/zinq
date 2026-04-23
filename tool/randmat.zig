const std = @import("std");
const zinq = @import("zinq");

const RealMatrix = zinq.real_matrix.RealMatrix;

const exportRealMatrix = zinq.device_write.exportRealMatrix;
const print = zinq.device_write.print;

pub fn help() !void {
    try print("USAGE: zinq-randmat [M] [N] [-s SEED] [-o OUTPUT] [-h]\n", .{});
}

pub fn parse(m: *usize, n: *usize, seed: *usize, output: *[]const u8, symmetric: *bool, allocator: std.mem.Allocator, h: *bool) !std.process.ArgIterator {
    var argc: usize = 0;
    var argv = try std.process.argsWithAllocator(allocator);
    _ = argv.next();

    while (argv.next()) |arg| {
        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {
            h.* = true;
            try help();
            return argv;
        }

        if (argc == 0) m.* = try std.fmt.parseInt(usize, arg, 10);
        if (argc == 1) n.* = try std.fmt.parseInt(usize, arg, 10);

        if (std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--output")) {
            output.* = argv.next() orelse return error.InvalidArgument;
            argc += 1;
        }
        if (std.mem.eql(u8, arg, "-s") or std.mem.eql(u8, arg, "--seed")) {
            seed.* = try std.fmt.parseInt(usize, argv.next() orelse return error.InvalidArgument, 10);
            argc += 1;
        }

        if (std.mem.eql(u8, arg, "--symmetric")) symmetric.* = true;

        argc += 1;
    }

    return argv;
}

pub fn main() !void {
    var m: usize = 0;
    var n: usize = 0;
    var seed: usize = @intCast(std.time.milliTimestamp());
    var output: []const u8 = "A.mat";
    var symmetric: bool = false;

    var timer_total = try std.time.Timer.start();
    var h = false;

    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const allocator = gpa.allocator();

    defer {
        if (gpa.deinit() == .leak) std.log.err("MEMORY LEAK DETECTED IN THE ALLOCATOR\n", .{});
    }

    var argv = try parse(&m, &n, &seed, &output, &symmetric, allocator, &h);
    defer argv.deinit();
    if (h) return;

    try print("RANDOM {s}MATRIX GENERATION - DIM: {d}x{d}, SEED: {d}", .{ if (symmetric) "SYMMETRIC " else "", m, n, seed });

    if (output.len > 0) try print(", OUTPUT: {s}", .{output});
    try print("\n", .{});

    if (m == 0 or n == 0) {
        std.log.err("YOU NEED TO PROVIDE VALID DIMENSIONS FOR THE MATRIX\n", .{});

        return error.InvalidArgument;
    }

    {
        var A = try RealMatrix(f64).init(m, n, allocator);
        defer A.deinit(allocator);

        var timer_randomize = try std.time.Timer.start();
        try print("\nRANDOMIZING THE MATRIX: ", .{});

        A.randomize(seed);
        if (symmetric) try A.symmetrize();

        try print("{D}\n", .{timer_randomize.read()});

        if (output.len > 0) try exportRealMatrix(f64, output, A);
    }

    try print("\nTOTAL EXECUTION TIME: {D}\n", .{timer_total.read()});
}
