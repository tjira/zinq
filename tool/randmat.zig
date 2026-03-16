const std = @import("std");
const zinq = @import("zinq");

const RealMatrix = zinq.real_matrix.RealMatrix;

const exportRealMatrix = zinq.device_write.exportRealMatrix;
const print = zinq.device_write.print;

pub fn help() !void {
    try print("USAGE: randmat [M] [N] [-s SEED] [-o OUTPUT]\n", .{});
}

pub fn parse(m: *usize, n: *usize, seed: *usize, output: *[]const u8, symmetric: *bool, allocator: std.mem.Allocator) !void {
    var argc: usize = 0; var argv = try std.process.argsWithAllocator(allocator); defer argv.deinit(); _ = argv.next();

    while (argv.next()) |arg| {

        if (argc == 0) m.* = try std.fmt.parseInt(usize, arg, 10);
        if (argc == 1) n.* = try std.fmt.parseInt(usize, arg, 10);

        if (std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--output")) {output.* = argv.next() orelse return error.InvalidArgument; argc += 1;}
        if (std.mem.eql(u8, arg, "-s") or std.mem.eql(u8, arg, "--seed"  )) {seed.*   = try std.fmt.parseInt(usize, argv.next() orelse return error.InvalidArgument, 10); argc += 1;}

        if (std.mem.eql(u8, arg, "--symmetric")) symmetric.* = true;

        argc += 1;
    }
}

pub fn main() !void {
    errdefer help() catch {};

    var m: usize = 0; var n: usize = 0; var seed: usize = @intCast(std.time.milliTimestamp()); var output: []const u8 = "A.mat"; var symmetric: bool = false;

    var timer_total = try std.time.Timer.start();

    var gpa = std.heap.GeneralPurposeAllocator(.{}){}; const allocator = gpa.allocator();

    defer {
        if (gpa.deinit() == .leak) std.log.err("MEMORY LEAK DETECTED IN THE ALLOCATOR\n", .{});
    }

    try parse(&m, &n, &seed, &output, &symmetric, allocator);

    if (m == 0 or n == 0 or (m != n and symmetric)) return error.InvalidArgument;

    {
        var timer_alloc = try std.time.Timer.start(); try print("ALLOCATING {d}x{d} REAL MATRIX: ", .{m, n});

        var A = try RealMatrix(f64).init(m, n, allocator); defer A.deinit(allocator);

        try print("{D}\n", .{timer_alloc.read()});

        var timer_randomize = try std.time.Timer.start(); try print("RANDOMIZING REAL MATRIX WITH SEED {d}: ", .{seed});

        A.randomize(seed);

        try print("{D}\n", .{timer_randomize.read()});

        if (symmetric) {

            var timer_symmetrize = try std.time.Timer.start(); try print("SYMMETRIZING REAL MATRIX: ", .{});

            try A.symmetrize();

            try print("{D}\n", .{timer_symmetrize.read()});
        }

        var timer_export = try std.time.Timer.start(); try print("EXPORTING REAL MATRIX TO '{s}': ", .{output});

        try exportRealMatrix(f64, output, A);

        try print("{D}\n", .{timer_export.read()});
    }

    try print("TOTAL EXECUTION TIME: {D}\n", .{timer_total.read()});
}
