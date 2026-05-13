const std = @import("std");
const zinq = @import("zinq");

const RealMatrix = zinq.real_matrix.RealMatrix;

const exportRealMatrix = zinq.device_write.exportRealMatrix;
const print = zinq.device_write.print;

pub fn help(io: std.Io) !void {
    try print(io, "USAGE: zinq-randmat [M] [N] [-s SEED] [-o OUTPUT] [-h]\n", .{});
}

pub fn parse(io: std.Io, m: *usize, n: *usize, seed: *usize, output: *[]const u8, symmetric: *bool, args: []const []const u8, h: *bool) !void {
    var argc: usize = 0;
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
            output.* = if (i < args.len) args[i] else return error.InvalidArgument;
            continue;
        }
        if (std.mem.eql(u8, arg, "-s") or std.mem.eql(u8, arg, "--seed")) {
            i += 1;
            seed.* = try std.fmt.parseInt(usize, if (i < args.len) args[i] else return error.InvalidArgument, 10);
            continue;
        }

        if (std.mem.eql(u8, arg, "--symmetric")) {
            symmetric.* = true;
            continue;
        }

        if (argc == 0) m.* = try std.fmt.parseInt(usize, arg, 10);
        if (argc == 1) n.* = try std.fmt.parseInt(usize, arg, 10);

        argc += 1;
    }
}

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    const io = init.io;

    var m: usize = 0;
    var n: usize = 0;
    var seed: usize = @intCast(std.Io.Timestamp.now(io, .real).toMilliseconds());
    var output: []const u8 = "A.mat";
    var symmetric: bool = false;

    var timer_total = std.Io.Timestamp.now(io, .real);
    var h = false;

    const args = try init.minimal.args.toSlice(init.arena.allocator());

    try parse(io, &m, &n, &seed, &output, &symmetric, args[1..], &h);

    if (h) return;

    try print(io, "RANDOM {s}MATRIX GENERATION - DIM: {d}x{d}, SEED: {d}", .{ if (symmetric) "SYMMETRIC " else "", m, n, seed });

    if (output.len > 0) try print(io, ", OUTPUT: {s}", .{output});
    try print(io, "\n", .{});

    if (m == 0 or n == 0) {
        try print(io, "YOU NEED TO PROVIDE VALID DIMENSIONS FOR THE MATRIX\n", .{});

        return error.InvalidArgument;
    }

    {
        var A = try RealMatrix(f64).init(m, n, allocator);
        defer A.deinit(allocator);

        var timer_randomize = std.Io.Timestamp.now(io, .real);
        try print(io, "\nRANDOMIZING THE MATRIX: ", .{});

        A.randomize(seed);
        if (symmetric) try A.symmetrize();

        try print(io, "{f}\n", .{timer_randomize.untilNow(io, .real)});

        if (output.len > 0) try exportRealMatrix(f64, io, output, A);
    }

    try print(io, "\nTOTAL EXECUTION TIME: {f}\n", .{timer_total.untilNow(io, .real)});
}
