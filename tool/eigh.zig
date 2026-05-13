const std = @import("std");
const zinq = @import("zinq");

const RealMatrix = zinq.real_matrix.RealMatrix;

const exportRealMatrix = zinq.device_write.exportRealMatrix;
const readRealMatrix = zinq.device_read.readRealMatrix;
const eigensystemHermitian = zinq.eigenproblem_solver.eigensystemHermitian;
const print = zinq.device_write.print;

pub fn help(io: std.Io) !void {
    try print(io, "USAGE: zinq-eigh [-i INPUT] [-o EIGENVALUE_OUTPUT EIGENVECTOR_OUTPUT] [-h]\n", .{});
}

pub fn parse(io: std.Io, input: *[]const u8, eigenvalue_output: *[]const u8, eigenvector_output: *[]const u8, args: []const []const u8, h: *bool) !void {
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
            eigenvalue_output.* = if (i < args.len) args[i] else return error.InvalidArgument;
            i += 1;
            eigenvector_output.* = if (i < args.len) args[i] else return error.InvalidArgument;
        }
    }
}

pub fn main(init: std.process.Init) !void {
    const allocator = init.gpa;
    const io = init.io;

    var input: []const u8 = "A.mat";
    var eigenvalue_output: []const u8 = "J.mat";
    var eigenvector_output: []const u8 = "C.mat";

    var timer_total = std.Io.Timestamp.now(io, .real);
    var h = false;

    const args = try init.minimal.args.toSlice(init.arena.allocator());

    try parse(io, &input, &eigenvalue_output, &eigenvector_output, args[1..], &h);

    if (h) return;

    try print(io, "EIGENPROBLEM SOLVER - INPUT: {s}", .{input});

    if (eigenvalue_output.len > 0) try print(io, ", EIGENVALUE OUTPUT: {s}", .{eigenvalue_output});
    if (eigenvector_output.len > 0) try print(io, ", EIGENVECTOR OUTPUT: {s}", .{eigenvector_output});

    try print(io, "\n", .{});

    {
        var A = try readRealMatrix(f64, io, input, allocator);
        defer A.deinit(allocator);

        var J = try RealMatrix(f64).init(A.rows, A.cols, allocator);
        defer J.deinit(allocator);

        var C = try RealMatrix(f64).init(A.rows, A.cols, allocator);
        defer C.deinit(allocator);

        var timer_eigh = std.Io.Timestamp.now(io, .real);
        try print(io, "\nSOLVING THE EIGENPROBLEM: ", .{});

        try eigensystemHermitian(f64, &J, &C, A);

        try print(io, "{f}\n", .{timer_eigh.untilNow(io, .real)});

        if (eigenvalue_output.len > 0) try exportRealMatrix(f64, io, eigenvalue_output, J);
        if (eigenvalue_output.len > 0) try exportRealMatrix(f64, io, eigenvector_output, C);
    }

    try print(io, "\nTOTAL EXECUTION TIME: {f}\n", .{timer_total.untilNow(io, .real)});
}
