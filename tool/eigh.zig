const std = @import("std");
const zinq = @import("zinq");

const RealMatrix = zinq.real_matrix.RealMatrix;

const exportRealMatrix = zinq.device_write.exportRealMatrix;
const readRealMatrix = zinq.device_read.readRealMatrix;
const eigensystemHermitian = zinq.eigenproblem_solver.eigensystemHermitian;
const print = zinq.device_write.print;

pub fn help() !void {
    try print("USAGE: zinq-eigh [-i INPUT] [-o EIGENVALUE_OUTPUT EIGENVECTOR_OUTPUT] [-h]\n", .{});
}

pub fn parse(input: *[]const u8, eigenvalue_output: *[]const u8, eigenvector_output: *[]const u8, allocator: std.mem.Allocator, h: *bool) !std.process.ArgIterator {
    var argc: usize = 0; var argv = try std.process.argsWithAllocator(allocator); _ = argv.next();

    while (argv.next()) |arg| {

        if (std.mem.eql(u8, arg, "-h") or std.mem.eql(u8, arg, "--help")) {h.* = true; try help(); return argv;}

        if (std.mem.eql(u8, arg, "-i") or std.mem.eql(u8, arg, "--input")) {input.* = argv.next() orelse return error.InvalidArgument; argc += 1;}

        if (std.mem.eql(u8, arg, "-o") or std.mem.eql(u8, arg, "--output")) {
            eigenvalue_output.*  = argv.next() orelse return error.InvalidArgument; argc += 1;
            eigenvector_output.* = argv.next() orelse return error.InvalidArgument; argc += 1;
        }

        argc += 1;
    }

    return argv;
}

pub fn main() !void {
    var input: []const u8 = "A.mat"; var eigenvalue_output: []const u8 = "J.mat"; var eigenvector_output: []const u8 = "C.mat";

    var timer_total = try std.time.Timer.start(); var h = false;

    var gpa = std.heap.GeneralPurposeAllocator(.{}){}; const allocator = gpa.allocator();

    defer {
        if (gpa.deinit() == .leak) std.log.err("MEMORY LEAK DETECTED IN THE ALLOCATOR\n", .{});
    }

    var argv = try parse(&input, &eigenvalue_output, &eigenvector_output, allocator, &h); defer argv.deinit(); if (h) return;

    try print("STARTING EIGENPROBLEM SOLVER WITH INPUT '{s}' and OUTPUTS '{s}' (EIGENVALUES) and '{s}' (EIGENVECTORS)\n", .{input, eigenvalue_output, eigenvector_output});

    {
        var A = try readRealMatrix(f64, input, allocator); defer A.deinit(allocator);

        var J = try RealMatrix(f64).init(A.rows, A.cols, allocator); defer J.deinit(allocator);
        var C = try RealMatrix(f64).init(A.rows, A.cols, allocator); defer C.deinit(allocator);

        var timer_eigh = try std.time.Timer.start(); try print("\nSOLVING THE EIGENPROBLEM: ", .{});

        try eigensystemHermitian(f64, &J, &C, A);

        try print("{D}\n", .{timer_eigh.read()});

        try exportRealMatrix(f64, eigenvalue_output,  J);
        try exportRealMatrix(f64, eigenvector_output, C);
    }

    try print("\nTOTAL EXECUTION TIME: {D}\n", .{timer_total.read()});
}
