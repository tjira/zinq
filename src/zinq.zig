//! Main file of the program.

const builtin = @import("builtin");
const std = @import("std");

pub const classical_dynamics = @import("classical_dynamics.zig");
pub const classical_particle = @import("classical_particle.zig");
pub const complex_matrix = @import("complex_matrix.zig");
pub const device_write = @import("device_write.zig");
pub const electronic_potential = @import("electronic_potential.zig");
pub const global_variables = @import("global_variables.zig");
pub const harmonic_potential = @import("harmonic_potential.zig");
pub const linear_algebra = @import("linear_algebra.zig");
pub const math_functions = @import("math_functions.zig");
pub const real_matrix = @import("real_matrix.zig");
pub const real_vector = @import("real_vector.zig");

const Complex = std.math.complex.Complex;
const RealMatrix = real_matrix.RealMatrix;
const ComplexMatrix = complex_matrix.ComplexMatrix;

const print = device_write.print;
const exportRealMatrix = device_write.exportRealMatrix;
const printRealMatrix = device_write.printRealMatrix;
const mmRealAlloc = linear_algebra.mmRealAlloc;

const MAX_INPUT_FILE_BYTES = global_variables.MAX_INPUT_FILE_BYTES;

/// Structure to handle different targets.
const Handler = struct {
    Options: type, run: fn (comptime type, anytype, std.mem.Allocator) anyerror!void,
};

/// Array of handlers for different targets.
const handlers = [_]struct {key: []const u8, handler: Handler}{
    .{.key = "classical_dynamics", .handler = Handler{
        .Options = classical_dynamics.Options(f64),
        .run = classical_dynamics.run
    }}
};

/// Parse the input JSON file and run the corresponding target.
pub fn parse(path: []const u8, allocator: std.mem.Allocator) !void {
    const file_contents = try std.fs.cwd().readFileAlloc(allocator, path, MAX_INPUT_FILE_BYTES); defer allocator.free(file_contents);

    try print("\nPROCESSED FILE: {s}\n", .{path});

    const input_json = try std.json.parseFromSlice(std.json.Value, allocator, file_contents, .{}); defer input_json.deinit();

    for (input_json.value.object.get("zinq").?.array.items) |object| {

        const name = object.object.get("name") orelse return error.MissingTargetName;
        const options = object.object.get("options") orelse return error.MissingTargetOptions;

        inline for (handlers) |pair| if (std.mem.eql(u8, name.string, pair.key)) {

            const options_struct = try std.json.parseFromValue(pair.handler.Options, allocator, options, .{}); defer options_struct.deinit();

            try pair.handler.run(f64, options_struct.value, allocator);
        };
    }
}

/// Main function of the program.
pub fn main() !void {
    var timer = try std.time.Timer.start();

    var gpa = std.heap.GeneralPurposeAllocator(.{}){}; const allocator = gpa.allocator();

    defer {
        if (gpa.deinit() == .leak) std.debug.panic("MEMORY LEAK DETECTED IN THE ALLOCATOR\n", .{});
    }

    try print("ZIG VERSION: {d}.{d}.{d}\n", .{builtin.zig_version.major, builtin.zig_version.minor, builtin.zig_version.patch});

    var argc: usize = 0; var argv = try std.process.argsWithAllocator(allocator); defer argv.deinit(); _ = argv.next();

    while (argv.next()) |arg| {
        try parse(arg, allocator); argc += 1;
    }

    if (argc == 0) parse("input.json", allocator) catch |err| {

        if (err != error.FileNotFound) return err;

        try print("\nNO INPUT FILE PROVIDED AND THE DEFAULT \"input.json\" NOT FOUND\n", .{});
    };

    try print("\nTOTAL EXECUTION TIME: {D}\n", .{timer.read()});
}

test {
    std.testing.refAllDecls(@This());
}
