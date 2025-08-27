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

/// Parse the input JSON file and run the corresponding target.
pub fn parse(file_contents: []const u8, allocator: std.mem.Allocator) !void {
    const input_json = try std.json.parseFromSlice(std.json.Value, allocator, file_contents, .{}); defer input_json.deinit();

    for (input_json.value.object.get("zinq").?.array.items) |object| {

        const name = object.object.get("name").?.string;
        const options = object.object.get("options").?;

        if (std.mem.eql(u8, name, "classical_dynamics")) {
            const options_struct = try std.json.parseFromValue(classical_dynamics.Options(f64), allocator, options, .{}); defer options_struct.deinit();

            try classical_dynamics.run(f64, options_struct.value, allocator);
        }
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

    var argc: usize = 0; var argv = try std.process.argsWithAllocator(allocator); defer argv.deinit();

    _ = argv.next(); while (argv.next()) |arg| {

        const file_contents = try std.fs.cwd().readFileAlloc(allocator, arg, MAX_INPUT_FILE_BYTES); defer allocator.free(file_contents);

        try print("\nPROCESSED FILE: {s}\n", .{arg});

        try parse(file_contents, allocator); argc += 1;
    }

    default: {if (argc == 0) {

        const filebuf = std.fs.cwd().readFileAlloc(allocator, "input.json", MAX_INPUT_FILE_BYTES) catch break :default;

        try print("\nPROCESSED FILE: {s}\n", .{"input.json"});

        try parse(filebuf, allocator);

        allocator.free(filebuf); 
    }}

    try print("\nTOTAL EXECUTION TIME: {D}\n", .{timer.read()});
}

test {
    std.testing.refAllDecls(@This());
}
