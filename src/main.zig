//! Main file of the program.

const builtin = @import("builtin");
const std = @import("std");

pub const array_functions = @import("array_functions.zig");
pub const basis_set = @import("basis_set.zig");
pub const classical_dynamics = @import("classical_dynamics.zig");
pub const classical_particle = @import("classical_particle.zig");
pub const complex_matrix = @import("complex_matrix.zig");
pub const complex_vector = @import("complex_vector.zig");
pub const contracted_gaussian = @import("contracted_gaussian.zig");
pub const derivative_coupling = @import("derivative_coupling.zig");
pub const device_write = @import("device_write.zig");
pub const eigenproblem_solver = @import("eigenproblem_solver.zig");
pub const electronic_potential = @import("electronic_potential.zig");
pub const fourier_transform = @import("fourier_transform.zig");
pub const global_variables = @import("global_variables.zig");
pub const grid_generator = @import("grid_generator.zig");
pub const hammes_schiffer_tully = @import("hammes_schiffer_tully.zig");
pub const harmonic_potential = @import("harmonic_potential.zig");
pub const integral_functions = @import("integral_functions.zig");
pub const landau_zener = @import("landau_zener.zig");
pub const math_functions = @import("math_functions.zig");
pub const matrix_multiplication = @import("matrix_multiplication.zig");
pub const molecular_integrals = @import("molecular_integrals.zig");
pub const norm_preserving_interpolation = @import("norm_preserving_interpolation.zig");
pub const object_array = @import("object_array.zig");
pub const prime_numbers = @import("prime_numbers.zig");
pub const primitive_gaussian = @import("primitive_gaussian.zig");
pub const quantum_dynamics = @import("quantum_dynamics.zig");
pub const real_matrix = @import("real_matrix.zig");
pub const real_tensor_four = @import("real_tensor_four.zig");
pub const real_tensor_three = @import("real_tensor_three.zig");
pub const real_vector = @import("real_vector.zig");
pub const ring_buffer = @import("ring_buffer.zig");
pub const strided_complex_vector = @import("strided_complex_vector.zig");
pub const surface_hopping_algorithm = @import("surface_hopping_algorithm.zig");

/// Available targets in the program.
const Target = enum {
    classical_dynamics,
    molecular_integrals,
    prime_numbers,
    quantum_dynamics
};

/// Handle a specific module by parsing its options and running it.
fn handle(comptime T: type, comptime Module: type, options: std.json.Value, allocator: std.mem.Allocator) !void {
    var parsed = try std.json.parseFromValue(Module.Options(T), allocator, options, .{}); defer parsed.deinit();

    var output = try Module.run(T, parsed.value, true, allocator); defer output.deinit();
}

/// Parse the input JSON file and run the corresponding target.
pub fn parse(path: []const u8, allocator: std.mem.Allocator) !void {
    const file_contents = try std.fs.cwd().readFileAlloc(allocator, path, global_variables.MAX_INPUT_FILE_BYTES); defer allocator.free(file_contents);

    try device_write.print("\nPROCESSED FILE: {s}\n", .{path});

    const input_json = try std.json.parseFromSlice(std.json.Value, allocator, file_contents, .{}); defer input_json.deinit();

    for (input_json.value.object.get("zinq").?.array.items) |object| {

        const name = object.object.get("name") orelse return error.MissingTargetName;
        const options = object.object.get("options") orelse return error.MissingTargetOptions;

        const tag = std.meta.stringToEnum(Target, name.string) orelse return error.UnknownTarget;

        switch (tag) {
            .classical_dynamics => try handle(f64, classical_dynamics, options, allocator),
            .molecular_integrals => try handle(f64, molecular_integrals, options, allocator),
            .prime_numbers => try handle(u64, prime_numbers, options, allocator),
            .quantum_dynamics => try handle(f64, quantum_dynamics, options, allocator)
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

    try device_write.print("ZIG VERSION: {d}.{d}.{d}\n", .{builtin.zig_version.major, builtin.zig_version.minor, builtin.zig_version.patch});

    var argc: usize = 0; var argv = try std.process.argsWithAllocator(allocator); defer argv.deinit(); _ = argv.next();

    while (argv.next()) |arg| {
        try parse(arg, allocator); argc += 1;
    }

    if (argc == 0) parse("input.json", allocator) catch |err| {

        if (err != error.FileNotFound) return err;

        try device_write.print("\nNO INPUT FILE PROVIDED AND THE DEFAULT \"input.json\" NOT FOUND\n", .{});
    };

    try device_write.print("\nTOTAL EXECUTION TIME: {D}\n", .{timer.read()});
}

test {
    std.testing.refAllDecls(@This());
}
