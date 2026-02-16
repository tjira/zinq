//! Main file of the program.

const builtin = @import("builtin");
const config = @import("config");
const std = @import("std");

pub const andersen_thermostat = @import("andersen_thermostat.zig");
pub const array_functions = @import("array_functions.zig");
pub const avoided_crossing_potential = @import("avoided_crossing_potential.zig");
pub const basis_set = @import("basis_set.zig");
pub const berendsen_thermostat = @import("berendsen_thermostat.zig");
pub const bias_potential = @import("bias_potential.zig");
pub const classical_dynamics = @import("classical_dynamics.zig");
pub const classical_particle = @import("classical_particle.zig");
pub const complex_absorbing_potential = @import("complex_absorbing_potential.zig");
pub const complex_gaussian = @import("complex_gaussian.zig");
pub const complex_matrix = @import("complex_matrix.zig");
pub const complex_runge_kutta = @import("complex_runge_kutta.zig");
pub const complex_vector = @import("complex_vector.zig");
pub const contracted_gaussian = @import("contracted_gaussian.zig");
pub const custom_potential = @import("custom_potential.zig");
pub const derivative_coupling = @import("derivative_coupling.zig");
pub const determinant = @import("determinant.zig");
pub const device_read = @import("device_read.zig");
pub const device_write = @import("device_write.zig");
pub const eigenproblem_solver = @import("eigenproblem_solver.zig");
pub const electronic_potential = @import("electronic_potential.zig");
pub const embedded_files = @import("embedded_files.zig");
pub const energy_derivative = @import("energy_derivative.zig");
pub const error_context = @import("error_context.zig");
pub const error_handling = @import("error_handling.zig");
pub const expression_evaluator = @import("expression_evaluator.zig");
pub const fewest_switches = @import("fewest_switches.zig");
pub const file_potential = @import("file_potential.zig");
pub const fourier_transform = @import("fourier_transform.zig");
pub const fractal_generator = @import("fractal_generator.zig");
pub const frequency_analysis = @import("frequency_analysis.zig");
pub const global_variables = @import("global_variables.zig");
pub const grid_generator = @import("grid_generator.zig");
pub const grid_wavefunction = @import("grid_wavefunction.zig");
pub const hammes_schiffer_tully = @import("hammes_schiffer_tully.zig");
pub const harmonic_potential = @import("harmonic_potential.zig");
pub const hartree_fock = @import("hartree_fock.zig");
pub const hermite_quadrature_nodes = @import("hermite_quadrature_nodes.zig");
pub const image = @import("image.zig");
pub const integer_arithmetic = @import("integer_arithmetic.zig");
pub const integral_functions = @import("integral_functions.zig");
pub const integral_transform = @import("integral_transform.zig");
pub const jahn_teller_potential = @import("jahn_teller_potential.zig");
pub const landau_zener = @import("landau_zener.zig");
pub const langevin_thermostat = @import("langevin_thermostat.zig");
pub const linear_interpolation = @import("linear_interpolation.zig");
pub const linear_solve = @import("linear_solve.zig");
pub const mapping_approach = @import("mapping_approach.zig");
pub const math_functions = @import("math_functions.zig");
pub const matrix_inverse = @import("matrix_inverse.zig");
pub const matrix_multiplication = @import("matrix_multiplication.zig");
pub const molecular_integrals = @import("molecular_integrals.zig");
pub const moller_plesset = @import("moller_plesset.zig");
pub const morse_potential = @import("morse_potential.zig");
pub const multiconfigurational_gaussian = @import("multiconfigurational_gaussian.zig");
pub const nonadiabatic_coupling_vector = @import("nonadiabatic_coupling_vector.zig");
pub const norm_preserving_interpolation = @import("norm_preserving_interpolation.zig");
pub const nose_hoover_thermostat = @import("nose_hoover_thermostat.zig");
pub const object_array = @import("object_array.zig");
pub const orbit_fractal = @import("orbit_fractal.zig");
pub const particle_optimization = @import("particle_optimization.zig");
pub const potential_plot = @import("potential_plot.zig");
pub const prime_numbers = @import("prime_numbers.zig");
pub const primitive_gaussian = @import("primitive_gaussian.zig");
pub const process = @import("process.zig");
pub const quantum_dynamics = @import("quantum_dynamics.zig");
pub const real_matrix = @import("real_matrix.zig");
pub const real_tensor_four = @import("real_tensor_four.zig");
pub const real_tensor_three = @import("real_tensor_three.zig");
pub const real_vector = @import("real_vector.zig");
pub const reverse_polish_notation = @import("reverse_polish_notation.zig");
pub const rgb = @import("rgb.zig");
pub const ring_buffer = @import("ring_buffer.zig");
pub const shunting_yard = @import("shunting_yard.zig");
pub const single_set_of_mcg = @import("single_set_of_mcg.zig");
pub const strided_complex_vector = @import("strided_complex_vector.zig");
pub const strided_real_vector = @import("strided_real_vector.zig");
pub const string_manipulation = @import("string_manipulation.zig");
pub const surface_hopping_algorithm = @import("surface_hopping_algorithm.zig");
pub const thermostat = @import("thermostat.zig");
pub const time_linear_potential = @import("time_linear_potential.zig");
pub const timestamp = @import("timestamp.zig");
pub const trajectory_thermodynamics = @import("trajectory_thermodynamics.zig");
pub const tully_potential = @import("tully_potential.zig");
pub const vibronic_coupling_potential = @import("vibronic_coupling_potential.zig");

pub const ClassicalDynamicsOptions = classical_dynamics.Options;
pub const EigenproblemSolverOptions = eigenproblem_solver.Options;
pub const ExpressionEvaluatorOptions = expression_evaluator.Options;
pub const FractalGeneratorOptions = fractal_generator.Options;
pub const HartreeFockOptions = hartree_fock.Options;
pub const MolecularIntegralsOptions = molecular_integrals.Options;
pub const MollerPlessetOptions = moller_plesset.Options;
pub const MulticonfigurationalGaussianOptions = multiconfigurational_gaussian.Options;
pub const PotentialPlotOptions = potential_plot.Options;
pub const PrimeNumbersOptions = prime_numbers.Options;
pub const QuantumDynamicsOptions = quantum_dynamics.Options;

/// Available targets in the program.
const Target = enum {
    classical_dynamics,
    eigenproblem_solver,
    expression_evaluator,
    fractal_generator,
    hartree_fock,
    molecular_integrals,
    moller_plesset,
    multiconfigurational_gaussian,
    potential_plot,
    prime_numbers,
    quantum_dynamics
};

/// Find out if the input JSON file contains an unrecognized field and print the expected field.
fn checkForUnrecognizedFields(comptime Struct: type, options: std.json.Value, err: anyerror) !void {
    inline for (std.meta.fields(Struct)) |field| if (@typeInfo(field.type) == .@"struct") {
        try checkForUnrecognizedFields(field.type, options.object.get(field.name).?, err);
    };

    for (options.object.keys()) |provided| {

        var found = false;

        for (std.meta.fieldNames(Struct)) |expected| if (std.mem.eql(u8, provided, expected)) {
            found = true; break;
        };

        if (!found) return error_handling.throwSpecific(void, "UNRECOGNIZED FIELD '{s}' IN INPUT", .{provided}, err);
    }
}

/// Handle a specific module by parsing its options and running it.
fn handle(comptime T: type, comptime Module: type, options: std.json.Value, allocator: std.mem.Allocator) !void {
    var parsed = std.json.parseFromValue(Module.Options(T), allocator, options, .{}) catch |err| {

        try checkForUnrecognizedFields(Module.Options(T), options, err);

        return error_handling.throwSpecific(void, "UNHANDLED ERROR WHILE PARSING INPUT", .{}, err);
    };

    var output = try Module.run(T, parsed.value, true, allocator); defer output.deinit(allocator);

    parsed.deinit();
}

/// Parse the input JSON file and run the corresponding target.
pub fn parse(path: []const u8, allocator: std.mem.Allocator) !void {
    const file_contents = std.fs.cwd().readFileAlloc(allocator, path, global_variables.MAX_INPUT_FILE_BYTES) catch return error.InputFileNotFound;

    try device_write.print("\nPROCESSED FILE: {s}\n", .{path});

    var scanner = std.json.Scanner.initCompleteInput(allocator, file_contents); defer scanner.deinit();
    var diagnostics = std.json.Diagnostics{}; scanner.enableDiagnostics(&diagnostics);

    const input_json = std.json.parseFromTokenSource(std.json.Value, allocator, &scanner, .{}) catch |err| {
        return error_handling.throwSpecific(void, "ERROR IN INPUT ON LINE {d}", .{diagnostics.line_number}, err);
    };

    for (input_json.value.object.get("zinq").?.array.items, 1..) |object, i| {

        try device_write.print("\nINPUT #{d}:\n", .{i});

        const name = object.object.get("name") orelse return error.MissingTargetName;
        const options = object.object.get("options") orelse return error.MissingTargetOptions;

        const tag = std.meta.stringToEnum(Target, name.string) orelse return error.UnknownTarget;

        switch (tag) {
            .classical_dynamics => try handle(f64, classical_dynamics, options, allocator),
            .eigenproblem_solver => try handle(f64, eigenproblem_solver, options, allocator),
            .expression_evaluator => try handle(f64, expression_evaluator, options, allocator),
            .fractal_generator => try handle(f64, fractal_generator, options, allocator),
            .hartree_fock => try handle(f64, hartree_fock, options, allocator),
            .molecular_integrals => try handle(f64, molecular_integrals, options, allocator),
            .moller_plesset => try handle(f64, moller_plesset, options, allocator),
            .multiconfigurational_gaussian => try handle(f64, multiconfigurational_gaussian, options, allocator),
            .potential_plot => try handle(f64, potential_plot, options, allocator),
            .prime_numbers => try handle(f64, prime_numbers, options, allocator),
            .quantum_dynamics => try handle(f64, quantum_dynamics, options, allocator)
        }
    }

    if (input_json.value.object.get("command")) |command| {

        try device_write.print("\nRUNNING COMMAND: {s}\n", .{command.string});

        const output = try process.executeCommand(command.string, allocator); defer allocator.free(output);

        if (output.len > 0) try device_write.print("COMMAND OUTPUTS:\n{s}", .{output});
    }

    allocator.free(file_contents); input_json.deinit();
}

/// Main function of the program.
pub fn main() !void {
    var timer = try std.time.Timer.start();

    var gpa = std.heap.GeneralPurposeAllocator(.{}){}; const allocator = gpa.allocator();

    defer {
        if (gpa.deinit() == .leak) std.debug.panic("MEMORY LEAK DETECTED IN THE ALLOCATOR\n", .{});
    }

    try device_write.print("ZIG VERSION: {d}.{d}.{d}, ZINQ VERSION: {s}", .{builtin.zig_version.major, builtin.zig_version.minor, builtin.zig_version.patch, config.version});

    {
        const ts = try timestamp.timestamp(allocator); defer allocator.free(ts);

        try device_write.print(", TIMESTAMP: {s}\n", .{ts});
    }

    var argc: usize = 0; var argv = try std.process.argsWithAllocator(allocator); defer argv.deinit(); _ = argv.next();

    while (argv.next()) |arg| {
        try parse(arg, allocator); argc += 1;
    }

    if (argc == 0) parse("input.json", allocator) catch |err| {

        if (err != error.InputFileNotFound) return err;

        try device_write.print("\nNO INPUT FILE PROVIDED AND THE DEFAULT \"input.json\" NOT FOUND\n", .{});
    };

    try device_write.print("\nTOTAL EXECUTION TIME: {D}\n", .{timer.read()});
}

test {
    std.testing.refAllDecls(@This());
}
