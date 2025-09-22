//! Target for performing quantum dynamics simulations.

const std = @import("std");

const device_write = @import("device_write.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const electronic_potential = @import("electronic_potential.zig");
const error_handling = @import("error_handling.zig");
const grid_generator = @import("grid_generator.zig");
const real_matrix = @import("real_matrix.zig");

const ElectronicPotential = electronic_potential.ElectronicPotential;
const RealMatrix = real_matrix.RealMatrix;

const exportRealMatrix = device_write.exportRealMatrix;
const positionAtRow = grid_generator.positionAtRow;
const print = device_write.print;
const printJson = device_write.printJson;
const throw = error_handling.throw;

/// The potential plot options struct.
pub fn Options(comptime T: type) type {
    return struct {
        pub const Grid = struct {
            limits: []const []const T,
            points: u32
        };
        pub const Write = struct {
            population: ?[]const u8 = null,
            wavefunction: ?[]const u8 = null
        };

        potential: ElectronicPotential(T),
        grid: Grid,
        output: []const u8,

        adiabatic: bool = false,
    };
}

/// The quantum dynamics output struct.
pub fn Output(comptime T: type) type {
    return struct {
        potential: RealMatrix(T),

        /// Deallocate the output structure.
        pub fn deinit(self: @This()) void {
            self.potential.deinit();
        }
    };
}

/// Run the potential plot target with the given options.
pub fn run(comptime T: type, options: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(options);

    var potential = options.potential;
    const ndim = potential.ndim();
    const nstate = potential.nstate();

    const file_potential = if (potential == .file) try potential.file.init(allocator) else null; defer if (file_potential) |U| U.deinit();

    if (options.grid.limits.len != ndim) return throw(Output(T), "LIMITS LENGTH MUST BE EQUAL TO NUMBER OF DIMENSIONS", .{});
    for (0..ndim) |i| if (options.grid.limits[i].len != 2) return throw(Output(T), "EACH LIMIT MUST HAVE A LENGTH OF 2", .{});

    const potential_matrix = try RealMatrix(T).init(std.math.pow(usize, @intCast(options.grid.points), ndim), ndim + nstate * nstate, allocator);

    for (0..potential_matrix.rows) |i| {

        var point = potential_matrix.row(i).slice(0, ndim);
        var value = potential_matrix.row(i).slice(ndim, ndim + nstate * nstate).asMatrix();

        try value.reshape(nstate, nstate);

        positionAtRow(T, &point, i, ndim, @intCast(options.grid.points), options.grid.limits);

        if (options.adiabatic) try potential.evaluateAdiabatic(&value, point, 0) else potential.evaluateDiabatic(&value, point, 0);
    }

    try exportRealMatrix(T, options.output, potential_matrix);

    return .{.potential = potential_matrix};
}
