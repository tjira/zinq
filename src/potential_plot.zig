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
        pub fn deinit(self: @This(), allocator: std.mem.Allocator) void {
            self.potential.deinit(allocator);
        }
    };
}

/// Run the potential plot target with the given opt.
pub fn run(comptime T: type, opt: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(opt);

    const ndim = opt.potential.ndim();
    const nstate = opt.potential.nstate();

    var custom_potential = if (opt.potential == .custom) try opt.potential.custom.init(allocator) else null; defer if (custom_potential) |*cp| cp.deinit(allocator);
    var file_potential = if (opt.potential == .file) try opt.potential.file.init(allocator) else null; defer if (file_potential) |*fp| fp.deinit(allocator);

    if (opt.grid.limits.len != ndim) return throw(Output(T), "LIMITS LENGTH MUST BE EQUAL TO NUMBER OF DIMENSIONS", .{});
    for (0..ndim) |i| if (opt.grid.limits[i].len != 2) return throw(Output(T), "EACH LIMIT MUST HAVE A LENGTH OF 2", .{});

    const potential_matrix = try RealMatrix(T).init(std.math.pow(usize, @intCast(opt.grid.points), ndim), ndim + nstate * nstate, allocator);

    for (0..potential_matrix.rows) |i| {

        var point = potential_matrix.row(i).slice(0, ndim);
        var value = potential_matrix.row(i).slice(ndim, ndim + nstate * nstate).asMatrix();

        try value.reshape(nstate, nstate);

        positionAtRow(T, &point, i, ndim, @intCast(opt.grid.points), opt.grid.limits);

        if (opt.adiabatic) try opt.potential.evaluateAdiabatic(&value, point, 0, null) else try opt.potential.evaluateDiabatic(&value, point, 0);
    }

    try exportRealMatrix(T, opt.output, potential_matrix);

    return .{.potential = potential_matrix};
}
