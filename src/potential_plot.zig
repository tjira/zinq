//! Target for performing quantum dynamics simulations.

const std = @import("std");

const device_write = @import("device_write.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const electronic_potential = @import("electronic_potential.zig");
const grid_generator = @import("grid_generator.zig");
const real_matrix = @import("real_matrix.zig");

const ElectronicPotential = electronic_potential.ElectronicPotential;
const RealMatrix = real_matrix.RealMatrix;

const exportRealMatrix = device_write.exportRealMatrix;
const positionAtRow = grid_generator.positionAtRow;
const print = device_write.print;
const printJson = device_write.printJson;

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
pub fn run(comptime T: type, raw_options: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(raw_options);

    var opt = raw_options; try opt.potential.init(allocator); defer opt.potential.deinit(allocator);

    if (opt.potential == .ab_initio) {

        std.log.err("AB INITIO POTENTIAL IS NOT SUPPORTED IN THE POTENTIAL PLOT TARGET", .{});

        return error.InvalidInput;
    }

    const ndim = try opt.potential.ndim();
    const nstate = opt.potential.nstate();

    if (opt.grid.limits.len != ndim) {

        std.log.err("GRID LIMITS MUST HAVE THE SAME LENGTH AS THE NUMBER OF DIMENSIONS, BUT GOT {d} AND {d}", .{opt.grid.limits.len, ndim});

        return error.InvalidInput;
    }


    for (0..ndim) |i| if (opt.grid.limits[i].len != 2) {

        std.log.err("EACH GRID LIMIT MUST HAVE EXACTLY 2 ELEMENTS, BUT GOT {d} FOR DIMENSION {d}", .{opt.grid.limits[i].len, i});

        return error.InvalidInput;
    };

    const potential_matrix = try RealMatrix(T).init(std.math.pow(usize, @intCast(opt.grid.points), ndim), ndim + nstate * nstate, allocator);

    for (0..potential_matrix.rows) |i| {

        var point = potential_matrix.row(i).slice(0, ndim);
        var value = potential_matrix.row(i).slice(ndim, ndim + nstate * nstate).asMatrix();

        try value.reshape(nstate, nstate);

        positionAtRow(T, &point, i, ndim, @intCast(opt.grid.points), opt.grid.limits);

        if (opt.adiabatic) try opt.potential.evaluateAdiabatic(&value, point, 0) else try opt.potential.evaluateDiabatic(&value, point, 0);
    }

    try exportRealMatrix(T, opt.output, potential_matrix);

    return .{.potential = potential_matrix};
}
