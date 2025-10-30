//! Main file for the vMCG program.

const std = @import("std");

const complex_gaussian = @import("complex_gaussian.zig");
const complex_matrix = @import("complex_matrix.zig");
const complex_vector = @import("complex_vector.zig");
const device_write = @import("device_write.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const electronic_potential = @import("electronic_potential.zig");
const global_variables = @import("global_variables.zig");
const grid_generator = @import("grid_generator.zig");
const grid_wavefunction = @import("grid_wavefunction.zig");
const harmonic_potential = @import("harmonic_potential.zig");
const math_functions = @import("math_functions.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const object_array = @import("object_array.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");
const tully_potential = @import("tully_potential.zig");

const Complex = std.math.Complex;
const ComplexGaussian = complex_gaussian.ComplexGaussian;
const ComplexMatrix = complex_matrix.ComplexMatrix;
const ComplexVector = complex_vector.ComplexVector;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const GridWavefunction = grid_wavefunction.GridWavefunction;
const HarmonicPotential = harmonic_potential.HarmonicPotential;
const RealMatrix = real_matrix.RealMatrix;
const RealMatrixArray = object_array.RealMatrixArray;
const RealVector = real_vector.RealVector;
const TullyPotential1 = tully_potential.TullyPotential1;

const exportRealMatrix = device_write.exportRealMatrix;
const exportRealMatrixWithLinspacedLeftColumn = device_write.exportRealMatrixWithLinspacedLeftColumn;
const fixGauge = eigenproblem_solver.fixGauge;
const mm = matrix_multiplication.mm;
const momentumGridAlloc = grid_generator.momentumGridAlloc;
const positionAtRow = grid_generator.positionAtRow;
const positionGridAlloc = grid_generator.positionGridAlloc;
const powi = math_functions.powi;
const print = device_write.print;
const printJson = device_write.printJson;
const printRealMatrix = device_write.printRealMatrix;
const printComplexMatrix = device_write.printComplexMatrix;

/// The vMCG dynamics opt struct.
pub fn Options(comptime T: type) type {
    return struct {
        pub const IntegrationGrid = struct {
            limits: []const []const T,
            points: u32
        };
        pub const LogIntervals = struct {
            iteration: u32 = 1
        };
        pub const Write = struct {
            population: ?[]const u8 = null,
            wavefunction: ?[]const u8 = null
        };
        pub const Method = union(enum) {
            vMCG: struct {
                position: []const T,
                gamma: []const T,
                momentum: []const T
            }
        };

        potential: ElectronicPotential(T),
        integration_grid: IntegrationGrid,
        method: Method,

        initial_state: u32,
        iterations: u32,
        time_step: T,
        mass: []const T,

        log_intervals: LogIntervals = .{},
        write: Write = .{},

        adiabatic: bool = false,
    };
}

/// The vMCG dynamics output struct.
pub fn Output(comptime T: type) type {
    return struct {
        kinetic_energy: T = undefined,
        population: RealMatrix(T),
        potential_energy: T = undefined,

        /// Allocate the output structure.
        pub fn init(nstate: usize, iterations: usize, allocator: std.mem.Allocator) !@This() {
            return @This(){
                .population = try RealMatrix(T).init(iterations + 1, nstate, allocator),
            };
        }

        /// Deallocate the output structure.
        pub fn deinit(self: @This()) void {
            self.population.deinit();
        }
    };
}

/// Container for custom structs related to the vMCG.
pub fn Custom(comptime T: type) type {
    return struct {

        /// Structure to hold information about each iteration used for logging.
        pub const IterationInfo = struct {
            iteration: usize,
            kinetic_energy: T,
            coefs: ComplexVector(T),
            momentum: RealVector(T),
            position: RealVector(T),
            potential_energy: T,
            time: T
        };
    };
}

/// Run quantum dynamics simulation.
pub fn run(comptime T: type, opt: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(opt);

    switch (opt.method) {
        .vMCG => return try runSingleGaussian(T, opt, enable_printing, allocator)
    }
}

/// Perform the simulation with a single gaussian.
pub fn runSingleGaussian(comptime T: type, opt: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(opt);

    var potential = opt.potential;
    const ndim = potential.ndim();
    const nstate = potential.nstate();

    const file_potential = if (potential == .file) try potential.file.init(allocator) else null; defer if (file_potential) |U| U.deinit();

    var output = try Output(T).init(nstate, @intCast(opt.iterations), allocator); var opt_copy = opt; opt_copy.potential = potential;

    var wavefunction_dynamics: ?RealMatrix(T) = if (opt.write.wavefunction) |_| try initializeWavefunctionDynamicsContainer(T, opt.integration_grid, nstate, opt.iterations, allocator) else null;

    var gaussian = try complex_gaussian.ComplexGaussian(T).init(opt.method.vMCG.position, opt.method.vMCG.gamma, opt.method.vMCG.momentum, allocator); defer gaussian.deinit();

    var coefs = try ComplexVector(T).init(nstate, allocator); defer coefs.deinit(); coefs.ptr(opt.initial_state).* = Complex(T).init(1, 0);

    var position = try RealVector(T).init(ndim, allocator); defer position.deinit();
    var momentum = try RealVector(T).init(ndim, allocator); defer momentum.deinit();

    try printIterationHeader(ndim, nstate);

    var timer = try std.time.Timer.start();

    for (0..opt.iterations + 1) |i| {

        const time = @as(T, @floatFromInt(i)) * opt.time_step;

        if (i > 0) try propagateSingleGaussian(T, &gaussian, &coefs, opt_copy, allocator);

        const V = try gaussian.potential(gaussian, potential, opt.integration_grid.limits, opt.integration_grid.points, 0, allocator); defer V.deinit();

        const potential_energy = try potentialEnergySingleGaussian(T, V, coefs);
        const kinetic_energy = (try gaussian.kinetic(gaussian, opt.mass)).re;

        for (0..ndim) |j| {
            position.ptr(j).* = gaussian.center[j]; momentum.ptr(j).* = gaussian.momentum[j];
        }

        for (0..nstate) |j| {
            output.population.ptr(i, j).* = coefs.at(j).magnitude() * coefs.at(j).magnitude();
        }

        if (opt.write.wavefunction != null) try assignWavefunctionStepSingleGaussian(T, &wavefunction_dynamics.?, gaussian, coefs, opt.integration_grid, i, allocator);

        if (i == opt.iterations) {
            output.kinetic_energy = kinetic_energy;
            output.potential_energy = potential_energy;
        }

        if (!enable_printing or (i > 0 and i % opt.log_intervals.iteration != 0)) continue;

        const info = Custom(T).IterationInfo{
            .iteration = i,
            .kinetic_energy = kinetic_energy,
            .coefs = coefs,
            .momentum = momentum,
            .position = position,
            .potential_energy = potential_energy,
            .time = time
        };

        try printIterationInfo(T, info, &timer);
    }

    if (enable_printing) for (0..nstate) |i| {
        try print("{s}FINAL POPULATION OF STATE {d}: {d:.6}\n", .{if (i == 0) "\n" else "", i, output.population.at(opt.iterations, i)});
    };

    const end_time = @as(T, @floatFromInt(opt.iterations)) * opt.time_step;

    if (opt.write.population) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.population, 0, end_time, opt.iterations + 1);
    if (opt.write.wavefunction) |path| try exportRealMatrix(T, path, wavefunction_dynamics.?, );

    if (wavefunction_dynamics != null) wavefunction_dynamics.?.deinit();

    return output;
}

/// Assigns the single gaussian to the wavefunction dynamics container.
pub fn assignWavefunctionStepSingleGaussian(comptime T: type, wfn: *RealMatrix(T), gaussian: ComplexGaussian(T), coefs: ComplexVector(T), grid: Options(T).IntegrationGrid, iteration: usize, allocator: std.mem.Allocator) !void {
    const ndim = gaussian.center.len;
    const nstate = coefs.len;

    var position_at_row = try RealVector(T).init(ndim, allocator); defer position_at_row.deinit();

    for (0..wfn.rows) |i| {

        positionAtRow(T, &position_at_row, i, ndim, grid.points, grid.limits);

        for (0..nstate) |j| {

            const value = gaussian.evaluate(position_at_row.data).mul(coefs.at(j));

            wfn.ptr(i, ndim + 2 * nstate * iteration + 2 * j + 0).* = value.re;
            wfn.ptr(i, ndim + 2 * nstate * iteration + 2 * j + 1).* = value.im;
        }
    }
}

/// Initialize the container for the wavefunction dynamics.
pub fn initializeWavefunctionDynamicsContainer(comptime T: type, grid: Options(T).IntegrationGrid, nstate: usize, iterations: usize, allocator: std.mem.Allocator) !RealMatrix(T) {
    const grid_points = powi(grid.points, grid.limits.len);

    var wavefunction_dynamics: RealMatrix(T) = try RealMatrix(T).init(grid_points, grid.limits.len + 2 * nstate * (iterations + 1), allocator);

    var position_at_row = try RealVector(T).init(grid.limits.len, allocator); defer position_at_row.deinit();

    for (0..grid_points) |i| {

        positionAtRow(T, &position_at_row, i, grid.limits.len, grid.points, grid.limits);

        for (0..grid.limits.len) |j| {
            wavefunction_dynamics.ptr(i, j).* = position_at_row.at(j);
        }
    }

    return wavefunction_dynamics;
}

/// Propagates the gaussian and coefficients using the 4th order Runge-Kutta method.
pub fn propagateSingleGaussian(comptime T: type, gaussian: *ComplexGaussian(T), coefs: *ComplexVector(T), opt: Options(T), allocator: std.mem.Allocator) !void {
    const k1p = try parameterDerivativesSingleGaussian(T, gaussian.*, opt, coefs.*, allocator); defer k1p.dq.deinit(); defer k1p.dp.deinit();
    const k1c = try coefficientDerivativesSingleGaussian(T, gaussian.*, opt, coefs.*, allocator); defer k1c.deinit();
    const g1 = try gaussian.clone(); defer g1.deinit();
    for (0..opt.potential.ndim()) |i| {
        g1.center[i] += 0.5 * opt.time_step * k1p.dq.at(i);
        g1.momentum[i] += 0.5 * opt.time_step * k1p.dp.at(i);
    }
    const c1 = try coefs.clone(); defer c1.deinit();
    for (0..coefs.len) |i| {
        c1.ptr(i).* = c1.at(i).add(k1c.at(i).mul(Complex(T).init(0.5 * opt.time_step, 0)));
    }

    const k2p = try parameterDerivativesSingleGaussian(T, g1, opt, c1, allocator); defer k2p.dq.deinit(); defer k2p.dp.deinit();
    const k2c = try coefficientDerivativesSingleGaussian(T, g1, opt, c1, allocator); defer k2c.deinit();
    const g2 = try gaussian.clone(); defer g2.deinit();
    for (0..opt.potential.ndim()) |i| {
        g2.center[i] += 0.5 * opt.time_step * k2p.dq.at(i);
        g2.momentum[i] += 0.5 * opt.time_step * k2p.dp.at(i);
    }
    const c2 = try coefs.clone(); defer c2.deinit();
    for (0..coefs.len) |i| {
        c2.ptr(i).* = c2.at(i).add(k2c.at(i).mul(Complex(T).init(0.5 * opt.time_step, 0)));
    }

    const k3p = try parameterDerivativesSingleGaussian(T, g2, opt, c2, allocator); defer k3p.dq.deinit(); defer k3p.dp.deinit();
    const k3c = try coefficientDerivativesSingleGaussian(T, g2, opt, c2, allocator); defer k3c.deinit();
    const g3 = try gaussian.clone(); defer g3.deinit();
    for (0..opt.potential.ndim()) |i| {
        g3.center[i] += opt.time_step * k3p.dq.at(i);
        g3.momentum[i] += opt.time_step * k3p.dp.at(i);
    }
    const c3 = try coefs.clone(); defer c3.deinit();
    for (0..coefs.len) |i| {
        c3.ptr(i).* = c3.at(i).add(k3c.at(i).mul(Complex(T).init(opt.time_step, 0)));
    }

    const k4p = try parameterDerivativesSingleGaussian(T, g3, opt, c3, allocator); defer k4p.dq.deinit(); defer k4p.dp.deinit();
    const k4c = try coefficientDerivativesSingleGaussian(T, g3, opt, c3, allocator); defer k4c.deinit();

    for (0..opt.potential.ndim()) |i| {
        gaussian.center[i] += opt.time_step * (k1p.dq.at(i) + 2 * k2p.dq.at(i) + 2 * k3p.dq.at(i) + k4p.dq.at(i)) / 6;
        gaussian.momentum[i] += opt.time_step * (k1p.dp.at(i) + 2 * k2p.dp.at(i) + 2 * k3p.dp.at(i) + k4p.dp.at(i)) / 6;
    }

    for (0..coefs.len) |i| {
        coefs.ptr(i).* = coefs.at(i).add(k1c.at(i).add(k2c.at(i).mul(Complex(T).init(2, 0))).add(k3c.at(i).mul(Complex(T).init(2, 0))).add(k4c.at(i)).mul(Complex(T).init(opt.time_step / 6, 0)));
    }
}

/// Returns the derivative of coefficients for a single gaussian given the potential energy matrix and coefficients.
pub fn coefficientDerivativesSingleGaussian(comptime T: type, g: ComplexGaussian(T), opt: Options(T), coefs: ComplexVector(T), allocator: std.mem.Allocator) !ComplexVector(T) {
    var dc = try ComplexVector(T).initZero(coefs.len, allocator);

    const V = try g.potential(g, opt.potential, opt.integration_grid.limits, opt.integration_grid.points, 0, allocator); defer V.deinit();

    var rho = try ComplexMatrix(T).initZero(V.rows, V.cols, dc.allocator); defer rho.deinit();
    var VC = try ComplexVector(T).init(V.cols, dc.allocator); defer VC.deinit();

    var VC_matrix = VC.asMatrix();
    const coefs_matrix = coefs.asMatrix();

    for (0..coefs.len) |i| for (0..coefs.len) |j| {
        rho.ptr(i, j).* = coefs.at(i).mul(coefs.at(j).conjugate());
    };

    for (0..rho.cols) |i| for (0..rho.rows) |j| {
        rho.ptr(i, j).* = if (i == j) Complex(T).init(1, 0).sub(rho.at(i, i)) else rho.at(i, j).neg();
    };

    try mm(T, &VC_matrix, V, false, coefs_matrix, false);

    for (0..dc.len) |i| for (0..coefs.len) |j| {
        dc.ptr(i).* = dc.at(i).add(rho.at(i, j).mul(VC.at(j)).mulbyi().neg());
    };

    return dc;
}

/// Returns the force acting on a gaussian due to the potential energy matrix and coefficients.
pub fn forceSingleGaussian(comptime T: type, dV: []const ComplexMatrix(T), coefs: ComplexVector(T), allocator: std.mem.Allocator) !RealVector(T) {
    var force = try RealVector(T).initZero(dV.len, allocator);

    for (0..force.len) |i| for (0..coefs.len) |j| for (0..coefs.len) |k| {
        force.ptr(i).* -= coefs.at(j).conjugate().mul(dV[i].at(j, k)).mul(coefs.at(k)).re;
    };

    return force;
}

/// Returns the potential energy of a single complex Gaussian given its matrix element of the potential operator and coefficients.
pub fn potentialEnergySingleGaussian(comptime T: type, V: ComplexMatrix(T), coefs: ComplexVector(T)) !T {
    var potential_energy: Complex(T) = Complex(T).init(0, 0);

    for (0..coefs.len) |i| for (0..coefs.len) |j| {
        potential_energy = potential_energy.add(coefs.at(i).conjugate().mul(V.at(i, j)).mul(coefs.at(j)));
    };

    return potential_energy.re;
}

/// Calculates the derivative of the parameters of a single gaussian given the potential derivative matrix and coefficients.
pub fn parameterDerivativesSingleGaussian(comptime T: type, g: ComplexGaussian(T), opt: Options(T), coefs: ComplexVector(T), allocator: std.mem.Allocator) !struct {dq: RealVector(T), dp: RealVector(T)} {
    var dq = try RealVector(T).initZero(g.center.len, allocator);
    var dp = try RealVector(T).initZero(g.center.len, allocator);

    const dV = try allocator.alloc(ComplexMatrix(T), opt.potential.ndim()); defer allocator.free(dV); defer for (0..opt.potential.ndim()) |j| dV[j].deinit();

    for (0..dV.len) |j| dV[j] = try g.potentialDerivative(g, opt.potential, j, opt.integration_grid.limits, opt.integration_grid.points, 0, allocator);

    const force = try forceSingleGaussian(T, dV, coefs, allocator); defer force.deinit();

    for (0..dV.len) |i| {
        dq.ptr(i).* = g.momentum[i] / opt.mass[i];
        dp.ptr(i).* = force.at(i);
    }

    return .{.dq = dq, .dp = dp};
}

/// Print header for iteration info.
pub fn printIterationHeader(ndim: usize, nstate: usize) !void {
    var buffer: [128]u8 = undefined;

    var writer = std.io.Writer.fixed(&buffer);

    try writer.print("\n{s:8} ", .{"ITER"});
    try writer.print("{s:12} {s:12} {s:12} ", .{"KINETIC", "POTENTIAL", "TOTAL"});
    try writer.print("{[value]s:[width]} ", .{.value = "POSITION", .width = 9 * ndim + 2 * (ndim - 1) + 2});
    try writer.print("{[value]s:[width]} ", .{.value = "MOMENTUM", .width = 9 * ndim + 2 * (ndim - 1) + 2});
    try writer.print("{[value]s:[width]} ", .{.value = "POPULATION", .width = 9 * nstate + 2 * (nstate - 1) + 2});
    try writer.print("{s:4}", .{"TIME"});

    try print("{s}\n", .{writer.buffered()});
}

/// Prints the iteration info to standard output.
pub fn printIterationInfo(comptime T: type, info: Custom(T).IterationInfo, timer: *std.time.Timer) !void {
    var buffer: [128]u8 = undefined;

    var writer = std.io.Writer.fixed(&buffer);

    try writer.print("{d:8} ", .{info.iteration});
    try writer.print("{d:12.6} {d:12.6} {d:12.6} ", .{info.kinetic_energy, info.potential_energy, info.kinetic_energy + info.potential_energy});

    try writer.print("[", .{});

    for (0..info.position.len) |i| {
        try writer.print("{d:9.4}{s}", .{info.position.at(i), if (i == info.position.len - 1) "" else ", "});
    }

    try writer.print("] [", .{});

    for (0..info.momentum.len) |i| {
        try writer.print("{d:9.4}{s}", .{info.momentum.at(i), if (i == info.momentum.len - 1) "" else ", "});
    }

    try writer.print("] [", .{});

    for (0..info.coefs.len) |i| {
        try writer.print("{d:9.4}{s}", .{info.coefs.at(i).magnitude() * info.coefs.at(i).magnitude(), if (i == info.coefs.len - 1) "" else ", "});
    }

    try writer.print("] {D}", .{timer.read()}); timer.reset();

    try print("{s}\n", .{writer.buffered()});
}
