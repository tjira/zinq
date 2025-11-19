//! Main file for the vMCG program.

const std = @import("std");

const complex_gaussian = @import("complex_gaussian.zig");
const complex_runge_kutta = @import("complex_runge_kutta.zig");
const complex_vector = @import("complex_vector.zig");
const device_write = @import("device_write.zig");
const electronic_potential = @import("electronic_potential.zig");
const error_handling = @import("error_handling.zig");
const hermite_quadrature_nodes = @import("hermite_quadrature_nodes.zig");
const global_variables = @import("global_variables.zig");
const grid_generator = @import("grid_generator.zig");
const harmonic_potential = @import("harmonic_potential.zig");
const math_functions = @import("math_functions.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");
const tully_potential = @import("tully_potential.zig");

const Complex = std.math.Complex;
const ComplexGaussian = complex_gaussian.ComplexGaussian;
const ComplexRungeKutta = complex_runge_kutta.ComplexRungeKutta;
const ComplexVector = complex_vector.ComplexVector;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const HarmonicPotential = harmonic_potential.HarmonicPotential;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;
const TullyPotential1 = tully_potential.TullyPotential1;

const exportRealMatrix = device_write.exportRealMatrix;
const exportRealMatrixWithLinspacedLeftColumn = device_write.exportRealMatrixWithLinspacedLeftColumn;
const mm = matrix_multiplication.mm;
const positionAtRow = grid_generator.positionAtRow;
const powi = math_functions.powi;
const print = device_write.print;
const printJson = device_write.printJson;
const throw = error_handling.throw;

const TEST_TOLERANCE = global_variables.TEST_TOLERANCE;
const WRITE_BUFFER_SIZE = global_variables.WRITE_BUFFER_SIZE;
const HERMITE_NODES = hermite_quadrature_nodes.HERMITE_NODES;

/// The vMCG dynamics opt struct.
pub fn Options(comptime T: type) type {
    return struct {
        pub const InitialConditions = struct {
            adiabatic: bool = false,
            mass: []const T,
            state: u32
        };
        pub const LogIntervals = struct {
            iteration: u32 = 1
        };
        pub const Write = struct {
            population: ?[]const u8 = null,
            wavefunction: ?[]const u8 = null
        };
        pub const WavefunctionGrid = struct {
            limits: []const []const T,
            points: u32
        };
        pub const Method = union(enum) {
            vMCG: struct {
                position: []const T,
                gamma: []const T,
                momentum: []const T,
                frozen: bool = true
            }
        };

        initial_conditions: InitialConditions,
        method: Method,
        potential: ElectronicPotential(T),

        iterations: u32,
        time_step: T,

        log_intervals: LogIntervals = .{},
        wavefunction_grid: ?WavefunctionGrid = null,
        write: Write = .{},

        finite_differences_step: T = 1e-6,
        integration_nodes: u32 = 32,
        adiabatic: bool = false,
        imaginary: bool = false
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

    if (opt.integration_nodes < 1 or opt.integration_nodes > 256) return throw(Output(T), "INTEGRATION NODES MUST BE BETWEEN 1 AND {d}", .{HERMITE_NODES.len});

    const ndim = opt.potential.ndim();
    const nstate = opt.potential.nstate();

    var custom_potential = if (opt.potential == .custom) try opt.potential.custom.init(allocator) else null; defer if (custom_potential) |*cp| cp.deinit();
    var file_potential = if (opt.potential == .file) try opt.potential.file.init(allocator) else null; defer if (file_potential) |*fp| fp.deinit();

    var output = try Output(T).init(nstate, @intCast(opt.iterations), allocator);

    var wavefunction_dynamics: ?RealMatrix(T) = if (opt.write.wavefunction) |_| try initializeWavefunctionDynamicsContainer(T, opt.wavefunction_grid, nstate, opt.iterations, allocator) else null;

    var gaussian = try complex_gaussian.ComplexGaussian(T).init(opt.method.vMCG.position, opt.method.vMCG.gamma, opt.method.vMCG.momentum, allocator); defer gaussian.deinit();

    var coefs = try ComplexVector(T).init(nstate, allocator); defer coefs.deinit(); coefs.ptr(opt.initial_conditions.state).* = Complex(T).init(1, 0);

    var propagator = try ComplexRungeKutta(T).init(3 * gaussian.position.len + coefs.len, allocator); defer propagator.deinit();

    if (enable_printing) try printIterationHeader(ndim, nstate);

    var timer = try std.time.Timer.start();

    for (0..opt.iterations + 1) |i| {

        const time = @as(T, @floatFromInt(i)) * opt.time_step;

        if (i > 0) try propagateSingleGaussian(T, &gaussian, &coefs, opt, &propagator, time, allocator);

        const position = RealVector(T){.data = gaussian.position, .len = gaussian.position.len, .allocator = null};
        const momentum = RealVector(T){.data = gaussian.momentum, .len = gaussian.momentum.len, .allocator = null};

        var coefs_adia: ?ComplexVector(T) = null; defer if (coefs_adia) |c| c.deinit();

        if (i == 0 and opt.initial_conditions.adiabatic) {
            coefs_adia = try transformCoefficients(T, coefs, opt.potential, position, time, false, allocator); try coefs_adia.?.copyTo(&coefs);
        }

        if (opt.adiabatic) {
            if (coefs_adia) |c| c.deinit(); coefs_adia = try transformCoefficients(T, coefs, opt.potential, position, time, true, allocator);
        }

        const potential_energy = try gaussian.potentialEnergy(opt.potential, coefs, opt.integration_nodes, time);
        const kinetic_energy = try gaussian.kineticEnergy(opt.initial_conditions.mass);

        for (0..nstate) |j| {
            output.population.ptr(i, j).* = if (opt.adiabatic) coefs_adia.?.at(j).magnitude() * coefs_adia.?.at(j).magnitude() else coefs.at(j).magnitude() * coefs.at(j).magnitude();
        }

        if (wavefunction_dynamics) |*wfn_container| {
            try assignWavefunctionStepSingleGaussian(T, wfn_container, gaussian, if (opt.adiabatic) coefs_adia.? else coefs, i);
        }

        if (i == opt.iterations) {
            output.kinetic_energy = kinetic_energy;
            output.potential_energy = potential_energy;
        }

        if (!enable_printing or (i > 0 and i % opt.log_intervals.iteration != 0)) continue;

        const info = Custom(T).IterationInfo{
            .iteration = i,
            .kinetic_energy = kinetic_energy,
            .coefs = if (opt.adiabatic) coefs_adia.? else coefs,
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

/// Transforms the coefficients between the diabatic and adiabatic representations.
pub fn transformCoefficients(comptime T: type, coefs: ComplexVector(T), potential: ElectronicPotential(T), position: RealVector(T), time: T, to_adiabatic: bool, allocator: std.mem.Allocator) !ComplexVector(T) {
    const coefs_adia = try ComplexVector(T).initZero(coefs.len, allocator); var coefs_adia_matrix = coefs_adia.asMatrix();

    var diabatic_potential = try RealMatrix(T).init(coefs.len, coefs.len, allocator); defer diabatic_potential.deinit();
    var adiabatic_potential = try RealMatrix(T).init(coefs.len, coefs.len, allocator); defer adiabatic_potential.deinit();
    var adiabatic_eigenvectors = try RealMatrix(T).init(coefs.len, coefs.len, allocator); defer adiabatic_eigenvectors.deinit();

    try potential.evaluateEigensystem(&diabatic_potential, &adiabatic_potential, &adiabatic_eigenvectors, position, time);

    if (to_adiabatic) {
        try mm(T, &coefs_adia_matrix, adiabatic_eigenvectors, true, coefs.asMatrix(), false);
    } else {
        try mm(T, &coefs_adia_matrix, adiabatic_eigenvectors, false, coefs.asMatrix(), false);
    }

    return coefs_adia;
}

/// Assigns the single gaussian to the wavefunction dynamics container.
pub fn assignWavefunctionStepSingleGaussian(comptime T: type, wfn_container: *RealMatrix(T), gaussian: ComplexGaussian(T), coefs: ComplexVector(T), iteration: usize) !void {
    for (0..wfn_container.rows) |i| {

        const gv = gaussian.evaluate(wfn_container.row(i).slice(0, gaussian.position.len).data);

        for (0..coefs.len) |j| {

            const value = gv.mul(coefs.at(j));

            wfn_container.ptr(i, gaussian.position.len + 2 * coefs.len * iteration + 2 * j + 0).* = value.re;
            wfn_container.ptr(i, gaussian.position.len + 2 * coefs.len * iteration + 2 * j + 1).* = value.im;
        }
    }
}

/// Initialize the container for the wavefunction dynamics.
pub fn initializeWavefunctionDynamicsContainer(comptime T: type, grid: ?Options(T).WavefunctionGrid, nstate: usize, iterations: usize, allocator: std.mem.Allocator) !RealMatrix(T) {
    if (grid == null) return throw(RealMatrix(T), "WAVEFUNCTION GRID MUST BE PROVIDED WHEN WRITING THE WAVEFUNCTION DYNAMICS", .{});

    const grid_points = powi(grid.?.points, grid.?.limits.len);

    var wavefunction_dynamics: RealMatrix(T) = try RealMatrix(T).init(grid_points, grid.?.limits.len + 2 * nstate * (iterations + 1), allocator);

    var position_at_row = try RealVector(T).init(grid.?.limits.len, allocator); defer position_at_row.deinit();

    for (0..grid_points) |i| {

        positionAtRow(T, &position_at_row, i, grid.?.limits.len, grid.?.points, grid.?.limits);

        for (0..grid.?.limits.len) |j| {
            wavefunction_dynamics.ptr(i, j).* = position_at_row.at(j);
        }
    }

    return wavefunction_dynamics;
}

/// Propagates the gaussian and coefficients using the 4th order Runge-Kutta method.
pub fn propagateSingleGaussian(comptime T: type, gaussian: *ComplexGaussian(T), coefs: *ComplexVector(T), opt: Options(T), rk: *ComplexRungeKutta(T), time: T, allocator: std.mem.Allocator) !void {
    var vars = try ComplexVector(T).init(3 * gaussian.position.len + coefs.len, allocator); defer vars.deinit();

    for (0..gaussian.position.len) |i| {
        vars.ptr(i).* = Complex(T).init(gaussian.position[i], 0); vars.ptr(gaussian.position.len + i).* = Complex(T).init(gaussian.momentum[i], 0); vars.ptr(2 * gaussian.position.len + i).* = gaussian.gamma[i];
    } defer {
        for (0..gaussian.position.len) |i| {
            gaussian.position[i] = vars.at(i).re; gaussian.momentum[i] = vars.at(gaussian.position.len + i).re; gaussian.gamma[i] = vars.at(2 * gaussian.position.len + i);
        }
    }

    for (0..coefs.len) |i| {
        vars.ptr(3 * gaussian.position.len + i).* = coefs.at(i);
    } defer {
        for (0..coefs.len) |i| coefs.ptr(i).* = vars.at(3 * gaussian.position.len + i);
    }

    const derivative = struct {
        pub fn call(k: *ComplexVector(T), v: ComplexVector(T), params: anytype) !void {
            const g = params.gaussian; const c = v.slice(3 * g.position.len, v.len);

            for (0..g.position.len) |i| {
                g.position[i] = v.at(i).re; g.momentum[i] = v.at(g.position.len + i).re; g.gamma[i] = v.at(2 * g.position.len + i);
            }

            var kq: RealVector(T) = undefined; defer kq.deinit();
            var kp: RealVector(T) = undefined; defer kp.deinit();
            var kg: ComplexVector(T) = undefined; defer kg.deinit();
            var kc: ComplexVector(T) = undefined; defer kc.deinit();

            if (params.opt.imaginary) {

                kp = try g.momentumDerivativeImaginary(params.opt.initial_conditions.mass);
                kq = try g.positionDerivativeImaginary(params.opt.potential, c, params.opt.integration_nodes, params.time, params.opt.finite_differences_step);
                kc = try g.coefficientDerivativeImaginary(c, params.opt.potential, params.opt.integration_nodes, params.time);

            } else {

                kp = try g.momentumDerivative(params.opt.potential, c, params.opt.integration_nodes, params.time, params.opt.finite_differences_step);
                kq = try g.positionDerivative(params.opt.initial_conditions.mass);
                kc = try g.coefficientDerivative(c, params.opt.potential, params.opt.integration_nodes, params.time);
            }

            if (!params.opt.method.vMCG.frozen) {
                if (params.opt.imaginary) {
                    kg = try g.gammaDerivativeImaginary(params.opt.potential, c, params.opt.initial_conditions.mass, params.opt.integration_nodes, params.time, params.opt.finite_differences_step);
                } else {
                    kg = try g.gammaDerivative(params.opt.potential, c, params.opt.initial_conditions.mass, params.opt.integration_nodes, params.time, params.opt.finite_differences_step);
                }
            }

            for (0..g.position.len) |i| {
                k.ptr(i).* = Complex(T).init(kq.at(i), 0); k.ptr(g.position.len + i).* = Complex(T).init(kp.at(i), 0); k.ptr(2 * g.position.len + i).* = kg.at(i);
            }

            for (0..c.len) |i| k.ptr(3 * g.position.len + i).* = kc.at(i);
        }
    }.call;

    try rk.rk4(&vars, derivative, .{.gaussian = gaussian, .opt = opt, .time = time, .allocator = allocator}, opt.time_step);
}

/// Print header for iteration info.
pub fn printIterationHeader(ndim: usize, nstate: usize) !void {
    var buffer: [WRITE_BUFFER_SIZE]u8 = undefined;

    const ndim_header_width = 9 * @as(usize, @min(ndim, 3)) + 2 * (@as(usize, @min(ndim, 3)) - 1) + @as(usize, if (ndim > 3) 7 else 2);
    const nstate_header_width = 7 * @as(usize, @min(4, nstate)) + 2 * (@as(usize, @min(4, nstate)) - 1) + @as(usize, if (nstate > 4) 7 else 2);

    var writer = std.io.Writer.fixed(&buffer);

    try writer.print("\n{s:8} ", .{"ITER"});
    try writer.print("{s:12} {s:12} {s:12} ", .{"KINETIC", "POTENTIAL", "TOTAL"});
    try writer.print("{[value]s:[width]} ", .{.value = "POSITION", .width = ndim_header_width});
    try writer.print("{[value]s:[width]} ", .{.value = "MOMENTUM", .width = ndim_header_width});
    try writer.print("{[value]s:[width]} ", .{.value = "POPULATION", .width = nstate_header_width});
    try writer.print("{s:4}", .{"TIME"});

    try print("{s}\n", .{writer.buffered()});
}

/// Prints the iteration info to standard output.
pub fn printIterationInfo(comptime T: type, info: Custom(T).IterationInfo, timer: *std.time.Timer) !void {
    var buffer: [WRITE_BUFFER_SIZE]u8 = undefined;

    var writer = std.io.Writer.fixed(&buffer);

    try writer.print("{d:8} ", .{info.iteration});
    try writer.print("{d:12.6} {d:12.6} {d:12.6} ", .{info.kinetic_energy, info.potential_energy, info.kinetic_energy + info.potential_energy});

    try writer.print("[", .{});

    for (0..@min(3, info.position.len)) |i| {
        try writer.print("{d:9.4}{s}", .{info.position.at(i), if (i == info.position.len - 1) "" else ", "});
    }

    if (info.position.len > 3) try writer.print("...", .{});

    try writer.print("] [", .{});

    for (0..@min(3, info.momentum.len)) |i| {
        try writer.print("{d:9.4}{s}", .{info.momentum.at(i), if (i == info.momentum.len - 1) "" else ", "});
    }

    if (info.momentum.len > 3) try writer.print("...", .{});

    try writer.print("] [", .{});

    for (0..@min(4, info.coefs.len)) |i| {
        try writer.print("{d:7.4}{s}", .{info.coefs.at(i).magnitude() * info.coefs.at(i).magnitude(), if (i == info.coefs.len - 1) "" else ", "});
    }

    if (info.coefs.len > 4) try writer.print("...", .{});

    try writer.print("] {D}", .{timer.read()}); timer.reset();

    try print("{s}\n", .{writer.buffered()});
}

test "vMCG on Tully's First Potential" {
    const opt = Options(f64){
        .initial_conditions = .{
            .mass = &.{2000},
            .state = 1,
        },
        .method = .{
            .vMCG = .{
                .position = &.{-10},
                .gamma = &.{2},
                .momentum = &.{15}
            }
        },
        .potential = .{
            .tully_1 = TullyPotential1(f64){}
        },
        .finite_differences_step = 1e-8,
        .integration_nodes = 32,
        .iterations = 350,
        .time_step = 10
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expect(@abs(output.kinetic_energy - 0.06574552356661) < TEST_TOLERANCE);
    try std.testing.expect(@abs(output.population.at(opt.iterations, 0) - 0.53765879822849) < TEST_TOLERANCE);
    try std.testing.expect(@abs(output.population.at(opt.iterations, 1) - 0.46233829327311) < TEST_TOLERANCE);
    try std.testing.expect(@abs(output.potential_energy - 0.00075320504955) < TEST_TOLERANCE);
}
