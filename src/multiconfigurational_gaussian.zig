//! Main file for the vMCG program.

const std = @import("std");

const complex_gaussian = @import("complex_gaussian.zig");
const complex_matrix = @import("complex_matrix.zig");
const complex_runge_kutta = @import("complex_runge_kutta.zig");
const complex_vector = @import("complex_vector.zig");
const device_write = @import("device_write.zig");
const electronic_potential = @import("electronic_potential.zig");
const error_handling = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");
const grid_generator = @import("grid_generator.zig");
const harmonic_potential = @import("harmonic_potential.zig");
const hermite_quadrature_nodes = @import("hermite_quadrature_nodes.zig");
const math_functions = @import("math_functions.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const quantum_dynamics = @import("quantum_dynamics.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");
const single_set_of_mcg = @import("single_set_of_mcg.zig");
const tully_potential = @import("tully_potential.zig");

const Complex = std.math.Complex;
const ComplexGaussian = complex_gaussian.ComplexGaussian;
const ComplexMatrix = complex_matrix.ComplexMatrix;
const ComplexRungeKutta = complex_runge_kutta.ComplexRungeKutta;
const ComplexVector = complex_vector.ComplexVector;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const HarmonicPotential = harmonic_potential.HarmonicPotential;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;
const SingleSetOfMCG = single_set_of_mcg.SingleSetOfMCG;
const TullyPotential1 = tully_potential.TullyPotential1;

const exportComplexMatrixWithLinspacedLeftColumn = device_write.exportComplexMatrixWithLinspacedLeftColumn;
const exportRealMatrix = device_write.exportRealMatrix;
const exportRealMatrixWithLinspacedLeftColumn = device_write.exportRealMatrixWithLinspacedLeftColumn;
const mm = matrix_multiplication.mm;
const energySpectrum = quantum_dynamics.energySpectrum;
const positionAtRow = grid_generator.positionAtRow;
const powi = math_functions.powi;
const print = device_write.print;
const printRealMatrix = device_write.printRealMatrix;
const printComplexMatrix = device_write.printComplexMatrix;
const printComplexVector = device_write.printComplexVector;
const printRealVector = device_write.printRealVector;
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
            bf_spread: []const u32,
            state: u32
        };
        pub const LogIntervals = struct {
            iteration: u32 = 1
        };
        pub const Spectrum = struct {
            window: enum {gaussian} = .gaussian,
            padding_order: u32 = 1
        };
        pub const Write = struct {
            population: ?[]const u8 = null,
            wavefunction: ?[]const u8 = null,
            spectrum: ?[]const u8 = null,
            autocorrelation_function: ?[]const u8 = null
        };
        pub const WavefunctionGrid = struct {
            limits: []const []const T,
            points: u32
        };
        pub const Method = union(enum) {
            vMCG: struct {
                position: []const []const T,
                gamma: []const []const T,
                momentum: []const []const T,
                frozen: bool = true
            }
        };

        initial_conditions: InitialConditions,
        method: Method,
        potential: ElectronicPotential(T),

        iterations: u32,
        time_step: T,

        log_intervals: LogIntervals = .{},
        spectrum: Spectrum = .{},
        wavefunction_grid: ?WavefunctionGrid = null,
        write: Write = .{},

        finite_differences_step: T = 1e-6,
        integration_nodes: u32 = 32,
        renormalize: bool = true,
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
        pub fn deinit(self: @This(), allocator: std.mem.Allocator) void {
            self.population.deinit(allocator);
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

        /// Matrices needed for the vMCG equations of motion.
        pub const MatrixEOM = struct {
            S: ComplexMatrix(T),
            T: ComplexMatrix(T),
            V: ComplexMatrix(T),
            tau: ComplexMatrix(T)
        };

        /// Structure to hold the variables used in the equations of motion.
        pub const Variables = struct {
            all: ComplexVector(T),
            q: RealVector(T),
            p: RealVector(T),
            g: ComplexVector(T),
            c: ComplexVector(T)
        };

        /// Temporary matrices and vectors used in the propagation.
        pub const Temporaries = struct {
            nstatev: ComplexMatrix(T),
            dimv: RealVector(T),
        };
    };
}

/// Perform the simulation.
pub fn run(comptime T: type, opt: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(opt);

    if (opt.integration_nodes < 1 or opt.integration_nodes > 256) return throw(Output(T), "INTEGRATION NODES MUST BE BETWEEN 1 AND {d}", .{256});

    const ndim = opt.potential.ndim();
    const nstate = opt.potential.nstate();
    const ngauss = opt.method.vMCG.position.len;
    const nparams = 3 * ngauss * ndim + nstate * ngauss;

    var custom_potential = if (opt.potential == .custom) try opt.potential.custom.init(allocator) else null; defer if (custom_potential) |*cp| cp.deinit(allocator);
    var file_potential = if (opt.potential == .file) try opt.potential.file.init(allocator) else null; defer if (file_potential) |*fp| fp.deinit(allocator);

    var output = try Output(T).init(nstate, @intCast(opt.iterations), allocator);
    
    var mcg = try SingleSetOfMCG(T).init(
        opt.method.vMCG.position,
        opt.method.vMCG.gamma,
        opt.method.vMCG.momentum,
        opt.initial_conditions.state,
        opt.initial_conditions.bf_spread,
        nstate,
        allocator
    ); defer mcg.deinit(allocator);

    const mcg_init = try mcg.clone(allocator); defer mcg_init.deinit(allocator);

    var matrix_eom = Custom(T).MatrixEOM{
        .S = try ComplexMatrix(T).initZero(ngauss * nstate, ngauss * nstate, allocator),
        .T = try ComplexMatrix(T).initZero(ngauss * nstate, ngauss * nstate, allocator),
        .V = try ComplexMatrix(T).initZero(ngauss * nstate, ngauss * nstate, allocator),
        .tau = try ComplexMatrix(T).initZero(ngauss * nstate, ngauss * nstate, allocator)
    };

    defer inline for (std.meta.fields(@TypeOf(matrix_eom))) |field| @as(field.type, @field(matrix_eom, field.name)).deinit(allocator);

    var vars = Custom(T).Variables{
        .all = try ComplexVector(T).initZero(nparams, allocator),
        .q = try RealVector(T).initZero(ngauss * ndim, allocator),
        .p = try RealVector(T).initZero(ngauss * ndim, allocator),
        .g = try ComplexVector(T).initZero(ngauss * ndim, allocator),
        .c = try ComplexVector(T).initZero(nstate * ngauss, allocator)
    };

    defer inline for (std.meta.fields(@TypeOf(vars))) |field| @as(field.type, @field(vars, field.name)).deinit(allocator);

    var temps = Custom(T).Temporaries{
        .nstatev = try ComplexMatrix(T).initZero(nstate, nstate, allocator),
        .dimv = try RealVector(T).initZero(ndim, allocator)
    };

    defer inline for (std.meta.fields(@TypeOf(temps))) |field| @as(field.type, @field(temps, field.name)).deinit(allocator);

    var propagator = try ComplexRungeKutta(T).init(nparams, allocator); defer propagator.deinit(allocator);

    var wavefunction_dynamics: ?RealMatrix(T) = if (opt.write.wavefunction) |_|
        try initializeWavefunctionDynamicsContainer(T, opt.wavefunction_grid, nstate, opt.iterations, allocator)
    else null;

    var acf = try ComplexVector(T).init(opt.iterations + 1, allocator); defer acf.deinit(allocator);

    if (enable_printing) try printIterationHeader(ndim, mcg.coefs.len);

    var timer = try std.time.Timer.start();

    for (0..opt.iterations + 1) |i| {

        const time = @as(T, @floatFromInt(i)) * opt.time_step;

        if (i > 0) try propagate(T, &mcg, opt, &propagator, time, &vars, &matrix_eom, &temps, allocator);

        try mcg.overlap(&matrix_eom.S);
        try mcg.kinetic(&matrix_eom.T, opt.initial_conditions.mass);
        try mcg.potential(&matrix_eom.V, opt.potential, opt.integration_nodes, time, &temps.nstatev, &temps.dimv);

        if (i > 0 and opt.renormalize) try mcg.normalize(matrix_eom.S);

        acf.ptr(i).* = try mcg_init.selfOverlap(mcg, allocator);

        var coefs_adia: ?ComplexVector(T) = null; defer if (coefs_adia) |c| c.deinit(allocator);

        if (i == 0 and opt.initial_conditions.adiabatic) {
            coefs_adia = try mcg.transformedCoefs(opt.potential, time, false, allocator); try coefs_adia.?.copyTo(&mcg.coefs);
        }

        if (opt.adiabatic) {
            if (coefs_adia) |c| c.deinit(allocator); coefs_adia = try mcg.transformedCoefs(opt.potential, time, true, allocator);
        }

        const position = try mcg.coordinateExpectation(.position, matrix_eom.S, allocator); defer position.deinit(allocator);
        const momentum = try mcg.coordinateExpectation(.momentum, matrix_eom.S, allocator); defer momentum.deinit(allocator);

        const kinetic_energy = try mcg.kineticEnergy(matrix_eom.T, matrix_eom.S);
        const potential_energy = try mcg.potentialEnergy(matrix_eom.V, matrix_eom.S);

        for (0..nstate) |j| {
            output.population.ptr(i, j).* = try mcg.population(matrix_eom.S, j, if (opt.adiabatic) coefs_adia.? else mcg.coefs);
        }

        if (wavefunction_dynamics) |*wfn_container| {
            try assignWavefunction(T, wfn_container, mcg.gaussians, if (opt.adiabatic) coefs_adia.? else mcg.coefs, i);
        }

        if (i == opt.iterations) {
            output.kinetic_energy = kinetic_energy;
            output.potential_energy = potential_energy;
        }

        if (!enable_printing or (i > 0 and i % opt.log_intervals.iteration != 0)) continue;

        const info = Custom(T).IterationInfo{
            .iteration = i,
            .kinetic_energy = kinetic_energy,
            .coefs = if (opt.adiabatic) coefs_adia.? else mcg.coefs,
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

    if (opt.write.spectrum) |path| {

        if (enable_printing) try print("\nCALCULATING ENERGY SPECTRUM: ", .{});

        const spectrum = try energySpectrum(T, acf, opt.spectrum, allocator);

        if (enable_printing) try print("DONE\n", .{});

        const nyquist_frequency = std.math.pi / opt.time_step;

        try exportRealMatrixWithLinspacedLeftColumn(T, path, spectrum.asMatrix(), 0, nyquist_frequency);
    }


    if (opt.write.autocorrelation_function) |path| try exportComplexMatrixWithLinspacedLeftColumn(T, path, acf.asMatrix(), 0, end_time);
    if (opt.write.population) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.population, 0, end_time);
    if (opt.write.wavefunction) |path| try exportRealMatrix(T, path, wavefunction_dynamics.?);

    if (wavefunction_dynamics != null) wavefunction_dynamics.?.deinit(allocator);

    return output;
}

/// Assigns the single gaussian to the wavefunction dynamics container.
pub fn assignWavefunction(comptime T: type, wfn_container: *RealMatrix(T), gaussians: []const ComplexGaussian(T), coefs: ComplexVector(T), iteration: usize) !void {
    for (0..wfn_container.rows) |i| {

        for (0..coefs.len / gaussians.len) |j| {

            var value = Complex(T).init(0, 0);

            for (gaussians, 0..) |gaussian, k| {
                value = value.add(gaussian.evaluate(wfn_container.row(i).slice(0, gaussian.position.len).data).mul(coefs.at(j * gaussians.len + k)));
            }

            wfn_container.ptr(i, gaussians[0].position.len + 2 * coefs.len / gaussians.len * iteration + 2 * j + 0).* = value.re;
            wfn_container.ptr(i, gaussians[0].position.len + 2 * coefs.len / gaussians.len * iteration + 2 * j + 1).* = value.im;
        }
    }
}

/// Fix gauge for Tau matrix.
pub fn fixGaugeTauMatrix(comptime T: type, tau: *ComplexMatrix(T), S: ComplexMatrix(T)) void {
    for (0..tau.rows) |j| {

        const norm_rate = tau.at(j, j).re;

        for (0..tau.cols) |i| {
            tau.ptr(i, j).* = tau.at(i, j).sub(S.at(i, j).mul(Complex(T).init(norm_rate, 0)));
        }
    }
}

/// Initialize the container for the wavefunction dynamics.
pub fn initializeWavefunctionDynamicsContainer(comptime T: type, grid: ?Options(T).WavefunctionGrid, nstate: usize, iterations: usize, allocator: std.mem.Allocator) !RealMatrix(T) {
    if (grid == null) return throw(RealMatrix(T), "WAVEFUNCTION GRID MUST BE PROVIDED WHEN WRITING THE WAVEFUNCTION DYNAMICS", .{});

    const grid_points = powi(grid.?.points, grid.?.limits.len);

    var wavefunction_dynamics: RealMatrix(T) = try RealMatrix(T).init(grid_points, grid.?.limits.len + 2 * nstate * (iterations + 1), allocator);

    var position_at_row = try RealVector(T).init(grid.?.limits.len, allocator); defer position_at_row.deinit(allocator);

    for (0..grid_points) |i| {

        positionAtRow(T, &position_at_row, i, grid.?.limits.len, grid.?.points, grid.?.limits);

        for (0..grid.?.limits.len) |j| {
            wavefunction_dynamics.ptr(i, j).* = position_at_row.at(j);
        }
    }

    return wavefunction_dynamics;
}

/// Propagates the MCG object using the 4th order Runge-Kutta method.
pub fn propagate(comptime T: type, mcg: *SingleSetOfMCG(T), opt: Options(T), rk: *ComplexRungeKutta(T), time: T, vars: *Custom(T).Variables, matrix_eom: *Custom(T).MatrixEOM, temps: *Custom(T).Temporaries, allocator: std.mem.Allocator) !void {
    mcg.exportParameterVector(&vars.all);

    const derivative = struct {
        pub fn call(k: *ComplexVector(T), v: ComplexVector(T), params: anytype) anyerror!void {
            const mass = params.opt.initial_conditions.mass;
            const pot = params.opt.potential;
            const n_nodes = params.opt.integration_nodes;
            const fdiff_step = params.opt.finite_differences_step;

            var kq = params.kq;
            var kp = params.kp;
            var kg = params.kg;
            var kc = params.kc;

            try params.mcg.loadParameterVector(v);

            try params.mcg.overlap(&params.matrix_eom.S);
            try params.mcg.kinetic(&params.matrix_eom.T, mass);
            try params.mcg.potential(&params.matrix_eom.V, pot, n_nodes, params.time, &params.temps.nstatev, &params.temps.dimv);

            if (params.opt.imaginary) {

                try params.mcg.positionDerivativeImaginaryEhrenfest(&kq, pot, n_nodes, params.time, fdiff_step, params.allocator);
                try params.mcg.momentumDerivativeImaginaryEhrenfest(&kp, mass);

                if (!params.opt.method.vMCG.frozen) {
                    try params.mcg.gammaDerivativeImaginaryEhrenfest(&kg, pot, mass, n_nodes, params.time, fdiff_step, &params.temps.nstatev, &params.temps.dimv);
                }

                try params.mcg.overlapDiffTime(&params.matrix_eom.tau, kq, kp, kg); fixGaugeTauMatrix(T, &params.matrix_eom.tau, params.matrix_eom.S);

                switch (params.mcg.gaussians.len) {
                    1 => try params.mcg.gaussians[0].coefficientDerivativeImaginary(&kc, params.mcg.coefs, pot, n_nodes, params.time, params.allocator),
                    else => try params.mcg.coefficientDerivativeImaginary(&kc, params.matrix_eom, params.allocator)
                }

            } else {

                try params.mcg.positionDerivativeEhrenfest(&kq, mass);
                try params.mcg.momentumDerivativeEhrenfest(&kp, pot, n_nodes, params.time, fdiff_step, params.allocator);

                if (!params.opt.method.vMCG.frozen) {
                    try params.mcg.gammaDerivativeEhrenfest(&kg, pot, mass, n_nodes, params.time, fdiff_step, &params.temps.nstatev, &params.temps.dimv);
                }

                try params.mcg.overlapDiffTime(&params.matrix_eom.tau, kq, kp, kg); fixGaugeTauMatrix(T, &params.matrix_eom.tau, params.matrix_eom.S);

                switch (params.mcg.gaussians.len) {
                    1 => try params.mcg.gaussians[0].coefficientDerivative(&kc, params.mcg.coefs, pot, n_nodes, params.time, params.allocator),
                    else => try params.mcg.coefficientDerivative(&kc, params.matrix_eom, params.allocator)
                }
            }

            for (0..kq.len / params.mcg.gaussians[0].position.len) |i| for (0..params.mcg.gaussians[0].position.len) |j| {

                    k.ptr(3 * params.mcg.gaussians[0].position.len * i + j + 0 * params.mcg.gaussians[0].position.len).* = Complex(T).init(kq.at(params.mcg.gaussians[0].position.len * i + j), 0);
                    k.ptr(3 * params.mcg.gaussians[0].momentum.len * i + j + 1 * params.mcg.gaussians[0].momentum.len).* = Complex(T).init(kp.at(params.mcg.gaussians[0].momentum.len * i + j), 0);

                    k.ptr(3 * params.mcg.gaussians[0].gamma.len * i + j + 2 * params.mcg.gaussians[0].gamma.len).* = kg.at(params.mcg.gaussians[0].gamma.len * i + j);
            };

            for (0..params.mcg.coefs.len) |i| k.ptr(v.data.len - params.mcg.coefs.len + i).* = kc.at(i);
        }
    }.call;

    try rk.rk4(&vars.all, derivative, .{
        .mcg = mcg,
        .opt = opt,
        .time = time,
        .kq = vars.q,
        .kp = vars.p,
        .kg = vars.g,
        .kc = vars.c,
        .temps = temps,
        .matrix_eom = matrix_eom,
        .allocator = allocator
    }, opt.time_step);

    try mcg.loadParameterVector(vars.all);
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

test "frozen real-time vMCG on Tully's First Potential" {
    const opt = Options(f64){
        .initial_conditions = .{
            .bf_spread = &.{1},
            .mass = &.{2000},
            .state = 1
        },
        .method = .{
            .vMCG = .{
                .position = &.{&.{-10}},
                .gamma = &.{&.{2}},
                .momentum = &.{&.{15}},
                .frozen = true
            }
        },
        .potential = .{
            .tully_1 = TullyPotential1(f64){}
        },
        .finite_differences_step = 1e-8,
        .integration_nodes = 32,
        .iterations = 3500,
        .time_step = 1
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(output.kinetic_energy, 0.06574685612813, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.population.at(opt.iterations, 0), 0.53765759177167, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.population.at(opt.iterations, 1), 0.46234240819512, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.potential_energy, 0.00075315183576, TEST_TOLERANCE);
}

test "thawed real-time vMCG on Tully's First Potential" {
    const opt = Options(f64){
        .initial_conditions = .{
            .bf_spread = &.{1},
            .mass = &.{2000},
            .state = 1
        },
        .method = .{
            .vMCG = .{
                .position = &.{&.{-10}},
                .gamma = &.{&.{2}},
                .momentum = &.{&.{15}},
                .frozen = false
            }
        },
        .potential = .{
            .tully_1 = TullyPotential1(f64){}
        },
        .finite_differences_step = 1e-3,
        .integration_nodes = 32,
        .iterations = 3500,
        .time_step = 1
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(output.kinetic_energy, 0.06763705056010, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.population.at(opt.iterations, 0), 0.44325081289219, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.population.at(opt.iterations, 1), 0.55674918710781, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.potential_energy, -0.00113498374022, TEST_TOLERANCE);
}

test "thawed real-time ss-vMCG on Tully's First Potential" {
    const opt = Options(f64){
        .initial_conditions = .{
            .bf_spread = &.{1, 0},
            .mass = &.{2000},
            .state = 1
        },
        .method = .{
            .vMCG = .{
                .position = &.{&.{-10.5}, &.{-9.5}},
                .gamma = &.{&.{2}, &.{2}},
                .momentum = &.{&.{15}, &.{15}},
                .frozen = false
            }
        },
        .potential = .{
            .tully_1 = TullyPotential1(f64){}
        },
        .finite_differences_step = 1e-3,
        .integration_nodes = 32,
        .iterations = 3500,
        .time_step = 1
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(output.kinetic_energy, 0.06706660848155, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.population.at(opt.iterations, 0), 0.44575305875514, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.population.at(opt.iterations, 1), 0.55424694124486, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.potential_energy, -0.00108493882998, TEST_TOLERANCE);
}

test "frozen real-time ss-vMCG on Tully's First Potential" {
    const opt = Options(f64){
        .initial_conditions = .{
            .bf_spread = &.{1, 0},
            .mass = &.{2000},
            .state = 1
        },
        .method = .{
            .vMCG = .{
                .position = &.{&.{-10.5}, &.{-9.5}},
                .gamma = &.{&.{2}, &.{2}},
                .momentum = &.{&.{15}, &.{15}},
                .frozen = true
            }
        },
        .potential = .{
            .tully_1 = TullyPotential1(f64){}
        },
        .finite_differences_step = 1e-8,
        .integration_nodes = 32,
        .iterations = 3500,
        .time_step = 1
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(output.kinetic_energy, 0.07135428435411, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.population.at(opt.iterations, 0), 0.55681879770101, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.population.at(opt.iterations, 1), 0.44318120229899, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.potential_energy, 0.00113637595402, TEST_TOLERANCE);
}
