//! Code to propagate classical particles.

const std = @import("std");

const bias_potential = @import("bias_potential.zig");
const classical_particle = @import("classical_particle.zig");
const complex_runge_kutta = @import("complex_runge_kutta.zig");
const complex_vector = @import("complex_vector.zig");
const derivative_coupling = @import("derivative_coupling.zig");
const device_write = @import("device_write.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const electronic_potential = @import("electronic_potential.zig");
const error_handling = @import("error_handling.zig");
const fewest_switches = @import("fewest_switches.zig");
const global_variables = @import("global_variables.zig");
const hammes_schiffer_tully = @import("hammes_schiffer_tully.zig");
const harmonic_potential = @import("harmonic_potential.zig");
const landau_zener = @import("landau_zener.zig");
const mapping_approach = @import("mapping_approach.zig");
const math_functions = @import("math_functions.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const nonadiabatic_coupling_vector = @import("nonadiabatic_coupling_vector.zig");
const norm_preserving_interpolation = @import("norm_preserving_interpolation.zig");
const object_array = @import("object_array.zig");
const parallel_tools = @import("parallel_tools.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");
const ring_buffer = @import("ring_buffer.zig");
const surface_hopping_algorithm = @import("surface_hopping_algorithm.zig");
const tully_potential = @import("tully_potential.zig");

const BiasPotential = bias_potential.BiasPotential;
const ClassicalParticle = classical_particle.ClassicalParticle;
const Complex = std.math.complex.Complex;
const ComplexRungeKutta = complex_runge_kutta.ComplexRungeKutta;
const ComplexVector = complex_vector.ComplexVector;
const DerivativeCoupling = derivative_coupling.DerivativeCoupling;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const FewestSwitches = fewest_switches.FewestSwitches;
const HarmonicPotential = harmonic_potential.HarmonicPotential;
const LandauZener = landau_zener.LandauZener;
const MappingApproach = mapping_approach.MappingApproach;
const NormPreservingInterpolation = norm_preserving_interpolation.NormPreservingInterpolation;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;
const RingBufferArray = object_array.RingBufferArray;
const RealMatrixArray = object_array.RealMatrixArray;
const RealVectorArray = object_array.RealVectorArray;
const SurfaceHoppingAlgorithm = surface_hopping_algorithm.SurfaceHoppingAlgorithm;
const TullyPotential1 = tully_potential.TullyPotential1;

const binomialConfInt = math_functions.binomialConfInt;
const checkParallelError = parallel_tools.checkParallelError;
const exportRealMatrix = device_write.exportRealMatrix;
const exportRealMatrixWithLinspacedLeftColumn = device_write.exportRealMatrixWithLinspacedLeftColumn;
const fixGauge = eigenproblem_solver.fixGauge;
const mm = matrix_multiplication.mm;
const print = device_write.print;
const printJson = device_write.printJson;
const printRealMatrix = device_write.printRealMatrix;
const throw = error_handling.throw;

const MAX_POOL_SIZE = global_variables.MAX_POOL_SIZE;
const WRITE_BUFFER_SIZE = global_variables.WRITE_BUFFER_SIZE;

var PARALLEL_ERROR = &global_variables.PARALLEL_ERROR;

/// Classical dynamics options struct.
pub fn Options(comptime T: type) type {
    return struct {
        pub const InitialConditions = struct {
            mass: []const T,
            momentum_mean: []const T,
            momentum_std: []const T,
            position_mean: []const T,
            position_std: []const T,
            state: u32
        };
        pub const LogIntervals = struct {
            trajectory: u32 = 1,
            iteration: u32 = 1
        };
        pub const Write = struct {
            kinetic_energy_mean: ?[]const u8 = null,
            momentum_mean: ?[]const u8 = null,
            population_mean: ?[]const u8 = null,
            position_mean: ?[]const u8 = null,
            potential_energy_mean: ?[]const u8 = null,
            time_derivative_coupling_mean: ?[]const u8 = null,
            total_energy_mean: ?[]const u8 = null
        };

        potential: ElectronicPotential(T),
        initial_conditions: InitialConditions,

        trajectories: u32,
        iterations: u32,
        time_step: T,

        bias: ?BiasPotential(T) = null,
        derivative_coupling: ?DerivativeCoupling(T) = null,
        log_intervals: LogIntervals = .{},
        surface_hopping: ?SurfaceHoppingAlgorithm(T) = null,
        write: Write = .{},

        finite_differences_step: T = 1e-8,
        seed: u32 = 0,
        nthread: u32 = 1
    };
}

/// Structure that hold the output of the simulation.
pub fn Output(comptime T: type) type {
    return struct {
        kinetic_energy_mean: RealVector(T),
        momentum_mean: RealMatrix(T),
        population_mean: RealMatrix(T),
        position_mean: RealMatrix(T),
        potential_energy_mean: RealVector(T),
        time_derivative_coupling_mean: RealMatrix(T),
        total_energy_mean: RealVector(T),

        /// Allocate the output structure.
        pub fn init(nstate: usize, ndim: usize, iterations: usize, allocator: std.mem.Allocator) !@This() {
            return @This(){
                .kinetic_energy_mean = try RealVector(T).initZero(iterations + 1, allocator),
                .momentum_mean = try RealMatrix(T).initZero(iterations + 1, ndim, allocator),
                .population_mean = try RealMatrix(T).initZero(iterations + 1, nstate, allocator),
                .position_mean = try RealMatrix(T).initZero(iterations + 1, ndim, allocator),
                .potential_energy_mean = try RealVector(T).initZero(iterations + 1, allocator),
                .time_derivative_coupling_mean = try RealMatrix(T).initZero(iterations + 1, nstate * nstate, allocator),
                .total_energy_mean = try RealVector(T).initZero(iterations + 1, allocator)
            };
        }

        /// Free the memory allocated for the output structure.
        pub fn deinit(self: @This()) void {
            self.kinetic_energy_mean.deinit();
            self.momentum_mean.deinit();
            self.population_mean.deinit();
            self.position_mean.deinit();
            self.potential_energy_mean.deinit();
            self.time_derivative_coupling_mean.deinit();
            self.total_energy_mean.deinit();
        }
    };
}

/// Container for custom structs related to the classical dynamics.
pub fn Custom(comptime T: type) type {
    return struct {

        /// Structure to hold information about each iteration used for logging.
        pub const IterationInfo = struct {
            iteration: usize,
            kinetic_energy: T,
            potential_energy: T,
            state: usize,
            system: ClassicalParticle(T),
            time: T,
            trajectory: usize,
            amplitudes: ComplexVector(T)
        };

        /// Structure to hold a single trajectory output.
        pub const TrajectoryOutput = struct {
            kinetic_energy: RealVector(T),
            momentum: RealMatrix(T),
            population: RealMatrix(T),
            position: RealMatrix(T),
            potential_energy: RealVector(T),
            time_derivative_coupling: RealMatrix(T),
            total_energy: RealVector(T),

            /// Allocate the trajectory output structure.
            pub fn init(nstate: usize, ndim: usize, iterations: usize, allocator: std.mem.Allocator) !@This() {
                return @This(){
                    .kinetic_energy = try RealVector(T).initZero(iterations + 1, allocator),
                    .momentum = try RealMatrix(T).initZero(iterations + 1, ndim, allocator),
                    .population = try RealMatrix(T).initZero(iterations + 1, nstate, allocator),
                    .position = try RealMatrix(T).initZero(iterations + 1, ndim, allocator),
                    .potential_energy = try RealVector(T).initZero(iterations + 1, allocator),
                    .time_derivative_coupling = try RealMatrix(T).initZero(iterations + 1, nstate * nstate, allocator),
                    .total_energy = try RealVector(T).initZero(iterations + 1, allocator),
                };
            }

            /// Free the memory allocated for the trajectory output structure.
            pub fn deinit(self: @This()) void {
                self.kinetic_energy.deinit();
                self.momentum.deinit();
                self.population.deinit();
                self.position.deinit();
                self.potential_energy.deinit();
                self.time_derivative_coupling.deinit();
                self.total_energy.deinit();
            }
        };
    };
}

/// Run classical dynamics simulation.
pub fn run(comptime T: type, opt: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(opt);

    const ndim = opt.potential.ndim();
    const nstate = opt.potential.nstate();

    var custom_potential = if (opt.potential == .custom) try opt.potential.custom.init(allocator) else null; defer if (custom_potential) |*cp| cp.deinit();
    var file_potential = if (opt.potential == .file) try opt.potential.file.init(allocator) else null; defer if (file_potential) |*fp| fp.deinit();

    var output = try Output(T).init(nstate, ndim, opt.iterations, allocator);

    const parallel_results = .{
        .kinetic_energy_mean = try RealVectorArray(T).init(opt.nthread, .{.rows = opt.iterations + 1}, allocator),
        .momentum_mean = try RealMatrixArray(T).init(opt.nthread, .{.rows = opt.iterations + 1, .cols = ndim}, allocator),
        .population_mean = try RealMatrixArray(T).init(opt.nthread, .{.rows = opt.iterations + 1, .cols = nstate}, allocator),
        .position_mean = try RealMatrixArray(T).init(opt.nthread, .{.rows = opt.iterations + 1, .cols = ndim}, allocator),
        .potential_energy_mean = try RealVectorArray(T).init(opt.nthread, .{.rows = opt.iterations + 1}, allocator),
        .time_derivative_coupling_mean = try RealMatrixArray(T).init(opt.nthread, .{.rows = opt.iterations + 1, .cols = nstate * nstate}, allocator),
        .total_energy_mean = try RealVectorArray(T).init(opt.nthread, .{.rows = opt.iterations + 1}, allocator),
    };

    defer inline for (std.meta.fields(@TypeOf(parallel_results))) |field| @as(field.type, @field(parallel_results, field.name)).deinit();

    var split_mix = std.Random.SplitMix64.init(opt.seed);

    if (enable_printing) try printIterationHeader(T, ndim, nstate, opt.surface_hopping);

    var pool: std.Thread.Pool = undefined; var wait: std.Thread.WaitGroup = undefined;

    try pool.init(.{.n_jobs = opt.nthread, .track_ids = true, .allocator = allocator});

    for (0..(opt.trajectories + MAX_POOL_SIZE - 1) / MAX_POOL_SIZE) |i| {

        wait.reset();

        for (0..@min(opt.trajectories - i * MAX_POOL_SIZE, MAX_POOL_SIZE)) |j| {

            const rng = std.Random.DefaultPrng.init(split_mix.next());

            const params = .{opt, i * MAX_POOL_SIZE + j, enable_printing, rng, allocator};

            if (opt.nthread == 1) {
                runTrajectoryParallel(1, T, parallel_results, params);
            } else {
                pool.spawnWgId(&wait, runTrajectoryParallel, .{T, parallel_results, params});
            }
        }

        wait.wait();
    }

    pool.deinit(); try checkParallelError();

    try finalizeOutput(T, &output, opt, parallel_results);

    if (enable_printing) try printFinalDetails(T, opt, output);

    return output;
}

/// Run a single trajectory.
pub fn runTrajectory(comptime T: type, opt: Options(T), system: *ClassicalParticle(T), index: usize, enable_printing: bool, allocator: std.mem.Allocator) !Custom(T).TrajectoryOutput {
    const nstate = opt.potential.nstate();
    const ndim = opt.potential.ndim();

    var split_mix = std.Random.SplitMix64.init(opt.seed + index); var rng = std.Random.DefaultPrng.init(split_mix.next()); var random = rng.random();

    var output = try Custom(T).TrajectoryOutput.init(nstate, ndim, opt.iterations, allocator);

    var current_state: usize = @intCast(opt.initial_conditions.state);

    if (current_state >= nstate) return throw(Custom(T).TrajectoryOutput, "ACTIVE STATE MUST NOT BE HIGHER THAN THE TOTAL NUMBER OF STATES", .{});

    var diabatic_potential = try RealMatrix(T).init(nstate, nstate, allocator); defer diabatic_potential.deinit();
    var adiabatic_potential = try RealMatrix(T).init(nstate, nstate, allocator); defer adiabatic_potential.deinit();
    var adiabatic_eigenvectors = try RealMatrix(T).init(nstate, nstate, allocator); defer adiabatic_eigenvectors.deinit();
    var previous_eigenvectors = try RealMatrix(T).init(nstate, nstate, allocator); defer previous_eigenvectors.deinit();
    var eigenvector_overlap = try RealMatrix(T).init(nstate, nstate, allocator); defer eigenvector_overlap.deinit();
    var time_derivative_coupling = try RealMatrix(T).initZero(nstate, nstate, allocator); defer time_derivative_coupling.deinit();

    var energy_gaps = try RingBufferArray(T).init(nstate * (nstate - 1) / 2, .{.max_len = 5}, allocator); defer energy_gaps.deinit();
    var jump_probabilities = try RealVector(T).init(nstate, allocator); defer jump_probabilities.deinit();
    var runge_kutta_solver = try ComplexRungeKutta(T).init(nstate, allocator); defer runge_kutta_solver.deinit();
    var amplitudes = try ComplexVector(T).initZero(nstate, allocator); defer amplitudes.deinit();
    var bloch_vector = try RealVector(T).initZero(3, allocator); defer bloch_vector.deinit();

    amplitudes.ptr(current_state).* = Complex(T).init(1, 0); bloch_vector.ptr(2).* = if (current_state == 1) 1 else -1;

    if (opt.surface_hopping != null and opt.surface_hopping.? == .mapping_approach) {

        const phi = 2 * std.math.pi * random.float(T);

        const cos_theta = bloch_vector.at(2) * std.math.sqrt(random.float(T));
        const sin_theta = std.math.sqrt(1 - cos_theta * cos_theta);

        bloch_vector.ptr(0).* = sin_theta * std.math.cos(phi);
        bloch_vector.ptr(1).* = sin_theta * std.math.sin(phi);
        bloch_vector.ptr(2).* = cos_theta;
    }

    const fs_parameters: fewest_switches.Parameters(T) = .{
        .adiabatic_potential = adiabatic_potential,
        .amplitudes = &amplitudes,
        .derivative_coupling = time_derivative_coupling,
        .runge_kutta = runge_kutta_solver,
        .time_step = opt.time_step
    };

    const lz_parameters: landau_zener.Parameters(T) = .{
        .energy_gaps = energy_gaps,
        .time_step = opt.time_step
    };

    const ma_parameters: mapping_approach.Parameters(T) = .{
        .adiabatic_potential = adiabatic_potential,
        .bloch_vector = &bloch_vector,
        .derivative_coupling = time_derivative_coupling,
        .time_step = opt.time_step
    };

    const surface_hopping_parameters: surface_hopping_algorithm.Parameters(T) = .{
        .fs_parameters = fs_parameters,
        .lz_parameters = lz_parameters,
        .ma_parameters = ma_parameters
    };

    const hst_parameters: hammes_schiffer_tully.Parameters(T) = .{
        .eigenvector_overlap = eigenvector_overlap,
        .time_step = opt.time_step
    };

    const nacv_parameters: nonadiabatic_coupling_vector.Parameters(T) = .{
        .adiabatic_potential = adiabatic_potential,
        .diabatic_potential = diabatic_potential,
        .adiabatic_eigenvectors = adiabatic_eigenvectors,
        .electronic_potential = opt.potential,
        .position = system.position,
        .velocity = system.velocity,
        .time = undefined
    };

    const npi_parameters: norm_preserving_interpolation.Parameters(T) = .{
        .eigenvector_overlap = eigenvector_overlap,
        .time_step = opt.time_step
    };

    const derivative_coupling_parameters: derivative_coupling.Parameters(T) = .{
        .hst_parameters = hst_parameters,
        .nacv_parameters = nacv_parameters,
        .npi_parameters = npi_parameters
    };

    var timer = try std.time.Timer.start();

    for (0..opt.iterations + 1) |i| {

        const time = @as(T, @floatFromInt(i)) * opt.time_step; @constCast(&nacv_parameters.time).* = time;

        try adiabatic_eigenvectors.copyTo(&previous_eigenvectors);

        if (i > 0) {
            try system.propagateVelocityVerlet(opt.potential, &adiabatic_potential, time, current_state, opt.time_step, opt.finite_differences_step, opt.bias);
        }

        try opt.potential.evaluateEigensystem(&diabatic_potential, &adiabatic_potential, &adiabatic_eigenvectors, system.position, time);

        if (i > 0) {

            try fixGauge(T, &adiabatic_eigenvectors, previous_eigenvectors);
            try mm(T, &eigenvector_overlap, previous_eigenvectors, true, adiabatic_eigenvectors, false);

            if (opt.derivative_coupling) |time_derivative_coupling_algorithm| {
                try time_derivative_coupling_algorithm.evaluate(&time_derivative_coupling, derivative_coupling_parameters);
            }
        }

        for (0..nstate) |j| for (j + 1..nstate) |k| {
            energy_gaps.ptr(j * (2 * nstate - j - 1) / 2 + (k - j - 1)).append(adiabatic_potential.at(k, k) - adiabatic_potential.at(j, j));
        };

        if (opt.surface_hopping) |algorithm| if (i > 1) {
            current_state = try algorithm.jump(system, &jump_probabilities, surface_hopping_parameters, adiabatic_potential, current_state, &random);
        };

        const kinetic_energy = system.kineticEnergy();
        const potential_energy = adiabatic_potential.at(current_state, current_state);

        output.kinetic_energy.ptr(i).* = kinetic_energy;
        output.population.ptr(i, current_state).* = 1;
        output.potential_energy.ptr(i).* = potential_energy;
        output.total_energy.ptr(i).* = kinetic_energy + potential_energy;

        for (0..nstate * nstate) |j| output.time_derivative_coupling.ptr(i, j).* = time_derivative_coupling.at(j / nstate, j % nstate);

        for (0..ndim) |j| {
            output.position.ptr(i, j).* = system.position.at(j); output.momentum.ptr(i, j).* = system.velocity.at(j) * system.masses.at(j);
        }

        if (!enable_printing or (index > 0 and (index + 1) % opt.log_intervals.trajectory != 0) or (i > 0 and i % opt.log_intervals.iteration != 0)) continue;

        const iteration_info = Custom(T).IterationInfo{
            .iteration = i,
            .kinetic_energy = kinetic_energy,
            .potential_energy = potential_energy,
            .state = current_state,
            .system = system.*,
            .time = time,
            .trajectory = index,
            .amplitudes = amplitudes
        };

        try printIterationInfo(T, iteration_info, opt.surface_hopping, &timer);
    }

    return output;
}

/// Parallel function to run a trajectory.
pub fn runTrajectoryParallel(id: usize, comptime T: type, results: anytype, params: anytype) void {
    var system = ClassicalParticle(T).initZero(params[0].potential.ndim(), params[0].initial_conditions.mass, params[4]) catch |e| {
        if (PARALLEL_ERROR.* == null) PARALLEL_ERROR.* = e; return;
    }; defer system.deinit();

    var rng = params[3]; var random = rng.random();

    sampleInitialConditions(T, &system, params[0].initial_conditions, &random) catch |e| {
        if (PARALLEL_ERROR.* == null) PARALLEL_ERROR.* = e; return;
    };

    const trajectory_output = runTrajectory(T, params[0], &system, params[1], params[2], params[4]) catch |e| {
        if (PARALLEL_ERROR.* == null) PARALLEL_ERROR.* = e; return;
    }; defer trajectory_output.deinit();

    results.kinetic_energy_mean.ptr(id - 1).add(trajectory_output.kinetic_energy) catch |e| {
        if (PARALLEL_ERROR.* == null) PARALLEL_ERROR.* = e; return;
    };

    results.momentum_mean.ptr(id - 1).add(trajectory_output.momentum) catch |e| {
        if (PARALLEL_ERROR.* == null) PARALLEL_ERROR.* = e; return;
    };

    results.population_mean.ptr(id - 1).add(trajectory_output.population) catch |e| {
        if (PARALLEL_ERROR.* == null) PARALLEL_ERROR.* = e; return;
    };

    results.position_mean.ptr(id - 1).add(trajectory_output.position) catch |e| {
        if (PARALLEL_ERROR.* == null) PARALLEL_ERROR.* = e; return;
    };

    results.potential_energy_mean.ptr(id - 1).add(trajectory_output.potential_energy) catch |e| {
        if (PARALLEL_ERROR.* == null) PARALLEL_ERROR.* = e; return;
    };

    results.time_derivative_coupling_mean.ptr(id - 1).add(trajectory_output.time_derivative_coupling) catch |e| {
        if (PARALLEL_ERROR.* == null) PARALLEL_ERROR.* = e; return;
    };

    results.total_energy_mean.ptr(id - 1).add(trajectory_output.total_energy) catch |e| {
        if (PARALLEL_ERROR.* == null) PARALLEL_ERROR.* = e; return;
    };
}

/// Print header for iteration info.
pub fn printIterationHeader(comptime T: type, ndim: usize, nstate: usize, surface_hopping: ?SurfaceHoppingAlgorithm(T)) !void {
    var buffer: [WRITE_BUFFER_SIZE]u8 = undefined;

    const ndim_header_width = 9 * @as(usize, @min(ndim, 3)) + 2 * (@as(usize, @min(ndim, 3)) - 1) + @as(usize, if (ndim > 3) 7 else 2);
    const nstate_header_width = 7 * @as(usize, @min(4, nstate)) + 2 * (@as(usize, @min(4, nstate)) - 1) + @as(usize, if (nstate > 4) 7 else 2);

    var writer = std.io.Writer.fixed(&buffer);

    try writer.print("\n{s:8} {s:8} ", .{"TRAJ", "ITER"});
    try writer.print("{s:12} {s:12} {s:12} ", .{"KINETIC", "POTENTIAL", "TOTAL"});
    try writer.print("{s:5} ", .{"STATE"});
    try writer.print("{[value]s:[width]} ", .{.value = "POSITION", .width = ndim_header_width});
    try writer.print("{[value]s:[width]} ", .{.value = "MOMENTUM", .width = ndim_header_width});

    if (surface_hopping) |algorithm| switch (algorithm) {
        .fewest_switches => try writer.print("{[value]s:[width]} ", .{.value = "|COEFS|^2", .width = nstate_header_width}),
        else => {}
    };

    try writer.print("{s:4}", .{"TIME"});

    try print("{s}\n", .{writer.buffered()});
}

/// Add the partial results to the output vectors and export them if requested.
pub fn finalizeOutput(comptime T: type, output: *Output(T), opt: Options(T), parallel_results: anytype) !void {
    const output_kinetic_energy_mean = parallel_results.kinetic_energy_mean;
    const output_momentum_mean = parallel_results.momentum_mean;
    const output_population_mean = parallel_results.population_mean;
    const output_position_mean = parallel_results.position_mean;
    const output_potential_energy_mean = parallel_results.potential_energy_mean;
    const output_time_derivative_coupling_mean = parallel_results.time_derivative_coupling_mean;
    const output_total_energy_mean = parallel_results.total_energy_mean;

    for (0..opt.nthread) |i| {
        try output.kinetic_energy_mean.add(output_kinetic_energy_mean.at(i));
        try output.momentum_mean.add(output_momentum_mean.at(i));
        try output.population_mean.add(output_population_mean.at(i));
        try output.position_mean.add(output_position_mean.at(i));
        try output.potential_energy_mean.add(output_potential_energy_mean.at(i));
        try output.time_derivative_coupling_mean.add(output_time_derivative_coupling_mean.at(i));
        try output.total_energy_mean.add(output_total_energy_mean.at(i));
    }

    output.kinetic_energy_mean.divs(@as(T, @floatFromInt(opt.trajectories)));
    output.momentum_mean.divs(@as(T, @floatFromInt(opt.trajectories)));
    output.population_mean.divs(@as(T, @floatFromInt(opt.trajectories)));
    output.position_mean.divs(@as(T, @floatFromInt(opt.trajectories)));
    output.potential_energy_mean.divs(@as(T, @floatFromInt(opt.trajectories)));
    output.time_derivative_coupling_mean.divs(@as(T, @floatFromInt(opt.trajectories)));
    output.total_energy_mean.divs(@as(T, @floatFromInt(opt.trajectories)));

    const end_time = @as(T, @floatFromInt(opt.iterations)) * opt.time_step;

    if (opt.write.kinetic_energy_mean) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.kinetic_energy_mean.asMatrix(), 0, end_time, opt.iterations + 1);
    if (opt.write.momentum_mean) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.momentum_mean, 0, end_time, opt.iterations + 1);
    if (opt.write.population_mean) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.population_mean, 0, end_time, opt.iterations + 1);
    if (opt.write.position_mean) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.position_mean, 0, end_time, opt.iterations + 1);
    if (opt.write.potential_energy_mean) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.potential_energy_mean.asMatrix(), 0, end_time, opt.iterations + 1);
    if (opt.write.time_derivative_coupling_mean) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.time_derivative_coupling_mean, 0, end_time, opt.iterations + 1);
    if (opt.write.total_energy_mean) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.total_energy_mean.asMatrix(), 0, end_time, opt.iterations + 1);
}

/// Initialize random number generators for parallel environment.
pub fn initRandomParallel(nthread: u32, seed: u32, allocator: std.mem.Allocator) !struct{rng: []std.Random.DefaultPrng, random: []std.Random} {
    var split_mix = std.Random.SplitMix64.init(seed);

    var rng_parallel = try allocator.alloc(std.Random.DefaultPrng, nthread);
    var random_parallel = try allocator.alloc(std.Random, nthread);

    for (0..nthread) |i| {
        rng_parallel[i] = std.Random.DefaultPrng.init(split_mix.next());
        random_parallel[i] = rng_parallel[i].random();
    }

    return .{.rng = rng_parallel, .random = random_parallel};
}

/// Prints the fine details after simulation.
pub fn printFinalDetails(comptime T: type, opt: Options(T), output: Output(T)) !void {
    try print("\nFINAL AVERAGED TOTAL ENERGY: {d:.14}\n", .{output.total_energy_mean.at(opt.iterations)});

    for (0..output.population_mean.cols) |i| {

        const population_error = binomialConfInt(output.population_mean.at(opt.iterations, i), opt.trajectories);

        const print_payload = .{if (i == 0) "\n" else "", i, output.population_mean.at(opt.iterations, i), population_error};

        try print("{s}FINAL POPULATION OF STATE {d:2}: {d:.6} ± {:.6}\n", print_payload);
    }
}

/// Prints the iteration info to standard output.
pub fn printIterationInfo(comptime T: type, info: Custom(T).IterationInfo, surface_hopping: ?SurfaceHoppingAlgorithm(T), timer: *std.time.Timer) !void {
    var buffer: [WRITE_BUFFER_SIZE]u8 = undefined;

    var writer = std.io.Writer.fixed(&buffer);

    try writer.print("{d:8} {d:8} ", .{info.trajectory + 1, info.iteration});
    try writer.print("{d:12.6} {d:12.6} {d:12.6} ", .{info.kinetic_energy, info.potential_energy, info.kinetic_energy + info.potential_energy});
    try writer.print("{d:5} ", .{info.state});

    try writer.print("[", .{});

    for (0..@min(3, info.system.ndim)) |i| {
        try writer.print("{d:9.4}{s}", .{info.system.position.at(i), if (i == info.system.ndim - 1) "" else ", "});
    }

    if (info.system.ndim > 3) try writer.print("...", .{});
    
    try writer.print("] [", .{});

    for (0..@min(3, info.system.ndim)) |i| {
        try writer.print("{d:9.4}{s}", .{info.system.velocity.at(i) * info.system.masses.at(i), if (i == info.system.ndim - 1) "" else ", "});
    }

    if (info.system.ndim > 3) try writer.print("...", .{});

    if (surface_hopping) |algorithm| {

        if (algorithm == .fewest_switches) try writer.print("] [", .{});

        if (algorithm == .fewest_switches) {

            for (0..@min(4, info.amplitudes.len)) |i| {
                try writer.print("{d:7.4}{s}", .{std.math.pow(T, info.amplitudes.at(i).magnitude(), 2), if (i == info.amplitudes.len - 1) "" else ", "});
            }

            if (info.amplitudes.len > 4) try writer.print("...", .{});
        }
    }

    try writer.print("] {D}", .{timer.read()}); timer.reset();

    try print("{s}\n", .{writer.buffered()});
}

/// Samples the initial conditions.
pub fn sampleInitialConditions(comptime T: type, system: *ClassicalParticle(T), initial_conditions: Options(T).InitialConditions, random: *std.Random) !void {
    try system.setPositionRandn(initial_conditions.position_mean, initial_conditions.position_std, random);
    try system.setMomentumRandn(initial_conditions.momentum_mean, initial_conditions.momentum_std, random);
}

test "Fewest Switches Surface Hopping on Tully's First Potential" {
    const opt = Options(f64){
        .derivative_coupling = .{
            .npi = NormPreservingInterpolation(f64){}
        },
        .initial_conditions = .{
            .mass = &.{2000},
            .momentum_mean = &.{15},
            .momentum_std = &.{1},
            .position_mean = &.{-10},
            .position_std = &.{0.5},
            .state = 1
        },
        .potential = .{
            .tully_1 = TullyPotential1(f64){}
        },
        .surface_hopping = .{
            .fewest_switches = FewestSwitches(f64){}
        },
        .iterations = 3500,
        .time_step = 1,
        .trajectories = 1000
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expectEqual(output.population_mean.at(opt.iterations, 0), 0.371);
    try std.testing.expectEqual(output.population_mean.at(opt.iterations, 1), 0.629);
}

test "Landau-Lener Surface Hopping on Tully's First Potential" {
    const opt = Options(f64){
        .initial_conditions = .{
            .mass = &.{2000},
            .momentum_mean = &.{15},
            .momentum_std = &.{1},
            .position_mean = &.{-10},
            .position_std = &.{0.5},
            .state = 1
        },
        .potential = .{
            .tully_1 = TullyPotential1(f64){}
        },
        .surface_hopping = .{
            .landau_zener = LandauZener(f64){}
        },
        .iterations = 3500,
        .time_step = 1,
        .trajectories = 1000
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expectEqual(output.population_mean.at(opt.iterations, 0), 0.498);
    try std.testing.expectEqual(output.population_mean.at(opt.iterations, 1), 0.502);
}

test "Mapping Approach to Surface Hopping on Tully's First Potential" {
    const opt = Options(f64){
        .derivative_coupling = .{
            .npi = NormPreservingInterpolation(f64){}
        },
        .initial_conditions = .{
            .mass = &.{2000},
            .momentum_mean = &.{15},
            .momentum_std = &.{1},
            .position_mean = &.{-10},
            .position_std = &.{0.5},
            .state = 1
        },
        .potential = .{
            .tully_1 = TullyPotential1(f64){}
        },
        .surface_hopping = .{
            .mapping_approach = MappingApproach(f64){}
        },
        .iterations = 3500,
        .time_step = 1,
        .trajectories = 1000
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit();

    try std.testing.expectEqual(output.population_mean.at(opt.iterations, 0), 0.385);
    try std.testing.expectEqual(output.population_mean.at(opt.iterations, 1), 0.615);
}
