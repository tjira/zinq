//! Code to propagate classical particles.

const std = @import("std");

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
const harmonic_potential = @import("harmonic_potential.zig");
const landau_zener = @import("landau_zener.zig");
const math_functions = @import("math_functions.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const norm_preserving_interpolation = @import("norm_preserving_interpolation.zig");
const object_array = @import("object_array.zig");
const parallel_tools = @import("parallel_tools.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");
const ring_buffer = @import("ring_buffer.zig");
const surface_hopping_algorithm = @import("surface_hopping_algorithm.zig");
const tully_potential = @import("tully_potential.zig");

const ClassicalParticle = classical_particle.ClassicalParticle;
const Complex = std.math.complex.Complex;
const ComplexRungeKutta = complex_runge_kutta.ComplexRungeKutta;
const ComplexVector = complex_vector.ComplexVector;
const DerivativeCoupling = derivative_coupling.DerivativeCoupling;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const FewestSwitches = fewest_switches.FewestSwitches;
const HarmonicPotential = harmonic_potential.HarmonicPotential;
const LandauZener = landau_zener.LandauZener;
const NormPreservingInterpolation = norm_preserving_interpolation.NormPreservingInterpolation;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;
const RingBufferArray = object_array.RingBufferArray;
const RealMatrixArray = object_array.RealMatrixArray;
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
            population_mean: ?[]const u8 = null
        };

        potential: ElectronicPotential(T),
        initial_conditions: InitialConditions,

        trajectories: u32,
        iterations: u32,
        time_step: T,

        derivative_coupling: ?DerivativeCoupling(T) = null,
        log_intervals: LogIntervals = .{},
        surface_hopping: ?SurfaceHoppingAlgorithm(T) = null,
        write: Write = .{},

        seed: u32 = 0,
        nthread: u32 = 1
    };
}

/// Structure that hold the output of the simulation.
pub fn Output(comptime T: type) type {
    return struct {
        population_mean: RealMatrix(T),

        /// Allocate the output structure.
        pub fn init(nstate: usize, iterations: usize, allocator: std.mem.Allocator) !@This() {
            return @This(){
                .population_mean = try RealMatrix(T).initZero(iterations + 1, nstate, allocator)
            };
        }

        /// Free the memory allocated for the output structure.
        pub fn deinit(self: @This()) void {
            self.population_mean.deinit();
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
            population: RealMatrix(T),

            /// Allocate the trajectory output structure.
            pub fn init(nstate: usize, iterations: usize, allocator: std.mem.Allocator) !@This() {
                return @This(){
                    .population = try RealMatrix(T).initZero(iterations + 1, nstate, allocator)
                };
            }

            /// Free the memory allocated for the trajectory output structure.
            pub fn deinit(self: @This()) void {
                self.population.deinit();
            }
        };
    };
}

/// Run classical dynamics simulation.
pub fn run(comptime T: type, opt: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(opt);

    var potential = opt.potential;
    const ndim = potential.ndim();
    const nstate = potential.nstate();

    const file_potential = if (potential == .file) try potential.file.init(allocator) else null; defer if (file_potential) |U| U.deinit();

    var output = try Output(T).init(nstate, opt.iterations, allocator);

    var output_population_mean = try RealMatrixArray(T).init(opt.nthread, .{.rows = opt.iterations + 1, .cols = nstate}, allocator); defer output_population_mean.deinit();

    var split_mix = std.Random.SplitMix64.init(opt.seed);

    if (enable_printing) try printIterationHeader(T, ndim, nstate, opt.surface_hopping);

    var pool: std.Thread.Pool = undefined; try pool.init(.{.n_jobs = opt.nthread, .track_ids = true, .allocator = allocator});
    var wait: std.Thread.WaitGroup = undefined; wait.reset();

    for (0..(opt.trajectories + MAX_POOL_SIZE - 1) / MAX_POOL_SIZE) |i| {

        wait.reset();

        for (0..@min(opt.trajectories - i * MAX_POOL_SIZE, MAX_POOL_SIZE)) |j| {

            const rng = std.Random.DefaultPrng.init(split_mix.next());
            var opt_copy = opt; opt_copy.potential = potential;

            const params = .{opt_copy, i * MAX_POOL_SIZE + j, enable_printing, rng, allocator};
            const results = .{&output_population_mean};

            if (opt.nthread == 1) runTrajectoryParallel(0, T, results, params) else pool.spawnWgId(&wait, runTrajectoryParallel, .{T, results, params});
        }

        wait.wait();
    }

    pool.deinit(); try checkParallelError();

    for (0..opt.nthread) |i| {
        try output.population_mean.add(output_population_mean.at(i));
    }

    output.population_mean.divs(@as(T, @floatFromInt(opt.trajectories)));

    if (enable_printing) for (0..nstate) |i| {

        const population_error = binomialConfInt(output.population_mean.at(opt.iterations, i), opt.trajectories);

        try print("{s}FINAL POPULATION OF STATE {d:2}: {d:.6} ± {:.6}\n", .{if (i == 0) "\n" else "", i, output.population_mean.at(opt.iterations, i), population_error});
    };

    const end_time = @as(T, @floatFromInt(opt.iterations)) * opt.time_step;

    if (opt.write.population_mean) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.population_mean, 0, end_time, opt.iterations + 1);

    return output;
}

/// Run a single trajectory.
pub fn runTrajectory(comptime T: type, opt: Options(T), system: *ClassicalParticle(T), index: usize, enable_printing: bool, allocator: std.mem.Allocator) !Custom(T).TrajectoryOutput {
    const nstate = opt.potential.nstate();

    var output = try Custom(T).TrajectoryOutput.init(nstate, opt.iterations, allocator);

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

    amplitudes.ptr(current_state).* = Complex(T).init(1, 0);

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

    const surface_hopping_parameters: surface_hopping_algorithm.Parameters(T) = .{
        .fs_parameters = fs_parameters,
        .lz_parameters = lz_parameters
    };

    var split_mix = std.Random.SplitMix64.init(opt.seed + index); var rng = std.Random.DefaultPrng.init(split_mix.next()); var random = rng.random();

    var timer = try std.time.Timer.start();

    for (0..opt.iterations + 1) |i| {

        const time = @as(T, @floatFromInt(i)) * opt.time_step;

        try adiabatic_eigenvectors.copyTo(&previous_eigenvectors);

        if (i > 0) {
            try system.propagateVelocityVerlet(opt.potential, &adiabatic_potential, time, current_state, opt.time_step);
        }

        try opt.potential.evaluateEigensystem(&diabatic_potential, &adiabatic_potential, &adiabatic_eigenvectors, system.position, time);

        if (i > 0) {

            try fixGauge(T, &adiabatic_eigenvectors, previous_eigenvectors);
            try mm(T, &eigenvector_overlap, previous_eigenvectors, true, adiabatic_eigenvectors, false);

            if (opt.derivative_coupling) |tdc_algorithm| {
                try tdc_algorithm.evaluate(&time_derivative_coupling, eigenvector_overlap, opt.time_step);
            }
        }

        for (0..nstate) |j| for (j + 1..nstate) |k| {
            energy_gaps.ptr(j * (2 * nstate - j - 1) / 2 + (k - j - 1)).append(adiabatic_potential.at(k, k) - adiabatic_potential.at(j, j));
        };

        if (opt.surface_hopping) |algorithm| if (i > 1) {
            current_state = algorithm.jump(system, &jump_probabilities, surface_hopping_parameters, adiabatic_potential, current_state, &random);
        };

        const kinetic_energy = system.kineticEnergy();
        const potential_energy = adiabatic_potential.at(current_state, current_state);

        output.population.ptr(i, current_state).* += 1;

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

    results[0].ptr(id - 1).add(trajectory_output.population) catch |e| {
        if (PARALLEL_ERROR.* == null) PARALLEL_ERROR.* = e; return;
    };
}

/// Print header for iteration info.
pub fn printIterationHeader(comptime T: type, ndim: usize, nstate: usize, surface_hopping: ?SurfaceHoppingAlgorithm(T)) !void {
    var buffer: [1024]u8 = undefined;

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

/// Prints the iteration info to standard output.
pub fn printIterationInfo(comptime T: type, info: Custom(T).IterationInfo, surface_hopping: ?SurfaceHoppingAlgorithm(T), timer: *std.time.Timer) !void {
    var buffer: [1024]u8 = undefined;

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

    try std.testing.expect(output.population_mean.at(opt.iterations, 0) == 0.371);
    try std.testing.expect(output.population_mean.at(opt.iterations, 1) == 0.629);
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

    try std.testing.expect(output.population_mean.at(opt.iterations, 0) == 0.498);
    try std.testing.expect(output.population_mean.at(opt.iterations, 1) == 0.502);
}
