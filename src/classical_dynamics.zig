//! Code to propagate classical particles.

const std = @import("std");

const andersen_thermostat = @import("andersen_thermostat.zig");
const berendsen_thermostat = @import("berendsen_thermostat.zig");
const bias_potential = @import("bias_potential.zig");
const classical_particle = @import("classical_particle.zig");
const complex_runge_kutta = @import("complex_runge_kutta.zig");
const complex_vector = @import("complex_vector.zig");
const derivative_coupling = @import("derivative_coupling.zig");
const device_write = @import("device_write.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const electronic_potential = @import("electronic_potential.zig");
const error_context = @import("error_context.zig");
const error_handling = @import("error_handling.zig");
const fewest_switches = @import("fewest_switches.zig");
const global_variables = @import("global_variables.zig");
const hammes_schiffer_tully = @import("hammes_schiffer_tully.zig");
const harmonic_potential = @import("harmonic_potential.zig");
const landau_zener = @import("landau_zener.zig");
const langevin_thermostat = @import("langevin_thermostat.zig");
const mapping_approach = @import("mapping_approach.zig");
const math_functions = @import("math_functions.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const nonadiabatic_coupling_vector = @import("nonadiabatic_coupling_vector.zig");
const norm_preserving_interpolation = @import("norm_preserving_interpolation.zig");
const nose_hoover_thermostat = @import("nose_hoover_thermostat.zig");
const object_array = @import("object_array.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");
const ring_buffer = @import("ring_buffer.zig");
const surface_hopping_algorithm = @import("surface_hopping_algorithm.zig");
const thermostat = @import("thermostat.zig");
const trajectory_thermodynamics = @import("trajectory_thermodynamics.zig");
const tully_potential = @import("tully_potential.zig");

const BiasPotential = bias_potential.BiasPotential;
const ClassicalParticle = classical_particle.ClassicalParticle;
const Complex = std.math.complex.Complex;
const ComplexRungeKutta = complex_runge_kutta.ComplexRungeKutta;
const ComplexVector = complex_vector.ComplexVector;
const DerivativeCoupling = derivative_coupling.DerivativeCoupling;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const ErrorContext = error_context.ErrorContext;
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
const Thermostat = thermostat.Thermostat;

const applyDecoherenceCorrection = surface_hopping_algorithm.applyDecoherenceCorrection;
const binomialConfInt = math_functions.binomialConfInt;
const exportRealMatrix = device_write.exportRealMatrix;
const exportRealMatrixWithLinspacedLeftColumn = device_write.exportRealMatrixWithLinspacedLeftColumn;
const fixGauge = eigenproblem_solver.fixGauge;
const h = math_functions.h;
const mm = matrix_multiplication.mm;
const print = device_write.print;
const printJson = device_write.printJson;
const schlitterEntropy = trajectory_thermodynamics.schlitterEntropy;
const sreEntropy = trajectory_thermodynamics.sreEntropy;
const throw = error_handling.throw;

const Eh = global_variables.Eh;
const Na = global_variables.Na;
const MAX_POOL_SIZE = global_variables.MAX_POOL_SIZE;
const TEST_TOLERANCE = global_variables.TEST_TOLERANCE;
const WRITE_BUFFER_SIZE = global_variables.WRITE_BUFFER_SIZE;
const A2AU = global_variables.A2AU;
const AU2K = global_variables.AU2K;

/// Classical dynamics options struct.
pub fn Options(comptime T: type) type {
    return struct {
        pub const InitialConditions = union(enum) {
            model: struct {
                mass: []const T,
                momentum_mean: []const T,
                momentum_std: []const T,
                position_mean: []const T,
                position_std: []const T,
                state: u32 = 0
            },
            molecule : struct {
                position: []const u8,
                velocity: ?[]const u8 = null,
                temperature: ?T = null,
                state: u32 = 0,
                charge: i32 = 0
            }
        };
        pub const LogIntervals = struct {
            trajectory: u32 = 1,
            iteration: u32 = 1
        };
        pub const Thermodynamics = struct {
            schlitter_entropy: bool = false,
            sre_entropy: bool = false
        };
        pub const Write = struct {
            bloch_vector_mean: ?[]const u8 = null,
            coefficient_mean: ?[]const u8 = null,
            kinetic_energy_mean: ?[]const u8 = null,
            momentum_mean: ?[]const u8 = null,
            population_mean: ?[]const u8 = null,
            position_mean: ?[]const u8 = null,
            potential_energy_mean: ?[]const u8 = null,
            state_potential_energy_mean: ?[]const u8 = null,
            time_derivative_coupling_mean: ?[]const u8 = null,
            temperature_mean: ?[]const u8 = null,
            total_energy_mean: ?[]const u8 = null,
            schlitter_entropy: ?[]const u8 = null,
            sre_entropy: ?[]const u8 = null
        };

        potential: ElectronicPotential(T),
        initial_conditions: InitialConditions,

        iterations: u32,
        time_step: T,

        bias: ?BiasPotential(T) = null,
        derivative_coupling: ?DerivativeCoupling(T) = null,
        log_intervals: LogIntervals = .{},
        surface_hopping: ?SurfaceHoppingAlgorithm(T) = null,
        thermostat: ?Thermostat(T) = null,
        thermodynamics: Thermodynamics = .{},
        write: Write = .{},

        equilibration_iterations: u32 = 0,
        finite_differences_step: T = 1e-6,
        trajectories: u32 = 1,
        seed: u32 = 0,
        nthread: u32 = 1
    };
}

/// Structure that hold the output of the simulation.
pub fn Output(comptime T: type) type {
    return struct {
        bloch_vector_mean: RealMatrix(T),
        coefficient_mean: RealMatrix(T),
        kinetic_energy_mean: RealVector(T),
        momentum_mean: RealMatrix(T),
        population_mean: RealMatrix(T),
        position_mean: RealMatrix(T),
        potential_energy_mean: RealVector(T),
        state_potential_energy_mean: RealMatrix(T),
        time_derivative_coupling_mean: RealMatrix(T),
        temperature_mean: RealVector(T),
        total_energy_mean: RealVector(T),
        schlitter_entropy: RealVector(T),
        sre_entropy: RealVector(T),

        /// Allocate the output structure.
        pub fn init(nstate: usize, ndim: usize, iterations: usize, trajectories: usize, allocator: std.mem.Allocator) !@This() {
            return @This(){
                .bloch_vector_mean = try RealMatrix(T).initZero(iterations + 1, 4, allocator),
                .coefficient_mean = try RealMatrix(T).initZero(iterations + 1, nstate, allocator),
                .kinetic_energy_mean = try RealVector(T).initZero(iterations + 1, allocator),
                .momentum_mean = try RealMatrix(T).initZero(iterations + 1, ndim, allocator),
                .population_mean = try RealMatrix(T).initZero(iterations + 1, nstate, allocator),
                .position_mean = try RealMatrix(T).initZero(iterations + 1, ndim, allocator),
                .potential_energy_mean = try RealVector(T).initZero(iterations + 1, allocator),
                .state_potential_energy_mean = try RealMatrix(T).initZero(iterations + 1, nstate, allocator),
                .time_derivative_coupling_mean = try RealMatrix(T).initZero(iterations + 1, nstate * nstate, allocator),
                .temperature_mean = try RealVector(T).initZero(iterations + 1, allocator),
                .total_energy_mean = try RealVector(T).initZero(iterations + 1, allocator),
                .schlitter_entropy = try RealVector(T).initZero(trajectories, allocator),
                .sre_entropy = try RealVector(T).initZero(trajectories, allocator)
            };
        }

        /// Free the memory allocated for the output structure.
        pub fn deinit(self: @This(), allocator: std.mem.Allocator) void {
            self.bloch_vector_mean.deinit(allocator);
            self.coefficient_mean.deinit(allocator);
            self.kinetic_energy_mean.deinit(allocator);
            self.momentum_mean.deinit(allocator);
            self.population_mean.deinit(allocator);
            self.position_mean.deinit(allocator);
            self.potential_energy_mean.deinit(allocator);
            self.state_potential_energy_mean.deinit(allocator);
            self.time_derivative_coupling_mean.deinit(allocator);
            self.temperature_mean.deinit(allocator);
            self.total_energy_mean.deinit(allocator);
            self.schlitter_entropy.deinit(allocator);
            self.sre_entropy.deinit(allocator);
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
            temperature: T,
            coefficient: ComplexVector(T)
        };

        /// Structure to hold a single trajectory output.
        pub const TrajectoryOutput = struct {
            pub const Thermodynamics = struct {
                schlitter_entropy: ?T = null,
                sre_entropy: ?T = null
            };

            bloch_vector: RealMatrix(T),
            coefficient: RealMatrix(T),
            kinetic_energy: RealVector(T),
            momentum: RealMatrix(T),
            population: RealMatrix(T),
            position: RealMatrix(T),
            potential_energy: RealVector(T),
            state_potential_energy: RealMatrix(T),
            time_derivative_coupling: RealMatrix(T),
            temperature: RealVector(T),
            total_energy: RealVector(T),
            thermodynamics: Thermodynamics = .{},

            /// Allocate the trajectory output structure.
            pub fn init(nstate: usize, ndim: usize, iterations: usize, allocator: std.mem.Allocator) !@This() {
                return @This(){
                    .bloch_vector = try RealMatrix(T).initZero(iterations + 1, 4, allocator),
                    .coefficient = try RealMatrix(T).initZero(iterations + 1, nstate, allocator),
                    .kinetic_energy = try RealVector(T).initZero(iterations + 1, allocator),
                    .momentum = try RealMatrix(T).initZero(iterations + 1, ndim, allocator),
                    .population = try RealMatrix(T).initZero(iterations + 1, nstate, allocator),
                    .position = try RealMatrix(T).initZero(iterations + 1, ndim, allocator),
                    .potential_energy = try RealVector(T).initZero(iterations + 1, allocator),
                    .state_potential_energy = try RealMatrix(T).initZero(iterations + 1, nstate, allocator),
                    .time_derivative_coupling = try RealMatrix(T).initZero(iterations + 1, nstate * nstate, allocator),
                    .temperature = try RealVector(T).initZero(iterations + 1, allocator),
                    .total_energy = try RealVector(T).initZero(iterations + 1, allocator),
                };
            }

            /// Free the memory allocated for the trajectory output structure.
            pub fn deinit(self: @This(), allocator: std.mem.Allocator) void {
                self.bloch_vector.deinit(allocator);
                self.coefficient.deinit(allocator);
                self.kinetic_energy.deinit(allocator);
                self.momentum.deinit(allocator);
                self.population.deinit(allocator);
                self.position.deinit(allocator);
                self.potential_energy.deinit(allocator);
                self.state_potential_energy.deinit(allocator);
                self.time_derivative_coupling.deinit(allocator);
                self.temperature.deinit(allocator);
                self.total_energy.deinit(allocator);
            }
        };
    };
}

/// Run classical dynamics simulation.
pub fn run(comptime T: type, opt: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(opt);

    const ndim = opt.potential.ndim();
    const nstate = opt.potential.nstate();

    if (nstate != 2 and opt.write.bloch_vector_mean != null) return throw(Output(T), "BLOCH VECTOR OUTPUT IS ONLY SUPPORTED FOR TWO-STATE SYSTEMS", .{});
    if (opt.potential == .ab_initio and opt.derivative_coupling != null and opt.derivative_coupling.? != .nacv) {
        return throw(Output(T), "ONLY NACV DERIVATIVE COUPLING IS SUPPORTED FOR AB INITIO POTENTIALS", .{});
    }
    if (opt.initial_conditions == .molecule and opt.initial_conditions.molecule.velocity != null and opt.initial_conditions.molecule.temperature != null) {
        return throw(Output(T), "INITIAL VELOCITY AND TEMPERATURE CONDITIONS ARE MUTUALLY EXCLUSIVE", .{});
    }

    var custom_potential = if (opt.potential == .custom) try opt.potential.custom.init(allocator) else null; defer if (custom_potential) |*cp| cp.deinit(allocator);
    var file_potential = if (opt.potential == .file) try opt.potential.file.init(allocator) else null; defer if (file_potential) |*fp| fp.deinit(allocator);

    var output = try Output(T).init(nstate, ndim, opt.iterations, opt.trajectories, allocator); errdefer output.deinit(allocator);

    const parallel_results = .{
        .bloch_vector_mean = try RealMatrixArray(T).initZero(opt.nthread, .{.rows = opt.iterations + 1, .cols = 4}, allocator),
        .coefficient_mean = try RealMatrixArray(T).initZero(opt.nthread, .{.rows = opt.iterations + 1, .cols = nstate}, allocator),
        .kinetic_energy_mean = try RealVectorArray(T).initZero(opt.nthread, .{.rows = opt.iterations + 1}, allocator),
        .momentum_mean = try RealMatrixArray(T).initZero(opt.nthread, .{.rows = opt.iterations + 1, .cols = ndim}, allocator),
        .population_mean = try RealMatrixArray(T).initZero(opt.nthread, .{.rows = opt.iterations + 1, .cols = nstate}, allocator),
        .position_mean = try RealMatrixArray(T).initZero(opt.nthread, .{.rows = opt.iterations + 1, .cols = ndim}, allocator),
        .potential_energy_mean = try RealVectorArray(T).initZero(opt.nthread, .{.rows = opt.iterations + 1}, allocator),
        .state_potential_energy_mean = try RealMatrixArray(T).initZero(opt.nthread, .{.rows = opt.iterations + 1, .cols = nstate}, allocator),
        .time_derivative_coupling_mean = try RealMatrixArray(T).initZero(opt.nthread, .{.rows = opt.iterations + 1, .cols = nstate * nstate}, allocator),
        .temperature_mean = try RealVectorArray(T).initZero(opt.nthread, .{.rows = opt.iterations + 1}, allocator),
        .total_energy_mean = try RealVectorArray(T).initZero(opt.nthread, .{.rows = opt.iterations + 1}, allocator),
    };

    var trajectory_based_results = .{
        .schlitter_entropy = try RealVector(T).initZero(opt.trajectories, allocator),
        .sre_entropy = try RealVector(T).initZero(opt.trajectories, allocator)
    };

    defer inline for (std.meta.fields(@TypeOf(parallel_results))) |field| @as(field.type, @field(parallel_results, field.name)).deinit(allocator);

    defer inline for (std.meta.fields(@TypeOf(trajectory_based_results))) |field| @as(field.type, @field(trajectory_based_results, field.name)).deinit(allocator);

    var split_mix = std.Random.SplitMix64.init(opt.seed);

    if (enable_printing) try printIterationHeader(T, ndim, nstate, opt.surface_hopping, opt.thermostat);

    var pool: std.Thread.Pool = undefined; var wait: std.Thread.WaitGroup = undefined;

    try pool.init(.{.n_jobs = opt.nthread, .track_ids = true, .allocator = allocator});

    var error_ctx = ErrorContext{};

    for (0..(opt.trajectories + MAX_POOL_SIZE - 1) / MAX_POOL_SIZE) |i| {

        wait.reset();

        for (0..@min(opt.trajectories - i * MAX_POOL_SIZE, MAX_POOL_SIZE)) |j| {

            const rng = std.Random.DefaultPrng.init(split_mix.next());

            const params = .{opt, i * MAX_POOL_SIZE + j, enable_printing, rng, allocator};

            if (opt.nthread == 1) {
                runTrajectoryParallel(1, T, parallel_results, &trajectory_based_results, params, &error_ctx);
            } else {
                pool.spawnWgId(&wait, runTrajectoryParallel, .{T, parallel_results, &trajectory_based_results, params, &error_ctx});
            }
        }

        wait.wait();
    }

    pool.deinit(); if (error_ctx.err) |err| return err;

    try finalizeOutput(T, &output, opt, parallel_results, trajectory_based_results);

    if (enable_printing) try printFinalDetails(T, opt, output);

    return output;
}

/// Run a single trajectory.
pub fn runTrajectory(comptime T: type, opt: Options(T), system: *ClassicalParticle(T), index: usize, equilibrate: bool, enable_printing: bool, allocator: std.mem.Allocator) !Custom(T).TrajectoryOutput {
    const nstate = opt.potential.nstate();
    const ndim = opt.potential.ndim();
    var dir = std.fs.cwd();

    if (opt.potential == .ab_initio) {

        const dir_name = try std.fmt.allocPrint(allocator, "TRAJ_{d}", .{index + 1}); defer allocator.free(dir_name);

        try std.fs.cwd().deleteTree(dir_name); try std.fs.cwd().makeDir(dir_name);

        dir = try std.fs.cwd().openDir(dir_name, .{});
    }

    var split_mix = std.Random.SplitMix64.init(opt.seed + index); var rng = std.Random.DefaultPrng.init(split_mix.next()); var random = rng.random();

    var output = try Custom(T).TrajectoryOutput.init(nstate, ndim, opt.iterations, allocator);

    var current_state: usize = switch (opt.initial_conditions) {
        .model => opt.initial_conditions.model.state,
        .molecule => opt.initial_conditions.molecule.state
    };

    if (current_state >= nstate) return throw(Custom(T).TrajectoryOutput, "ACTIVE STATE MUST NOT BE HIGHER THAN THE TOTAL NUMBER OF STATES", .{});

    var diabatic_potential = try RealMatrix(T).init(nstate, nstate, allocator); defer diabatic_potential.deinit(allocator);
    var adiabatic_potential = try RealMatrix(T).init(nstate, nstate, allocator); defer adiabatic_potential.deinit(allocator);
    var adiabatic_eigenvectors = try RealMatrix(T).init(nstate, nstate, allocator); defer adiabatic_eigenvectors.deinit(allocator);
    var previous_eigenvectors = try RealMatrix(T).init(nstate, nstate, allocator); defer previous_eigenvectors.deinit(allocator);
    var eigenvector_overlap = try RealMatrix(T).init(nstate, nstate, allocator); defer eigenvector_overlap.deinit(allocator);
    var time_derivative_coupling = try RealMatrix(T).initZero(nstate, nstate, allocator); defer time_derivative_coupling.deinit(allocator);

    var energy_gaps = try RingBufferArray(T).init(nstate * (nstate - 1) / 2, .{.max_len = 5}, allocator); defer energy_gaps.deinit(allocator);
    var jump_probabilities = try RealVector(T).init(nstate, allocator); defer jump_probabilities.deinit(allocator);
    var runge_kutta_solver = try ComplexRungeKutta(T).init(nstate, allocator); defer runge_kutta_solver.deinit(allocator);
    var coefficient = try ComplexVector(T).initZero(nstate, allocator); defer coefficient.deinit(allocator);
    var bloch_vector = try RealVector(T).initZero(3, allocator); defer bloch_vector.deinit(allocator);

    coefficient.ptr(current_state).* = Complex(T).init(1, 0); if (nstate == 2) bloch_vector.ptr(2).* = coefficient.at(1).squaredMagnitude() - coefficient.at(0).squaredMagnitude();

    var Wcp: T = 1; var Wpp: T = 1; var Sz0: T = bloch_vector.at(2); var xi: T = 0;

    if (opt.surface_hopping != null and opt.surface_hopping.? == .mapping_approach) {

        const phi = 2 * std.math.pi * random.float(T);

        const cos_theta = if (current_state == 1) random.float(T) else random.float(T) - 1;
        const sin_theta = std.math.sqrt(1 - cos_theta * cos_theta);

        bloch_vector.ptr(0).* = sin_theta * std.math.cos(phi);
        bloch_vector.ptr(1).* = sin_theta * std.math.sin(phi);
        bloch_vector.ptr(2).* = cos_theta;

        Sz0 = bloch_vector.at(2); Wcp = 2; Wpp = 2 * @abs(Sz0);
    }

    const fs_parameters: fewest_switches.Parameters(T) = .{
        .adiabatic_potential = adiabatic_potential,
        .coefficient = &coefficient,
        .derivative_coupling = time_derivative_coupling,
        .runge_kutta = runge_kutta_solver,
        .time_step = opt.time_step
    };

    const lz_parameters: landau_zener.Parameters(T) = .{
        .energy_gaps = energy_gaps,
        .time_step = opt.time_step,
        .coefficient = &coefficient
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
        .time = undefined,
        .dir = dir,
        .allocator = allocator
    };

    const npi_parameters: norm_preserving_interpolation.Parameters(T) = .{
        .eigenvector_overlap = eigenvector_overlap,
        .time_step = opt.time_step
    };

    var derivative_coupling_parameters: derivative_coupling.Parameters(T) = .{
        .hst_parameters = hst_parameters,
        .nacv_parameters = nacv_parameters,
        .npi_parameters = npi_parameters
    };

    const andersen_parameters: andersen_thermostat.Parameters(T) = .{
        .time_step = opt.time_step,
        .masses = system.masses,
        .random = &random,
        .molecule = opt.potential == .ab_initio
    };

    const berendsen_parameters: berendsen_thermostat.Parameters(T) = .{
        .time_step = opt.time_step,
        .temperature = undefined
    };

    const langevin_parameters: langevin_thermostat.Parameters(T) = .{
        .time_step = opt.time_step,
        .masses = system.masses,
        .random = &random
    };

    const nose_hoover_parameters: nose_hoover_thermostat.Parameters(T) = .{
        .time_step = opt.time_step,
        .ndof = @as(T, @floatFromInt(system.ndof)),
        .masses = system.masses,
        .xi = &xi
    };

    var thermostat_parameters: thermostat.Parameters(T) = .{
        .andersen_params = andersen_parameters,
        .berendsen_params = berendsen_parameters,
        .langevin_params = langevin_parameters,
        .nose_hoover_params = nose_hoover_parameters
    };

    var timer = try std.time.Timer.start();

    for (0..if (equilibrate) opt.equilibration_iterations + 1 else opt.iterations + 1) |i| {

        const time = @as(T, @floatFromInt(i)) * opt.time_step; derivative_coupling_parameters.nacv_parameters.time = time;

        try adiabatic_eigenvectors.copyTo(&previous_eigenvectors);

        if (i > 0 and opt.thermostat != null) {
            try opt.thermostat.?.apply(&system.velocity, thermostat_parameters, .Before);
        }

        if (i > 0) system.propagateVelocityVerletFirstHalf(opt.time_step);

        if (opt.potential == .ab_initio) try opt.potential.ab_initio.runElectronicStructureCalculation(system.*, dir, allocator);

        try opt.potential.evaluateEigensystem(&diabatic_potential, &adiabatic_potential, &adiabatic_eigenvectors, system.position, time, dir);

        var potential_energy = adiabatic_potential.at(current_state, current_state);

        for (0..nstate) |j| output.state_potential_energy.ptr(i, j).* = Wpp * adiabatic_potential.at(j, j);

        if (!equilibrate and i > 0) {

            if (opt.derivative_coupling) |time_derivative_coupling_algorithm| {

                if (time_derivative_coupling_algorithm != .nacv) {
                    try fixGauge(T, &adiabatic_eigenvectors, previous_eigenvectors);
                    try mm(T, &eigenvector_overlap, previous_eigenvectors, true, adiabatic_eigenvectors, false);
                }

                try time_derivative_coupling_algorithm.evaluate(&time_derivative_coupling, derivative_coupling_parameters);
            }
        }

        if (!equilibrate) for (0..nstate) |j| for (j + 1..nstate) |k| {
            energy_gaps.ptr(j * (2 * nstate - j - 1) / 2 + (k - j - 1)).append(adiabatic_potential.at(k, k) - adiabatic_potential.at(j, j));
        };

        if (opt.surface_hopping) |algorithm| if (!equilibrate and i > (if (algorithm == .landau_zener) @as(usize, 1) else @as(usize, 0))) {

            const new_state = try algorithm.jump(system, &jump_probabilities, surface_hopping_parameters, adiabatic_potential, current_state, &random);

            if (new_state != current_state) {
                current_state = new_state; potential_energy = adiabatic_potential.at(current_state, current_state);
            }

            if (algorithm == .fewest_switches and algorithm.fewest_switches.decoh_alpha != null) {
                applyDecoherenceCorrection(T, &coefficient, adiabatic_potential, system.kineticEnergy(), current_state, opt.time_step, algorithm.fewest_switches.decoh_alpha.?);
            }
        };

        try system.calculateAcceleration(opt.potential, &adiabatic_potential, time, current_state, opt.finite_differences_step, opt.bias, dir);

        if (i > 0) try system.propagateVelocityVerletSecondHalf(opt.time_step);

        if (i > 0 and opt.thermostat != null) {

            thermostat_parameters.berendsen_params.temperature = system.kineticTemperature();

            try opt.thermostat.?.apply(&system.velocity, thermostat_parameters, .After);
        }

        const kinetic_energy = system.kineticEnergy(); const temperature = system.kineticTemperature();

        output.kinetic_energy.ptr(i).* = Wpp * kinetic_energy;
        output.population.ptr(i, current_state).* = Wpp;
        output.potential_energy.ptr(i).* = Wpp * potential_energy;
        output.total_energy.ptr(i).* = Wpp * (kinetic_energy + potential_energy);
        output.temperature.ptr(i).* = Wpp * temperature;

        for (0..nstate * nstate) |j| output.time_derivative_coupling.ptr(i, j).* = Wpp * time_derivative_coupling.at(j / nstate, j % nstate);

        for (0..ndim) |j| {
            output.position.ptr(i, j).* = Wpp * system.position.at(j); output.momentum.ptr(i, j).* = Wpp * system.velocity.at(j) * system.masses.at(j);
        }

        if (!equilibrate and opt.surface_hopping != null and (opt.surface_hopping.? == .fewest_switches or opt.surface_hopping.? == .landau_zener)) {

            for (0..nstate) |j| output.coefficient.ptr(i, j).* = coefficient.at(j).squaredMagnitude();

            output.bloch_vector.ptr(i, 0).* = 2 * coefficient.at(0).mul(coefficient.at(1).conjugate()).re;
            output.bloch_vector.ptr(i, 1).* = 2 * coefficient.at(0).mul(coefficient.at(1).conjugate()).im;

            output.bloch_vector.ptr(i, 2).* = coefficient.at(1).squaredMagnitude() - coefficient.at(0).squaredMagnitude();

            output.bloch_vector.ptr(i, 3).* = std.math.sqrt(std.math.pow(T, output.bloch_vector.at(i, 0), 2) + std.math.pow(T, output.bloch_vector.at(i, 1), 2));
        }

        if (!equilibrate and opt.surface_hopping != null and opt.surface_hopping.? == .mapping_approach) {

            for (0..3) |j| output.bloch_vector.ptr(i, j).* = Wcp * bloch_vector.at(j);

            output.bloch_vector.ptr(i, 2).* = Wpp * std.math.sign(bloch_vector.at(2));

            output.coefficient.ptr(i, 0).* = (1 - output.bloch_vector.at(i, 2)) / 2;
            output.coefficient.ptr(i, 1).* = (1 + output.bloch_vector.at(i, 2)) / 2;

            output.bloch_vector.ptr(i, 3).* = std.math.sqrt(std.math.pow(T, output.bloch_vector.at(i, 0), 2) + std.math.pow(T, output.bloch_vector.at(i, 1), 2));
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
            .temperature = temperature,
            .coefficient = coefficient
        };

        try printIterationInfo(T, iteration_info, opt.surface_hopping, opt.thermostat, equilibrate, &timer);
    }

    if (!equilibrate) {

        try getThermodynamicProperties(T, &output, system.*, opt, allocator);

        if (enable_printing and opt.trajectories == 1) try printThermodynamicProperties(T, output.thermodynamics);
    }

    if (opt.potential == .ab_initio) dir.close();

    return output;
}

/// Sample initial conditions for a trajectory from Boltzmann distribution.
pub fn sampleFromBoltzmann(comptime T: type, system: *ClassicalParticle(T), temperature: T, random: *std.Random) void {
    const kBT = temperature / AU2K;

    for (0..system.velocity.len) |i| {
        system.velocity.ptr(i).* = std.math.sqrt(kBT / system.masses.at(i)) * random.floatNorm(T);
    }

    return;
}

/// Parallel function to run a trajectory.
pub fn runTrajectoryParallel(id: usize, comptime T: type, results: anytype, trajectory_based_results: anytype, params: anytype, error_ctx: *ErrorContext) void {
    if (error_ctx.err != null) return;

    var system: ClassicalParticle(T) = undefined; var rng = params[3]; var random = rng.random();

    if (params[0].potential == .ab_initio) {

        const position = params[0].initial_conditions.molecule.position;
        const velocity = params[0].initial_conditions.molecule.velocity;
        const charge = params[0].initial_conditions.molecule.charge;

        system = classical_particle.read(T, position, charge, params[1], params[4]) catch |err| {
            error_ctx.capture(err); return;
        };

        if (velocity) |vel| {

            const velocity_system = classical_particle.read(T, vel, charge, params[1], params[4]) catch |err| {
                error_ctx.capture(err); return;
            }; defer velocity_system.deinit(params[4]);

            for (0..system.velocity.len) |i| system.velocity.ptr(i).* = velocity_system.position.at(i) / A2AU;
        }

        if (params[0].initial_conditions.molecule.temperature) |temp| sampleFromBoltzmann(T, &system, temp, &random);
    }

    if (params[0].potential != .ab_initio) system = ClassicalParticle(T).initZero(params[0].potential.ndim(), params[0].initial_conditions.model.mass, params[4]) catch |err| {
        error_ctx.capture(err); return;
    };

    if (params[0].potential != .ab_initio) sampleInitialConditions(T, &system, params[0].initial_conditions.model, &random) catch |err| {
        error_ctx.capture(err); return;
    };

    if (params[0].equilibration_iterations > 0) {
        const equilibration_output = runTrajectory(T, params[0], &system, params[1], true, params[2], params[4]) catch |err| {
            error_ctx.capture(err); return;
        }; defer equilibration_output.deinit(params[4]);
    }

    const trajectory_output = runTrajectory(T, params[0], &system, params[1], false, params[2], params[4]) catch |err| {
        error_ctx.capture(err); return;
    }; defer trajectory_output.deinit(params[4]);

    inline for (std.meta.fields(@TypeOf(trajectory_output.thermodynamics))) |field| {
        if (@field(trajectory_output.thermodynamics, field.name)) |val| @field(trajectory_based_results, field.name).ptr(params[1]).* = val;
    }
    
    inline for (std.meta.fields(@TypeOf(results))) |field| {

        const result = @as(field.type, @field(results, field.name)).ptr(id - 1); const output = @field(trajectory_output, field.name[0..field.name.len - 5]);

        result.add(output) catch |err| {
            error_ctx.capture(err); return;
        };
    }

    system.deinit(params[4]);
}

/// Print header for iteration info.
pub fn printIterationHeader(comptime T: type, ndim: usize, nstate: usize, surface_hopping: ?SurfaceHoppingAlgorithm(T), thm: ?Thermostat(T)) !void {
    var buffer: [WRITE_BUFFER_SIZE]u8 = undefined;

    const ndim_header_width = 9 * @as(usize, @min(ndim, 3)) + 2 * (@as(usize, @min(ndim, 3)) - 1) + @as(usize, if (ndim > 3) 7 else 2);
    const nstate_header_width = 7 * @as(usize, @min(4, nstate)) + 2 * (@as(usize, @min(4, nstate)) - 1) + @as(usize, if (nstate > 4) 7 else 2);

    var writer = std.io.Writer.fixed(&buffer);

    try writer.print("\n{s:9} {s:8} ", .{"TRAJ", "ITER"});
    try writer.print("{s:12} {s:12} {s:12} ", .{"KIN (Eh)", "POT (Eh)", "TOT (Eh)"});

    if (thm) |_| try writer.print("{s:10} ", .{"TEMP (K)"});

    try writer.print("{s:5} ", .{"STATE"});
    try writer.print("{[value]s:[width]} ", .{.value = "POS (a0)", .width = ndim_header_width});
    try writer.print("{[value]s:[width]} ", .{.value = "MOM (a0)", .width = ndim_header_width});

    if (surface_hopping) |algorithm| switch (algorithm) {
        .fewest_switches => try writer.print("{[value]s:[width]} ", .{.value = "|COEFS|^2", .width = nstate_header_width}),
        else => {}
    };

    try writer.print("{s:4}", .{"TIME"});

    try print("{s}\n", .{writer.buffered()});
}

/// Add the partial results to the output vectors and export them if requested.
pub fn finalizeOutput(comptime T: type, output: *Output(T), opt: Options(T), parallel_results: anytype, trajectory_based_results: anytype) !void {
    const output_bloch_vector_mean = parallel_results.bloch_vector_mean;
    const output_coefficient_mean = parallel_results.coefficient_mean;
    const output_kinetic_energy_mean = parallel_results.kinetic_energy_mean;
    const output_momentum_mean = parallel_results.momentum_mean;
    const output_population_mean = parallel_results.population_mean;
    const output_position_mean = parallel_results.position_mean;
    const output_potential_energy_mean = parallel_results.potential_energy_mean;
    const output_state_potential_energy_mean = parallel_results.state_potential_energy_mean;
    const output_time_derivative_coupling_mean = parallel_results.time_derivative_coupling_mean;
    const output_temperature_mean = parallel_results.temperature_mean;
    const output_total_energy_mean = parallel_results.total_energy_mean;

    for (0..opt.nthread) |i| {
        try output.bloch_vector_mean.add(output_bloch_vector_mean.at(i));
        try output.coefficient_mean.add(output_coefficient_mean.at(i));
        try output.kinetic_energy_mean.add(output_kinetic_energy_mean.at(i));
        try output.momentum_mean.add(output_momentum_mean.at(i));
        try output.population_mean.add(output_population_mean.at(i));
        try output.position_mean.add(output_position_mean.at(i));
        try output.potential_energy_mean.add(output_potential_energy_mean.at(i));
        try output.state_potential_energy_mean.add(output_state_potential_energy_mean.at(i));
        try output.time_derivative_coupling_mean.add(output_time_derivative_coupling_mean.at(i));
        try output.temperature_mean.add(output_temperature_mean.at(i));
        try output.total_energy_mean.add(output_total_energy_mean.at(i));
    }

    output.bloch_vector_mean.divs(@as(T, @floatFromInt(opt.trajectories)));
    output.coefficient_mean.divs(@as(T, @floatFromInt(opt.trajectories)));
    output.kinetic_energy_mean.divs(@as(T, @floatFromInt(opt.trajectories)));
    output.momentum_mean.divs(@as(T, @floatFromInt(opt.trajectories)));
    output.population_mean.divs(@as(T, @floatFromInt(opt.trajectories)));
    output.position_mean.divs(@as(T, @floatFromInt(opt.trajectories)));
    output.potential_energy_mean.divs(@as(T, @floatFromInt(opt.trajectories)));
    output.state_potential_energy_mean.divs(@as(T, @floatFromInt(opt.trajectories)));
    output.time_derivative_coupling_mean.divs(@as(T, @floatFromInt(opt.trajectories)));
    output.temperature_mean.divs(@as(T, @floatFromInt(opt.trajectories)));
    output.total_energy_mean.divs(@as(T, @floatFromInt(opt.trajectories)));

    const end_time = @as(T, @floatFromInt(opt.iterations)) * opt.time_step;

    if (opt.write.bloch_vector_mean) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.bloch_vector_mean, 0, end_time);
    if (opt.write.coefficient_mean) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.coefficient_mean, 0, end_time);
    if (opt.write.kinetic_energy_mean) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.kinetic_energy_mean.asMatrix(), 0, end_time);
    if (opt.write.momentum_mean) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.momentum_mean, 0, end_time);
    if (opt.write.population_mean) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.population_mean, 0, end_time);
    if (opt.write.position_mean) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.position_mean, 0, end_time);
    if (opt.write.potential_energy_mean) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.potential_energy_mean.asMatrix(), 0, end_time);
    if (opt.write.state_potential_energy_mean) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.state_potential_energy_mean, 0, end_time);
    if (opt.write.time_derivative_coupling_mean) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.time_derivative_coupling_mean, 0, end_time);
    if (opt.write.temperature_mean) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.temperature_mean.asMatrix(), 0, end_time);
    if (opt.write.total_energy_mean) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.total_energy_mean.asMatrix(), 0, end_time);

    if (opt.write.schlitter_entropy) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, trajectory_based_results.schlitter_entropy.asMatrix(), 1, @as(T, @floatFromInt(opt.trajectories)));
    if (opt.write.sre_entropy) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, trajectory_based_results.sre_entropy.asMatrix(), 1, @as(T, @floatFromInt(opt.trajectories)));
}

/// Calculate thermodynamic properties from the trajectory output.
pub fn getThermodynamicProperties(comptime T: type, output: *Custom(T).TrajectoryOutput, system: ClassicalParticle(T), opt: Options(T), allocator: std.mem.Allocator) !void {
    const temp = if (opt.thermostat) |thm| switch (thm) {
        inline else => |field| field.temperature
    } else 0;

    if (opt.thermodynamics.schlitter_entropy or opt.write.schlitter_entropy != null) {
        output.thermodynamics.schlitter_entropy = try schlitterEntropy(T, output.position, system.masses, temp, allocator);
    }

    if (opt.thermodynamics.sre_entropy or opt.write.sre_entropy != null) {
        output.thermodynamics.sre_entropy = try sreEntropy(T, output.momentum, system.masses, temp, opt.time_step, system.ndof, allocator);
    }
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

        const print_payload = .{if (i == 0) "\n" else "", i, output.population_mean.at(opt.iterations, i), if (std.math.isNan(population_error)) 0 else population_error};

        try print("{s}FINAL POPULATION OF STATE {d:2}: {d:.6} ± {:.6}\n", print_payload);
    }
}

/// Prints the iteration info to standard output.
pub fn printIterationInfo(comptime T: type, info: Custom(T).IterationInfo, surface_hopping: ?SurfaceHoppingAlgorithm(T), thm: ?Thermostat(T), equilibrate: bool, timer: *std.time.Timer) !void {
    var buffer: [WRITE_BUFFER_SIZE]u8 = undefined;

    var writer = std.io.Writer.fixed(&buffer);

    if (equilibrate) try writer.print("{d:6}-EQ {d:8} ", .{info.trajectory + 1, info.iteration}) else try writer.print("{d:9} {d:8} ", .{info.trajectory + 1, info.iteration});

    try writer.print("{d:12.6} {d:12.6} {d:12.6} ", .{info.kinetic_energy, info.potential_energy, info.kinetic_energy + info.potential_energy});

    if (thm) |_| {
        if (info.temperature < 1e7) try writer.print("{d:10.3} ", .{info.temperature}) else try writer.print("{e:10.4} ", .{info.temperature});
    }

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

            for (0..@min(4, info.coefficient.len)) |i| {
                try writer.print("{d:7.4}{s}", .{std.math.pow(T, info.coefficient.at(i).magnitude(), 2), if (i == info.coefficient.len - 1) "" else ", "});
            }

            if (info.coefficient.len > 4) try writer.print("...", .{});
        }
    }

    try writer.print("] {D}", .{timer.read()}); timer.reset();

    try print("{s}\n", .{writer.buffered()});
}

/// Print the thermodynamic properties to standard output.
pub fn printThermodynamicProperties(comptime T: type, output: Custom(T).TrajectoryOutput.Thermodynamics) !void {
    inline for (std.meta.fields(@TypeOf(output))) |field| if (@field(output, field.name) != null) {
        try print("\n", .{}); break;
    };

    if (output.schlitter_entropy) |out| try print("SCHLITTER ENTROPY: {d:.8} Eh/K = {d:.8} J/MOL/K\n", .{out, out * Na * Eh});
    if (output.sre_entropy) |out| try print("SRE ENTROPY: {d:.8} Eh/K = {d:.8} J/MOL/K\n", .{out, out * Na * Eh});
}

/// Samples the initial conditions.
pub fn sampleInitialConditions(comptime T: type, system: *ClassicalParticle(T), initial_conditions: anytype, random: *std.Random) !void {
    try system.setPositionRandn(initial_conditions.position_mean, initial_conditions.position_std, random);
    try system.setMomentumRandn(initial_conditions.momentum_mean, initial_conditions.momentum_std, random);
}

test "Fewest Switches Surface Hopping on Tully's First Potential" {
    const opt = Options(f64){
        .derivative_coupling = .{
            .npi = NormPreservingInterpolation(f64){}
        },
        .initial_conditions = .{
            .model = .{
                .mass = &.{2000},
                .momentum_mean = &.{15},
                .momentum_std = &.{1},
                .position_mean = &.{-10},
                .position_std = &.{0.5},
                .state = 1
            }
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

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit(std.testing.allocator);

    try std.testing.expectEqual(output.population_mean.at(opt.iterations, 0), 0.372);
    try std.testing.expectEqual(output.population_mean.at(opt.iterations, 1), 0.628);
}

test "Landau-Lener Surface Hopping on Tully's First Potential" {
    const opt = Options(f64){
        .initial_conditions = .{
            .model = .{
                .mass = &.{2000},
                .momentum_mean = &.{15},
                .momentum_std = &.{1},
                .position_mean = &.{-10},
                .position_std = &.{0.5},
                .state = 1
            }
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

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit(std.testing.allocator);

    try std.testing.expectEqual(output.population_mean.at(opt.iterations, 0), 0.498);
    try std.testing.expectEqual(output.population_mean.at(opt.iterations, 1), 0.502);
}

test "Mapping Approach to Surface Hopping on Tully's First Potential" {
    const opt = Options(f64){
        .derivative_coupling = .{
            .npi = NormPreservingInterpolation(f64){}
        },
        .initial_conditions = .{
            .model = .{
                .mass = &.{2000},
                .momentum_mean = &.{15},
                .momentum_std = &.{1},
                .position_mean = &.{-10},
                .position_std = &.{0.5},
                .state = 1
            }
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

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(output.population_mean.at(opt.iterations, 0), 0.40697437617713, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.population_mean.at(opt.iterations, 1), 0.58171004312212, TEST_TOLERANCE);
}
