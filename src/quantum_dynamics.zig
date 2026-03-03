//! Target for performing quantum dynamics simulations.

const std = @import("std");

const array_functions = @import("array_functions.zig");
const complex_absorbing_potential = @import("complex_absorbing_potential.zig");
const complex_matrix = @import("complex_matrix.zig");
const complex_vector = @import("complex_vector.zig");
const device_write = @import("device_write.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const electronic_potential = @import("electronic_potential.zig");
const error_handling = @import("error_handling.zig");
const fourier_transform = @import("fourier_transform.zig");
const global_variables = @import("global_variables.zig");
const grid_generator = @import("grid_generator.zig");
const grid_wavefunction = @import("grid_wavefunction.zig");
const harmonic_potential = @import("harmonic_potential.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");
const tully_potential = @import("tully_potential.zig");

const ComplexAbsorbingPotential = complex_absorbing_potential.ComplexAbsorbingPotential;
const Complex = std.math.Complex;
const ComplexMatrix = complex_matrix.ComplexMatrix;
const ComplexVector = complex_vector.ComplexVector;
const ElectronicPotential = electronic_potential.ElectronicPotential;
const GridWavefunction = grid_wavefunction.GridWavefunction;
const HarmonicPotential = harmonic_potential.HarmonicPotential;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;
const TullyPotential1 = tully_potential.TullyPotential1;

const cfft1 = fourier_transform.cfft1;
const exportComplexMatrixWithLinspacedLeftColumn = device_write.exportComplexMatrixWithLinspacedLeftColumn;
const exportRealMatrix = device_write.exportRealMatrix;
const exportRealMatrixWithLinspacedLeftColumn = device_write.exportRealMatrixWithLinspacedLeftColumn;
const fixGauge = eigenproblem_solver.fixGauge;
const mm = matrix_multiplication.mm;
const momentumGridAlloc = grid_generator.momentumGridAlloc;
const positionAtRow = grid_generator.positionAtRow;
const positionGridAlloc = grid_generator.positionGridAlloc;
const print = device_write.print;
const printJson = device_write.printJson;
const prod = array_functions.prod;
const throw = error_handling.throw;

const TEST_TOLERANCE = global_variables.TEST_TOLERANCE;
const WRITE_BUFFER_SIZE = global_variables.WRITE_BUFFER_SIZE;
const MAX_PATH_LENGTH = global_variables.MAX_PATH_LENGTH;
const QD_MOMENTUM_TOLERANCE = global_variables.QD_MOMENTUM_TOLERANCE;

/// The quantum dynamics options struct.
pub fn Options(comptime T: type) type {
    return struct {
        pub const Grid = struct {
            limits: []const []const T,
            points: u32
        };
        pub const Imaginary = struct {
            states: u32 = 1
        };
        pub const InitialConditions = struct {
            pub const Spread = struct {
                momentum: ?struct {
                    end: []const T,
                    step: []const T
                } = null,
                position: ?struct {
                    end: []const T,
                    step: []const T
                } = null
            };

            adiabatic: bool = false,
            gamma: []const T,
            mass: T,
            momentum: []const T,
            position: []const T,
            state: u32,
            spread: ?Spread = null
        };
        pub const LogIntervals = struct {
            iteration: u32 = 1
        };
        pub const Spectrum = struct {
            window: enum {gaussian} = .gaussian,
            padding_order: u32 = 1
        };
        pub const Write = struct {
            spatial_bloch_vector: ?[]const u8 = null,
            bloch_vector: ?[]const u8 = null,
            final_wavefunction: ?[]const u8 = null,
            kinetic_energy: ?[]const u8 = null,
            momentum: ?[]const u8 = null,
            population: ?[]const u8 = null,
            position: ?[]const u8 = null,
            potential_energy: ?[]const u8 = null,
            wavefunction: ?[]const u8 = null,
            spectrum: ?[]const u8 = null,
            total_energy: ?[]const u8 = null,
            autocorrelation_function: ?[]const u8 = null,
            transition_probability: ?[]const u8 = null
        };

        potential: ElectronicPotential(T),
        grid: Grid,
        initial_conditions: InitialConditions,

        iterations: u32,
        time_step: T,

        cap: ?ComplexAbsorbingPotential(T) = null,
        imaginary: ?Imaginary = null,
        log_intervals: LogIntervals = .{},
        spectrum: Spectrum = .{},
        write: Write = .{},

        adiabatic: bool = false,
        fix_gauge: bool = true
    };
}

/// The quantum dynamics output struct.
pub fn Output(comptime T: type) type {
    return struct {
        bloch_vector: RealMatrix(T),
        kinetic_energy: RealVector(T),
        momentum: RealMatrix(T),
        population: RealMatrix(T),
        position: RealMatrix(T),
        potential_energy: RealVector(T),
        total_energy: RealVector(T),

        /// Allocate the output structure.
        pub fn init(nstate: usize, ndim: usize, iterations: usize, allocator: std.mem.Allocator) !@This() {
            return @This(){
                .bloch_vector = try RealMatrix(T).initZero(iterations + 1, 4, allocator),
                .kinetic_energy = try RealVector(T).initZero(iterations + 1, allocator),
                .momentum = try RealMatrix(T).initZero(iterations + 1, ndim, allocator),
                .population = try RealMatrix(T).initZero(iterations + 1, nstate, allocator),
                .position = try RealMatrix(T).initZero(iterations + 1, ndim, allocator),
                .potential_energy = try RealVector(T).initZero(iterations + 1, allocator),
                .total_energy = try RealVector(T).initZero(iterations + 1, allocator)
            };
        }

        /// Deallocate the output structure.
        pub fn deinit(self: @This(), allocator: std.mem.Allocator) void {
            self.bloch_vector.deinit(allocator);
            self.kinetic_energy.deinit(allocator);
            self.momentum.deinit(allocator);
            self.population.deinit(allocator);
            self.position.deinit(allocator);
            self.potential_energy.deinit(allocator);
            self.total_energy.deinit(allocator);
        }
    };
}

/// Container for custom structs related to the quantum dynamics.
pub fn Custom(comptime T: type) type {
    return struct {

        /// Structure to hold information about each iteration used for logging.
        pub const IterationInfo = struct {
            density_matrix: ComplexMatrix(T),
            population_loss: RealVector(T),
            iteration: usize,
            kinetic_energy: T,
            momentum: RealVector(T),
            position: RealVector(T),
            potential_energy: T,
            time: T
        };
    };
}

/// Run quantum dynamics simulations.
pub fn run(comptime T: type, opt: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(opt);

    if (opt.potential == .ab_initio) return throw(Output(T), "AB INITIO POTENTIAL IS NOT SUPPORTED IN QUANTUM DYNAMICS SIMULATIONS", .{});

    const nstate = opt.potential.nstate();
    const ndim = try opt.potential.ndim();

    if (nstate != 2 and (opt.write.bloch_vector != null or opt.write.spatial_bloch_vector != null)) return throw(Output(T), "BLOCH VECTOR OUTPUT IS ONLY SUPPORTED FOR TWO-STATE SYSTEMS", .{});

    var custom_potential = if (opt.potential == .custom) try opt.potential.custom.init(allocator) else null; defer if (custom_potential) |*cp| cp.deinit(allocator);
    var file_potential = if (opt.potential == .file) try opt.potential.file.init(allocator) else null; defer if (file_potential) |*fp| fp.deinit(allocator);

    if (opt.initial_conditions.spread) |ics| if (ics.position) |qstruct| if (qstruct.end.len != ndim or qstruct.step.len != ndim) {
        return throw(Output(T), "INITIAL CONDITIONS POSITION SPREAD MUST HAVE THE SAME LENGTH AS THE NUMBER OF DIMENSIONS", .{});
    };

    if (opt.initial_conditions.spread) |ics| if (ics.momentum) |pstruct| if (pstruct.end.len != ndim or pstruct.step.len != ndim) {
        return throw(Output(T), "INITIAL CONDITIONS MOMENTUM SPREAD MUST HAVE THE SAME LENGTH AS THE NUMBER OF DIMENSIONS", .{});
    };

    var output: Output(T) = undefined;

    const ics_pos = if (opt.initial_conditions.spread) |ics| ics.position else null;
    const ics_mom = if (opt.initial_conditions.spread) |ics| ics.momentum else null;

    const q0 = opt.initial_conditions.position; var q1: []const T = undefined;
    const p0 = opt.initial_conditions.momentum; var p1: []const T = undefined;

    if (ics_pos) |qstruct| q1 = qstruct.end else q1 = q0;
    if (ics_mom) |pstruct| p1 = pstruct.end else p1 = p0;

    var qsteps_i = try allocator.alloc(usize, ndim); defer allocator.free(qsteps_i); var qsteps: usize = 1; @memset(qsteps_i, 1);
    var psteps_i = try allocator.alloc(usize, ndim); defer allocator.free(psteps_i); var psteps: usize = 1; @memset(psteps_i, 1);

    if (ics_pos) |qstruct| for (0..qstruct.step.len) |i| {

        if (qstruct.step[i] == 0) return throw(Output(T), "INITIAL CONDITIONS POSITION SPREAD STEP CANNOT BE ZERO", .{});

        if (qstruct.end[i] < q0[i] and qstruct.step[i] > 0) return throw(Output(T), "INITIAL CONDITIONS POSITION SPREAD END MUST BE GREATER THAN START FOR POSITIVE STEP", .{});

        qsteps_i[i] = @as(usize, @intFromFloat(std.math.ceil((q1[i] - q0[i]) / qstruct.step[i]))) + 1;

        qsteps *= qsteps_i[i];
    };

    if (ics_mom) |pstruct| for (0..pstruct.step.len) |i| {

        if (pstruct.step[i] == 0) return throw(Output(T), "INITIAL CONDITIONS MOMENTUM SPREAD STEP CANNOT BE ZERO", .{});

        if (pstruct.end[i] < p0[i] and pstruct.step[i] > 0) return throw(Output(T), "INITIAL CONDITIONS MOMENTUM SPREAD END MUST BE GREATER THAN START FOR POSITIVE STEP", .{});

        psteps_i[i] = @as(usize, @intFromFloat(std.math.ceil((p1[i] - p0[i]) / pstruct.step[i]))) + 1;

        psteps *= psteps_i[i];
    };

    var q = try allocator.alloc(T, ndim); defer allocator.free(q);
    var p = try allocator.alloc(T, ndim); defer allocator.free(p);

    var transition_probability = try RealMatrix(T).initZero(qsteps * psteps, 2 * ndim + nstate, allocator); defer transition_probability.deinit(allocator);

    for (0..qsteps) |i| for (0..psteps) |j| {

        var temp_i = i; var temp_j = j;

        for (0..ndim) |k| {q[k] = q0[k] + @as(T, @floatFromInt(temp_i % qsteps_i[k])) * if (ics_pos) |qstruct| qstruct.step[k] else 0; temp_i /= qsteps_i[k];}
        for (0..ndim) |k| {p[k] = p0[k] + @as(T, @floatFromInt(temp_j % psteps_i[k])) * if (ics_mom) |pstruct| pstruct.step[k] else 0; temp_j /= psteps_i[k];}

        var options = opt; options.initial_conditions.position = q; options.initial_conditions.momentum = p;

        try renameOutputFilesWithPositionAndMomentum(T, &options, options.initial_conditions.position, options.initial_conditions.momentum);

        const result = try performDynamics(T, options, enable_printing, allocator);

        for (0..ndim) |k| transition_probability.ptr(i * psteps + j,        k).* = q[k];
        for (0..ndim) |k| transition_probability.ptr(i * psteps + j, ndim + k).* = p[k];

        for (0..nstate) |k| transition_probability.ptr(i * psteps + j, 2 * ndim + k).* = result.population.at(result.population.rows - 1, k);

        if (i == 0 and j == 0) output = result else result.deinit(allocator);
    };

    if (opt.write.transition_probability) |path| try exportRealMatrix(T, path, transition_probability);

    return output;
}

/// Appends the position and momentum to all output files.
pub fn renameOutputFilesWithPositionAndMomentum(comptime T: type, opt: *Options(T), q: []const T, p: []const T) !void {
    inline for (std.meta.fields(@TypeOf(opt.write))) |field| if (@as(field.type, @field(opt.write, field.name)) != null) {

        const path = &@field(opt.write, field.name).?; var new_path: [MAX_PATH_LENGTH]u8 = undefined;

        var fbs = std.io.fixedBufferStream(&new_path); const writer = fbs.writer();

        try writer.print("{s}_Q=", .{path.*[0 .. path.len - 4]});

        for (q, 0..) |val, i| {
            if (i > 0) try writer.writeAll(","); try writer.print("{d:.4}", .{val});
        }

        try writer.writeAll("_P=");

        for (p, 0..) |val, i| {
            if (i > 0) try writer.writeAll(","); try writer.print("{d:.4}", .{val});
        }

        try writer.writeAll(".mat"); path.* = fbs.getWritten();
    };
}

/// Perform the quantum dynamics simulation.
pub fn performDynamics(comptime T: type, opt: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    const ndim = try opt.potential.ndim();
    const nstate = opt.potential.nstate();

    var output = try Output(T).init(nstate, ndim, @intCast(opt.iterations), allocator); errdefer output.deinit(allocator);

    var wavefunction = try GridWavefunction(T).init(@intCast(opt.grid.points), nstate, ndim, opt.grid.limits, opt.initial_conditions.mass, allocator); defer wavefunction.deinit(allocator);

    var wavefunction_dynamics: ?RealMatrix(T) = if (opt.write.wavefunction) |_| try initializeTemporalGridContainer(T, wavefunction, opt.iterations, 2 * wavefunction.nstate, allocator) else null;
    var spatial_bloch: ?RealMatrix(T) = if (opt.write.spatial_bloch_vector) |_| try initializeTemporalGridContainer(T, wavefunction, opt.iterations, 3, allocator) else null;

    var optimized_wavefunctions = std.ArrayList(GridWavefunction(T)){}; defer optimized_wavefunctions.deinit(allocator);

    var density_matrix = try ComplexMatrix(T).init(nstate, nstate, allocator); defer density_matrix.deinit(allocator);
    var position = try RealVector(T).init(ndim, allocator); defer position.deinit(allocator);
    var momentum = try RealVector(T).init(ndim, allocator); defer momentum.deinit(allocator);
    var acf = try ComplexVector(T).init(opt.iterations + 1, allocator); defer acf.deinit(allocator);

    var timer = try std.time.Timer.start(); const n_dynamics = if (opt.imaginary) |f| f.states else 1;

    for (0..n_dynamics) |i| {

        try wavefunction.initialGaussian(opt.initial_conditions.position, opt.initial_conditions.momentum, opt.initial_conditions.state, opt.initial_conditions.gamma);

        if (opt.initial_conditions.adiabatic) try wavefunction.transformRepresentation(opt.potential, 0, false);

        const save_iwf = (opt.write.spectrum != null or opt.write.autocorrelation_function != null) and i == n_dynamics - 1;

        const initial_wavefunction: ?GridWavefunction(T) = if (save_iwf) try wavefunction.clone(allocator) else null; defer if (initial_wavefunction) |iwf| iwf.deinit(allocator);

        if (opt.imaginary != null) for (optimized_wavefunctions.items) |owf| wavefunction.orthogonalize(owf);

        if (enable_printing) {

            if (opt.imaginary != null) try print("\nIMAGINARY TIME PROPAGATION - STATE {d}/{d}", .{i + 1, n_dynamics});

            try printIterationHeader(ndim, nstate);
        }

        for (0..opt.iterations + 1) |j| {

            const time = @as(T, @floatFromInt(j)) * opt.time_step;

            if (j > 0) try wavefunction.propagate(opt.potential, opt.cap, time, opt.time_step, opt.adiabatic, opt.imaginary != null, opt.fix_gauge);

            if (j > 0 and opt.imaginary != null) for (optimized_wavefunctions.items) |owf| wavefunction.orthogonalize(owf);

            if (initial_wavefunction) |iwf| acf.ptr(j).* = iwf.overlap(wavefunction);

            try wavefunction.density(&density_matrix, opt.potential, time, opt.adiabatic, opt.fix_gauge);
            wavefunction.positionMean(&position);

            const kinetic_energy = try wavefunction.kineticEnergyAndMomentumMean(&momentum);
            const potential_energy = try wavefunction.potentialEnergy(opt.potential, time);

            if (j == 0) for (0..momentum.len) |k| if (@abs(momentum.at(k) - opt.initial_conditions.momentum[k]) > QD_MOMENTUM_TOLERANCE) {
                return throw(Output(T), "YOUR GRID IS NOT SUFFICIENTLY DENSE TO ACCURATELY REPRESENT THE MOMENTUM, CONSIDER INCREASING THE NUMBER OF GRID POINTS OR TIGHTENING THE GRID LIMITS", .{});
            };

            if (i == n_dynamics - 1) {
                output.kinetic_energy.ptr(j).* = kinetic_energy;
                output.potential_energy.ptr(j).* = potential_energy;
                output.total_energy.ptr(j).* = kinetic_energy + potential_energy;
                for (0..ndim) |k| output.momentum.ptr(j, k).* = momentum.at(k);
                for (0..ndim) |k| output.position.ptr(j, k).* = position.at(k);
                for (0..nstate) |k| output.population.ptr(j, k).* = density_matrix.at(k, k).re + wavefunction.population_loss.at(k);
                if (nstate == 2) output.bloch_vector.ptr(j, 0).* = 2 * density_matrix.at(0, 1).re;
                if (nstate == 2) output.bloch_vector.ptr(j, 1).* = 2 * density_matrix.at(0, 1).im;
                if (nstate == 2) output.bloch_vector.ptr(j, 2).* = density_matrix.at(1, 1).re - density_matrix.at(0, 0).re;
                if (nstate == 2) output.bloch_vector.ptr(j, 3).* = std.math.sqrt(std.math.pow(T, output.bloch_vector.at(j, 0), 2) + std.math.pow(T, output.bloch_vector.at(j, 1), 2));
            }

            if (nstate == 2 and opt.write.spatial_bloch_vector != null) try assignSpatialBlochStep(T, &spatial_bloch.?, &wavefunction, opt.potential, j, time, opt.adiabatic, opt.fix_gauge, allocator);
            if (opt.write.wavefunction != null and i == n_dynamics - 1) try assignWavefunctionStep(T, &wavefunction_dynamics.?, wavefunction, opt.potential, j, time, opt.adiabatic, opt.fix_gauge, allocator);

            if (!enable_printing or (j > 0 and j % opt.log_intervals.iteration != 0)) continue;

            const iteration_info = Custom(T).IterationInfo{
                .density_matrix = density_matrix,
                .population_loss = wavefunction.population_loss,
                .iteration = j,
                .kinetic_energy = kinetic_energy,
                .momentum = momentum,
                .position = position,
                .potential_energy = potential_energy,
                .time = time
            };

            try printIterationInfo(T, iteration_info, &timer);
        }

        if (opt.imaginary) |f| if (i < f.states - 1) {
            try optimized_wavefunctions.append(allocator, try wavefunction.clone(allocator));
        };
    }

    for (optimized_wavefunctions.items) |wf| wf.deinit(allocator);

    if (enable_printing and nstate > 1) for (0..nstate) |i| {

        const tag = if (opt.adiabatic) "ADIABATIC" else "DIABATIC";

        try print("{s}FINAL {s} POPULATION OF STATE {d}: {d:.6}\n", .{if (i == 0) "\n" else "", tag, i, output.population.at(output.population.rows - 1, i)});
    };

    const end_time = @as(T, @floatFromInt(opt.iterations)) * opt.time_step;

    if (opt.write.spectrum) |path| {

        if (enable_printing) try print("\nCALCULATING ENERGY SPECTRUM: ", .{});

        const spectrum = try energySpectrum(T, acf, opt.spectrum, allocator); defer spectrum.deinit(allocator);

        if (enable_printing) try print("DONE\n", .{});

        const nyquist_frequency = std.math.pi / opt.time_step;

        try exportRealMatrixWithLinspacedLeftColumn(T, path, spectrum.asMatrix(), 0, nyquist_frequency);
    }

    if (opt.write.autocorrelation_function) |path| try exportComplexMatrixWithLinspacedLeftColumn(T, path, acf.asMatrix(), 0, end_time);
    if (opt.write.bloch_vector) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.bloch_vector, 0, end_time);
    if (opt.write.final_wavefunction) |path| try wavefunction.write(path, allocator);
    if (opt.write.kinetic_energy) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.kinetic_energy.asMatrix(), 0, end_time);
    if (opt.write.momentum) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.momentum, 0, end_time);
    if (opt.write.population) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.population, 0, end_time);
    if (opt.write.position) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.position, 0, end_time);
    if (opt.write.potential_energy) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.potential_energy.asMatrix(), 0, end_time);
    if (opt.write.total_energy) |path| try exportRealMatrixWithLinspacedLeftColumn(T, path, output.total_energy.asMatrix(), 0, end_time);
    if (opt.write.wavefunction) |path| try exportRealMatrix(T, path, wavefunction_dynamics.?);
    if (opt.write.spatial_bloch_vector) |path| try exportRealMatrix(T, path, spatial_bloch.?);

    if (spatial_bloch != null) spatial_bloch.?.deinit(allocator);
    if (wavefunction_dynamics != null) wavefunction_dynamics.?.deinit(allocator);

    return output;
}

/// Assign current spatial Bloch vector to the spatial Bloch dynamics matrix.
pub fn assignSpatialBlochStep(comptime T: type, spatial_bloch: *RealMatrix(T), wavefunction: *GridWavefunction(T), potential: ElectronicPotential(T), iter: usize, time: T, adiabatic: bool, fix_gauge: bool, allocator: std.mem.Allocator) !void {
    var density_matrix_row = try ComplexMatrix(T).init(wavefunction.nstate, wavefunction.nstate, allocator); defer density_matrix_row.deinit(allocator);

    for (0..wavefunction.data.rows) |i| {

        try wavefunction.densityAtRow(&density_matrix_row, i, potential, time, adiabatic, fix_gauge);

        spatial_bloch.ptr(i, wavefunction.ndim + 3 * iter + 0).* = 2 * density_matrix_row.at(0, 1).re;
        spatial_bloch.ptr(i, wavefunction.ndim + 3 * iter + 1).* = 2 * density_matrix_row.at(0, 1).im;
        spatial_bloch.ptr(i, wavefunction.ndim + 3 * iter + 2).* = density_matrix_row.at(1, 1).re - density_matrix_row.at(0, 0).re;
    }
}

/// Assign current wavefunction to the wavefunction dynamics matrix.
pub fn assignWavefunctionStep(comptime T: type, wavefunction_dynamics: *RealMatrix(T), wavefunction: GridWavefunction(T), potential: ElectronicPotential(T), iter: usize, time: T, adiabatic: bool, fix_gauge: bool, allocator: std.mem.Allocator) !void {
    var diabatic_potential = try RealMatrix(T).init(wavefunction.nstate, wavefunction.nstate, allocator); defer diabatic_potential.deinit(allocator);
    var adiabatic_potential = try RealMatrix(T).init(wavefunction.nstate, wavefunction.nstate, allocator); defer adiabatic_potential.deinit(allocator);
    var adiabatic_eigenvectors = try RealMatrix(T).init(wavefunction.nstate, wavefunction.nstate, allocator); defer adiabatic_eigenvectors.deinit(allocator);
    var previous_eigenvectors = try RealMatrix(T).init(wavefunction.nstate, wavefunction.nstate, allocator); defer previous_eigenvectors.deinit(allocator);

    var position_at_row = try RealVector(T).init(wavefunction.ndim, allocator); defer position_at_row.deinit(allocator);

    var wavefunction_row = try ComplexMatrix(T).init(wavefunction.nstate, 1, allocator); defer wavefunction_row.deinit(allocator);

    for (0..wavefunction.data.rows) |i| {

        for (0..wavefunction.nstate) |j| wavefunction_row.ptr(j, 0).* = wavefunction.data.at(i, j);

        if (adiabatic) {

            try adiabatic_eigenvectors.copyTo(&previous_eigenvectors);

            positionAtRow(T, &position_at_row, i, wavefunction.ndim, wavefunction.npoint, wavefunction.limits);

            try potential.evaluateEigensystem(&diabatic_potential, &adiabatic_potential, &adiabatic_eigenvectors, position_at_row, time);

            if (i > 0 and fix_gauge) try fixGauge(T, &adiabatic_eigenvectors, previous_eigenvectors);

            try mm(T, &wavefunction_row, adiabatic_eigenvectors, true, wavefunction.data.row(i).asMatrix(), false);
        }

        for (0..wavefunction.data.cols) |j| {
            wavefunction_dynamics.ptr(i, wavefunction.ndim + 2 * wavefunction.nstate * iter + 2 * j + 0).* = wavefunction_row.at(j, 0).re;
            wavefunction_dynamics.ptr(i, wavefunction.ndim + 2 * wavefunction.nstate * iter + 2 * j + 1).* = wavefunction_row.at(j, 0).im;
        }
    }
}

/// Energy spectrum from autocorrelation function.
pub fn energySpectrum(comptime T: type, acf: ComplexVector(T), spectrum_options: anytype, allocator: std.mem.Allocator) !RealVector(T) {
    const spectrum_len = try std.math.powi(usize, 2, std.math.log2_int_ceil(usize, acf.len) + spectrum_options.padding_order);

    var spectrum_complex = try ComplexVector(T).initZero(spectrum_len, allocator); defer spectrum_complex.deinit(allocator);

    for (0..acf.len) |i| {
        spectrum_complex.ptr(spectrum_complex.len / 2 + i).* = acf.at(i);
        spectrum_complex.ptr(spectrum_complex.len / 2 - i).* = acf.at(i).conjugate();
    }

    const sigma = @as(T, @floatFromInt(acf.len)) / 4;

    if (spectrum_options.window == .gaussian) for (0..spectrum_complex.len) |i| {

        const dist = @as(T, @floatFromInt(i)) - @as(T, @floatFromInt(spectrum_complex.len / 2));

        const exponent = dist * dist / (2 * sigma * sigma);

        spectrum_complex.ptr(i).* = spectrum_complex.at(i).mul(Complex(T).init(std.math.exp(-exponent), 0));
    };

    for (0..spectrum_complex.len / 2) |i| {
        std.mem.swap(Complex(T), spectrum_complex.ptr(i), spectrum_complex.ptr(i + spectrum_complex.len / 2));
    }

    var spectrum_complex_strided = spectrum_complex.asStrided();

    try cfft1(T, &spectrum_complex_strided, 1);

    var spectrum = try RealVector(T).init(spectrum_complex.len / 2 + 1, allocator);

    for (0..spectrum.len) |i| {
        spectrum.ptr(i).* = spectrum_complex.at(i).re;
    }

    var max_value: T = 0;

    for (0..spectrum.len) |i| {
        if (spectrum.at(i) > max_value) {
            max_value = spectrum.at(i);
        }
    }

    spectrum.divs(max_value);

    return spectrum;
}

/// Initialize the container for the spatial Bloch vector.
pub fn initializeTemporalGridContainer(comptime T: type, wavefunction: GridWavefunction(T), iterations: usize, elements: usize, allocator: std.mem.Allocator) !RealMatrix(T) {
    var container: RealMatrix(T) = try RealMatrix(T).init(wavefunction.data.rows, wavefunction.ndim + elements * (iterations + 1), allocator);

    var position_at_row = try RealVector(T).init(wavefunction.ndim, allocator); defer position_at_row.deinit(allocator);

    for (0..wavefunction.data.rows) |i| {

        positionAtRow(T, &position_at_row, i, wavefunction.ndim, wavefunction.npoint, wavefunction.limits);

        for (0..wavefunction.ndim) |j| {
            container.ptr(i, j).* = position_at_row.at(j);
        }
    }

    return container;
}

/// Print header for iteration info.
pub fn printIterationHeader(ndim: usize, nstate: usize) !void {
    var buffer: [WRITE_BUFFER_SIZE]u8 = undefined;

    var writer = std.io.Writer.fixed(&buffer);

    try writer.print("\n{s:8} ", .{"ITER"});
    try writer.print("{s:12} {s:12} {s:12} ", .{"KINETIC", "POTENTIAL", "TOTAL"});
    try writer.print("{[value]s:[width]} ", .{.value = "POSITION", .width = 10 * ndim + 2 * (ndim - 1) + 2});
    try writer.print("{[value]s:[width]} ", .{.value = "MOMENTUM", .width = 10 * ndim + 2 * (ndim - 1) + 2});
    try writer.print("{[value]s:[width]} ", .{.value = "POPULATION", .width = 9 * nstate + 2 * (nstate - 1) + 2});
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

    for (0..info.position.len) |i| {
        try writer.print("{d:10.4}{s}", .{info.position.at(i), if (i == info.position.len - 1) "" else ", "});
    }

    try writer.print("] [", .{});

    for (0..info.momentum.len) |i| {
        try writer.print("{d:10.4}{s}", .{info.momentum.at(i), if (i == info.momentum.len - 1) "" else ", "});
    }

    try writer.print("] [", .{});

    for (0..info.density_matrix.rows) |i| {
        try writer.print("{d:9.4}{s}", .{info.density_matrix.at(i, i).re + info.population_loss.at(i), if (i == info.density_matrix.rows - 1) "" else ", "});
    }

    try writer.print("] {D}", .{timer.read()}); timer.reset();

    try print("{s}\n", .{writer.buffered()});
}

test "Exact Dynamics on 1D Harmonic Potential" {
    const opt = Options(f64){
        .grid = .{
            .limits = &.{&.{-8, 8}},
            .points = 64
        },
        .initial_conditions = .{
            .mass = 1,
            .momentum = &.{0},
            .position = &.{1},
            .state = 0,
            .gamma = &.{2}
        },
        .potential = .{
            .harmonic = HarmonicPotential(f64){}
        },
        .iterations = 1000,
        .time_step = 0.1
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(output.kinetic_energy.last(), 0.52726330098766, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.potential_energy.last(), 0.59766836993815, TEST_TOLERANCE);
}

test "Exact Dynamics on 2D Harmonic Potential" {
    const opt = Options(f64){
        .grid = .{
            .limits = &.{&.{-8, 8}, &.{-8, 8}},
            .points = 64
        },
        .initial_conditions = .{
            .mass = 1,
            .momentum = &.{0, 0},
            .position = &.{1, 1},
            .state = 0,
            .gamma = &.{2, 2}
        },
        .potential = .{
            .harmonic = HarmonicPotential(f64){
                .k = &.{1, 1}
            }
        },
        .iterations = 1000,
        .time_step = 0.1
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(output.kinetic_energy.last(), 1.05452660197613, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.potential_energy.last(), 1.19533673987723, TEST_TOLERANCE);
}

test "Exact Nonadiabatic Dynamics on Tully's First Potential" {
    const opt = Options(f64){
        .grid = .{
            .limits = &.{&.{-24, 32}},
            .points = 2048
        },
        .initial_conditions = .{
            .mass = 2000,
            .momentum = &.{15},
            .position = &.{-10},
            .state = 1,
            .gamma = &.{2}
        },
        .potential = .{
            .tully_1 = TullyPotential1(f64){}
        },
        .iterations = 350,
        .time_step = 10
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(output.kinetic_energy.last(), 0.06471011654226, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.population.at(opt.iterations, 0), 0.58949426578088, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.population.at(opt.iterations, 1), 0.41050573421579, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.potential_energy.last(), 0.00178988577130, TEST_TOLERANCE);
}

test "Imaginary Time Propagation on 1D Harmonic Potential" {
    const opt = Options(f64){
        .grid = .{
            .limits = &.{&.{-8, 8}},
            .points = 64
        },
        .initial_conditions = .{
            .mass = 1,
            .momentum = &.{0},
            .position = &.{1},
            .state = 0,
            .gamma = &.{2}
        },
        .potential = .{
            .harmonic = HarmonicPotential(f64){}
        },
        .iterations = 1000,
        .time_step = 0.1,
        .imaginary = .{}
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(output.kinetic_energy.last(), 0.25031230493126, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.potential_energy.last(), 0.24968808471946, TEST_TOLERANCE);
}

test "Imaginary Time Propagation on 2D Harmonic Potential" {
    const opt = Options(f64){
        .grid = .{
            .limits = &.{&.{-8, 8}, &.{-8, 8}},
            .points = 64
        },
        .initial_conditions = .{
            .mass = 1,
            .momentum = &.{0, 0},
            .position = &.{1, 0},
            .state = 0,
            .gamma = &.{2, 2}
        },
        .potential = .{
            .harmonic = HarmonicPotential(f64){
                .k = &.{1, 1}
            }
        },
        .iterations = 1000,
        .time_step = 0.1,
        .imaginary = .{}
    };

    const output = try run(f64, opt, false, std.testing.allocator); defer output.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(output.kinetic_energy.last(), 0.50062460986252, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.potential_energy.last(), 0.49937616943892, TEST_TOLERANCE);
}
