//! Wrapper around FFTW for performing fast Fourier transforms on multi-dimensional complex arrays.

const std = @import("std");

const fftw = @import("cimport.zig").fftw;

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;
const Matrix = @import("tensor.zig").Matrix;
const RandomDistribution = @import("tensor.zig").RandomDistribution;
const Tensor = @import("tensor.zig").Tensor;

const isComplex = @import("value.zig").isComplex;
const primType = @import("value.zig").primType;
const printf = @import("read_write.zig").printf;
const printComplexTensor = @import("read_write.zig").printComplexTensor;
const readComplexTensor = @import("read_write.zig").readComplexTensor;
const readTensor = @import("read_write.zig").readTensor;
const readTensorRank = @import("read_write.zig").readTensorRank;
const writeComplexTensor = @import("read_write.zig").writeComplexTensor;

var fftw_lock = std.atomic.Value(bool).init(false);

/// Configuration parameters for loading, executing, and saving discrete Fourier transforms.
pub const Options = struct {
    complex: bool = false,
    direction: Direction = .forward,
    log: FftLog = .{},
    tensor: TensorSource,
    write: FftWrite = .{},
};

/// Discrete Fourier transform direction specifier for forward or backward transformation.
pub const Direction = enum(i32) {
    forward = -1,
    backward = 1,
};

/// Flags specifying terminal logging verbosity for Fourier transformed tensors.
pub const FftLog = struct {
    tensor: bool = false,
};

/// Destination paths for persisting computed Fourier transformed tensors to disk.
pub const FftWrite = struct {
    tensor: ?[]const u8 = null,
};

/// Parameters for pseudorandom tensor sampling across arbitrary multi-dimensional shapes.
pub const RandomOptions = struct {
    distribution: RandomDistribution = .{ .normal = .{} },
    seed: u64 = 0,
    shape: []const usize,
};

/// Input data source for loading multi-dimensional arrays from disk or procedural generation.
pub const TensorSource = union(enum) {
    file: []const u8,
    random: RandomOptions,
};

/// Returns an FFT plan type that encapsulates execution details for discrete Fourier transforms.
pub fn FftPlan(comptime T: type) type {
    return struct {
        plan: fftw.fftw_plan,
        sign: i32,

        /// Initializes a discrete Fourier transform plan for the specified shape, direction, and planning mode.
        pub fn init(arr: []T, shape: []i32, sign: i32, mode: u32) !@This() {
            if (comptime primType(T) != f64) @compileError("FFTW NOW ONLY SUPPORTS F64 NUMBERS");

            const ptr = @as([*c]fftw.fftw_complex, @ptrCast(arr.ptr));

            lockFftw();

            const plan = fftw.fftw_plan_dft(@intCast(shape.len), shape.ptr, ptr, ptr, sign, mode);

            unlockFftw();

            return .{ .plan = plan orelse return error.PlanCreationFailed, .sign = sign };
        }

        /// Destroys the FFTW plan, releasing all internal resources and plans.
        pub fn deinit(self: @This()) void {
            lockFftw();

            fftw.fftw_destroy_plan(self.plan);

            unlockFftw();
        }

        /// Creates a duplicate of the existing FFTW plan with identical transform properties.
        pub fn clone(self: @This()) !@This() {
            lockFftw();

            const plan = fftw.fftw_copy_plan(self.plan);

            unlockFftw();

            return .{ .plan = plan orelse return error.PlanDuplicationFailed, .sign = self.sign };
        }

        /// Executes the DFT in-place and applies normalization (1/N) if computing an inverse transform.
        pub fn execute(self: @This(), arr: []T) void {
            const ptr = @as([*c]fftw.fftw_complex, @ptrCast(arr.ptr));

            fftw.fftw_execute_dft(self.plan, ptr, ptr);

            if (self.sign == 1) for (0..arr.len) |i| {
                arr[i].re /= @as(f64, @floatFromInt(arr.len));
                arr[i].im /= @as(f64, @floatFromInt(arr.len));
            };
        }
    };
}

/// Generic container holding the computed Fourier transformed tensor.
pub fn Result(comptime T: type) type {
    return struct {
        matrix: Matrix(Complex(primType(T))),

        /// Deallocates memory associated with the result matrix.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.matrix.deinit(gpa);
        }
    };
}

/// Locks the FFTW library for thread-safe operations, preventing concurrent access to shared resources.
pub fn lockFftw() void {
    while (fftw_lock.swap(true, .acquire)) {
        std.Thread.yield() catch {};
    }
}

/// Executes discrete Fourier transform on multi-dimensional tensors specified by options.
pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator) !Result(T) {
    const rank = switch (opt.tensor) {
        .file => |path| try readTensorRank(io, path),
        .random => |rand_opt| rand_opt.shape.len,
    };

    switch (rank) {
        inline 1...5 => |N| return try runRank(T, N, io, opt, log, gpa),
        else => return error.UnsupportedRank,
    }
}

/// Unlocks the FFTW library, allowing other threads to access shared resources after completing operations.
pub fn unlockFftw() void {
    fftw_lock.store(false, .release);
}

/// Initializes a complex tensor from an input file or procedurally sampled distribution.
fn initComplexTensor(comptime T: type, comptime N: usize, io: std.Io, source: TensorSource, complex: bool, gpa: Allocator) !Tensor(Complex(primType(T)), N) {
    switch (source) {
        .file => |path| {
            if (complex) {
                return try readComplexTensor(T, N, io, path, gpa);
            }

            var real_ten = try readTensor(primType(T), N, io, path, gpa);
            defer real_ten.deinit(gpa);

            var A = try Tensor(Complex(primType(T)), N).init(real_ten.shape, gpa);

            for (0..A.data.len) |i| {
                A.data[i] = Complex(primType(T)).init(real_ten.data[i], 0);
            }

            return A;
        },
        .random => |rand_opt| return try randomComplexTensor(T, N, io, rand_opt, complex, gpa),
    }
}

/// Samples random complex tensor values from Gaussian or uniform distributions across shape bounds.
fn randomComplexTensor(comptime T: type, comptime N: usize, io: std.Io, opt: RandomOptions, complex: bool, gpa: Allocator) !Tensor(Complex(primType(T)), N) {
    if (opt.shape.len != N) return error.InvalidRank;

    var shape: [N]usize = undefined;

    inline for (0..N) |i| shape[i] = opt.shape[i];

    var A = try Tensor(Complex(primType(T)), N).init(shape, gpa);
    errdefer A.deinit(gpa);

    const now = std.Io.Timestamp.now(io, .real);

    const seed = if (opt.seed == 0) @as(u64, @truncate(@as(u96, @bitCast(now.nanoseconds)))) else opt.seed;

    var split_mix = std.Random.SplitMix64.init(seed);
    var rng = std.Random.DefaultPrng.init(split_mix.next());

    const random = rng.random();

    switch (opt.distribution) {
        .normal => |norm_opt| {
            for (0..A.data.len) |i| {
                const re = norm_opt.mean + norm_opt.std * random.floatNorm(primType(T));

                const im = if (complex) norm_opt.mean + norm_opt.std * random.floatNorm(primType(T)) else 0;

                A.data[i] = Complex(primType(T)).init(re, im);
            }
        },
        .uniform => |unif_opt| {
            const diff = unif_opt.max - unif_opt.min;

            for (0..A.data.len) |i| {
                const re = unif_opt.min + diff * random.float(primType(T));

                const im = if (complex) unif_opt.min + diff * random.float(primType(T)) else 0;

                A.data[i] = Complex(primType(T)).init(re, im);
            }
        },
    }

    return A;
}

/// Executes the discrete Fourier transform pipeline for an N-dimensional tensor.
fn runRank(comptime T: type, comptime N: usize, io: std.Io, opt: Options, log: bool, gpa: Allocator) !Result(T) {
    if (log) {
        try printf(io, "\nINITIALIZE TENSOR: ", .{});
    }

    var timer = std.Io.Timestamp.now(io, .real);

    var A = try initComplexTensor(T, N, io, opt.tensor, opt.complex, gpa);
    errdefer A.deinit(gpa);

    if (log) {
        try printf(io, "{f}\n", .{timer.untilNow(io, .real)});
    }

    var shape_i32: [N]i32 = undefined;

    inline for (0..N) |i| shape_i32[i] = @intCast(A.shape[i]);

    const dir = @intFromEnum(opt.direction);

    var plan = try FftPlan(Complex(primType(T))).init(A.data, &shape_i32, dir, fftw.FFTW_ESTIMATE);
    defer plan.deinit();

    if (log) {
        try printf(io, "EXECUTE TRANSFORM: ", .{});
    }

    timer = std.Io.Timestamp.now(io, .real);

    plan.execute(A.data);

    if (log) {
        try printf(io, "{f}\n", .{timer.untilNow(io, .real)});
    }

    if (opt.write.tensor) |path| {
        if (log) {
            try printf(io, "WRITE TENSOR DATA: ", .{});
        }

        timer = std.Io.Timestamp.now(io, .real);

        try writeComplexTensor(primType(T), N, io, path, A);

        if (log) {
            try printf(io, "{f}\n", .{timer.untilNow(io, .real)});
        }
    }

    if (opt.log.tensor) {
        try printf(io, "\nFOURIER TRANSFORM:\n", .{});

        try printComplexTensor(primType(T), N, io, A);
    }

    return Result(T){ .matrix = A.asMatrix() };
}
