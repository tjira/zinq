//! Computes discrete Fourier transforms on multi-dimensional complex arrays using Cooley-Tukey algorithm.

const std = @import("std");

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

/// Computes the one-dimensional discrete Fourier transform using the radix-2 Cooley-Tukey algorithm.
pub fn fft(comptime T: type, output: []Complex(T), input: []const Complex(T), comptime exp: i32) void {
    std.debug.assert(output.len == input.len);

    if (output.len <= 1) return;

    if (output.ptr != input.ptr) {
        @memcpy(output, input);
    }

    std.debug.assert(std.math.isPowerOfTwo(output.len));

    const n, const shift = .{ output.len, @as(u6, @intCast(@bitSizeOf(usize) - @ctz(output.len))) };

    for (0..n) |i| {
        const j = @bitReverse(i) >> shift;

        if (i < j) {
            const temp = output[i];

            output[i], output[j] = .{ output[j], temp };
        }
    }

    var len: usize = 2;

    while (len <= n) : (len <<= 1) {
        const angle = (if (comptime exp == -1) -2.0 else 2.0) * std.math.pi / @as(T, @floatFromInt(len));

        var i: usize, const half = .{ 0, len >> 1 };

        const w_step = Complex(T).init(@cos(angle), @sin(angle));

        while (i < n) : (i += len) {
            var w = Complex(T).init(1, 0);

            for (0..half) |j| {
                const u, const v = .{ output[i + j], output[i + j + half].mul(w) };

                output[i + j], output[i + j + half] = .{ u.add(v), u.sub(v) };

                w = w.mul(w_step);
            }
        }
    }

    if (comptime exp == 1) {
        const inv_n = @as(T, 1) / @as(T, @floatFromInt(n));

        for (output) |*val| {
            val.re *= inv_n;
            val.im *= inv_n;
        }
    }
}

/// Computes the multi-dimensional discrete Fourier transform across arbitrary spatial axes using scratch line memory.
pub fn fftn(comptime T: type, output: []Complex(T), input: []const Complex(T), shape: []const usize, comptime exp: i32, scratch: []Complex(T)) void {
    std.debug.assert(output.len == input.len);

    if (shape.len == 0 or output.len == 0) return;

    if (output.ptr != input.ptr) {
        @memcpy(output, input);
    }

    if (shape.len == 1) {
        fft(T, output, output, exp);

        return;
    }

    var stride_after: usize, var d: usize = .{ 1, shape.len };

    while (d > 0) {
        d -= 1;

        const n_d = shape[d];

        if (n_d <= 1) continue;

        var stride_before: usize = 1;

        for (0..d) |k| stride_before *= shape[k];

        if (stride_after == 1) for (0..stride_before) |p| {
            const base = p * n_d;

            fft(T, output[base .. base + n_d], output[base .. base + n_d], exp);
        };

        if (stride_after > 1) {
            std.debug.assert(scratch.len >= n_d);

            for (0..stride_before) |p| {
                const p_offset = p * n_d * stride_after;

                for (0..stride_after) |s| {
                    const base = p_offset + s;

                    for (0..n_d) |k| scratch[k] = output[base + k * stride_after];

                    fft(T, scratch[0..n_d], scratch[0..n_d], exp);

                    for (0..n_d) |k| output[base + k * stride_after] = scratch[k];
                }
            }
        }

        stride_after *= n_d;
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

    for (opt.shape) |dim| if (!std.math.isPowerOfTwo(dim)) return error.InvalidDimension;

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

    for (A.shape) |dim| if (!std.math.isPowerOfTwo(dim)) {
        try printf(io, "TENSOR DIMENSIONS MUST BE POWERS OF TWO\n", .{});

        return error.InvalidDimension;
    };

    const dir = @backingInt(opt.direction);

    if (log) {
        try printf(io, "EXECUTE TRANSFORM: ", .{});
    }

    timer = std.Io.Timestamp.now(io, .real);

    var max_dim: usize = 0;

    for (A.shape) |dim| max_dim = @max(max_dim, dim);

    const scratch = try gpa.alloc(Complex(primType(T)), if (N > 1) max_dim else 0);
    defer gpa.free(scratch);

    if (dir == -1) {
        fftn(primType(T), A.data, A.data, &A.shape, -1, scratch);
    } else {
        fftn(primType(T), A.data, A.data, &A.shape, 1, scratch);
    }

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
