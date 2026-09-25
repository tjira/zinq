//! Multi-dimensional array (tensor, matrix, vector) containers for numerical physics and linear algebra.

const std = @import("std");

const cblas = @import("cimport.zig").cblas;

const Allocator = std.mem.Allocator;

const Value = @import("value.zig").Value;

const eigh = @import("linear_algebra.zig").eigh;
const mm = @import("linear_algebra.zig").mm;
const printf = @import("read_write.zig").printf;
const printMatrix = @import("read_write.zig").printMatrix;
const readMatrix = @import("read_write.zig").readMatrix;
const writeMatrix = @import("read_write.zig").writeMatrix;

/// Configuration options for executing tensor and matrix transformations from data files.
pub const Options = struct {
    operation: Operation,
};

/// Flags for printing computed eigenvalues and eigenvectors to terminal output.
pub const EighLog = struct {
    eigenvalues: bool = false,
    eigenvectors: bool = false,
};

/// Parameters, logging preferences, and output destinations for symmetric eigendecomposition.
pub const EighOptions = struct {
    matrix: []const u8,
    log: EighLog = .{},
    nthreads: u32 = 1,
    write: EighWrite = .{},
};

/// Output target file paths for saving computed eigenvalues and eigenvectors.
pub const EighWrite = struct {
    eigenvalues: ?[]const u8 = null,
    eigenvectors: ?[]const u8 = null,
};

/// Flags for printing computed matrix multiplication quantities to terminal output.
pub const MatmulLog = struct {
    product: bool = false,
};

/// Parameters and output destinations for matrix-matrix multiplication.
pub const MatmulOptions = struct {
    a: []const u8,
    b: []const u8,
    alpha: f64 = 1,
    beta: f64 = 0,
    log: MatmulLog = .{},
    nthreads: u32 = 1,
    trans_a: bool = false,
    trans_b: bool = false,
    write: MatmulWrite = .{},
};

/// Output target file paths for saving matrix multiplication products.
pub const MatmulWrite = struct {
    product: ?[]const u8 = null,
};

/// Mean and standard deviation parameters for Gaussian random number generation.
pub const NormalDistribution = struct {
    mean: f64 = 0,
    std: f64 = 1,
};

/// Tagged union specifying the linear algebra operation to execute.
pub const Operation = union(enum) {
    eigh: EighOptions,
    matmul: MatmulOptions,
    random: RandomOptions,
};

/// Probability distribution specifier for pseudorandom number sampling.
pub const RandomDistribution = union(enum) {
    normal: NormalDistribution,
    uniform: UniformDistribution,
};

/// Flags for printing generated random matrices to terminal output.
pub const RandomLog = struct {
    matrix: bool = false,
};

/// Parameters, distribution type, and output destinations for random matrix generation.
pub const RandomOptions = struct {
    distribution: RandomDistribution = .{ .normal = .{} },
    log: RandomLog = .{},
    seed: u64 = 0,
    shape: [2]usize,
    symmetric: bool = false,
    write: RandomWrite = .{},
};

/// Output target file paths for saving generated random matrices.
pub const RandomWrite = struct {
    matrix: ?[]const u8 = null,
};

/// Interval bounds for uniform pseudorandom number generation.
pub const UniformDistribution = struct {
    max: f64 = 1,
    min: f64 = 0,
};

/// Returns a 2D matrix type representing linear operators or grids in coordinate space.
pub fn Matrix(comptime T: type) type {
    return struct {
        data: []T,
        shape: [2]usize,

        /// Allocates memory for a matrix of size rows * cols.
        pub fn init(rows: usize, cols: usize, gpa: std.mem.Allocator) !@This() {
            return .{ .data = try gpa.alloc(T, rows * cols), .shape = .{ rows, cols } };
        }

        /// Frees the allocated memory of the matrix.
        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            gpa.free(self.data);
        }

        /// Promotes the 2D matrix to a rank-2 tensor representation.
        pub fn asTensor(self: @This()) Tensor(T, 2) {
            return .{ .data = self.data, .shape = self.shape };
        }

        /// Flattens the 2D matrix into a 1D vector representation.
        pub fn asVector(self: @This()) Vector(T) {
            return .{ .data = self.data, .shape = .{self.data.len} };
        }

        /// Returns the matrix element at row i and column j.
        pub fn at(self: @This(), i: usize, j: usize) T {
            std.debug.assert(i < self.shape[0]);
            std.debug.assert(j < self.shape[1]);

            return self.data[i * self.shape[1] + j];
        }

        /// Creates a deep copy of the matrix using the provided allocator.
        pub fn clone(self: @This(), gpa: std.mem.Allocator) !@This() {
            var A = try @This().init(self.shape[0], self.shape[1], gpa);

            for (0..self.data.len) |i| {
                A.data[i] = self.data[i];
            }

            return A;
        }

        /// Divides all elements of the matrix by a scalar in-place.
        pub fn divs(self: *@This(), scalar: T) void {
            for (0..self.data.len) |i| {
                self.data[i] = Value(T).init(self.data[i]).div(Value(T).init(scalar)).val;
            }
        }

        /// Fills the matrix with a uniform scalar value.
        pub fn fill(self: *@This(), scalar: T) void {
            for (0..self.data.len) |i| {
                self.data[i] = scalar;
            }
        }

        /// Wraps an existing slice into a 2D matrix view of size rows * cols.
        pub fn fromSlice(rows: usize, cols: usize, data: []T) @This() {
            std.debug.assert(data.len == rows * cols);

            return .{ .data = data, .shape = .{ rows, cols } };
        }

        /// Allocates a matrix and initializes all elements to zero.
        pub fn initZero(rows: usize, cols: usize, gpa: std.mem.Allocator) !@This() {
            var A = try @This().init(rows, cols, gpa);

            A.zero();

            return A;
        }

        /// Returns the maximum absolute value (infinity norm candidate) of the matrix elements.
        pub fn max(self: @This()) T {
            var max_val: T = 0;

            for (self.data) |x| {
                const abs_val = Value(T).init(x).abs().val;

                if (abs_val > max_val) {
                    max_val = abs_val;
                }
            }

            return max_val;
        }

        /// Returns the number of columns in the matrix.
        pub fn ncol(self: @This()) usize {
            return self.shape[1];
        }

        /// Returns the number of rows in the matrix.
        pub fn nrow(self: @This()) usize {
            return self.shape[0];
        }

        /// Returns a pointer to the matrix element at row i and column j.
        pub fn ptr(self: *@This(), i: usize, j: usize) *T {
            std.debug.assert(i < self.shape[0]);
            std.debug.assert(j < self.shape[1]);

            return &self.data[i * self.shape[1] + j];
        }

        /// Returns the rank (number of dimensions) of the matrix.
        pub fn rank(self: @This()) usize {
            return self.shape.len;
        }

        /// Computes the root-mean-square value of the matrix elements.
        pub fn rms(self: @This()) T {
            var sum_sq: T = 0;

            for (self.data) |x| {
                const abs_val = Value(T).init(x).abs().val;

                sum_sq += abs_val * abs_val;
            }

            return @sqrt(sum_sq / @as(T, @floatFromInt(self.data.len)));
        }

        /// Returns a 1D vector view of the specified row.
        pub fn row(self: @This(), i: usize) Vector(T) {
            std.debug.assert(i < self.shape[0]);

            return .{ .data = self.data[i * self.shape[1] .. (i + 1) * self.shape[1]], .shape = .{self.shape[1]} };
        }

        /// Returns a slice view of the specified row.
        pub fn rowSlice(self: @This(), i: usize) []T {
            std.debug.assert(i < self.shape[0]);

            return self.data[i * self.shape[1] .. (i + 1) * self.shape[1]];
        }

        /// Symmetrizes a square matrix in-place by averaging off-diagonal elements.
        pub fn symmetrize(self: *@This()) void {
            std.debug.assert(self.shape[0] == self.shape[1]);

            for (0..self.shape[0]) |i| for (i + 1..self.shape[1]) |j| {
                const avg = Value(T).init(self.at(i, j)).add(Value(T).init(self.at(j, i))).divs(2).val;

                self.ptr(i, j).* = avg;
                self.ptr(j, i).* = avg;
            };
        }

        /// Returns a submatrix view consisting of the first n rows.
        pub fn takeRows(self: @This(), n: usize) @This() {
            std.debug.assert(n <= self.shape[0]);

            return .{ .data = self.data[0 .. n * self.shape[1]], .shape = .{ n, self.shape[1] } };
        }

        /// Sets all elements of the matrix to zero.
        pub fn zero(self: *@This()) void {
            self.fill(std.mem.zeroes(T));
        }
    };
}

/// Generic container holding computed multi-dimensional tensors or matrices.
pub fn Result(comptime T: type) type {
    return struct {
        tensors: []Matrix(T),

        /// Deallocates memory associated with the result tensors.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            for (self.tensors) |*t| t.deinit(gpa);

            gpa.free(self.tensors);
        }
    };
}

/// Returns a multi-dimensional array type of rank N for physics tensors.
pub fn Tensor(comptime T: type, comptime N: usize) type {
    return struct {
        data: []T,
        shape: [N]usize,

        /// Allocates memory for a tensor with the given multi-index shape.
        pub fn init(shape: [N]usize, gpa: std.mem.Allocator) !@This() {
            var size: usize = 1;

            inline for (0..N) |i| {
                size *= shape[i];
            }

            return .{ .data = try gpa.alloc(T, size), .shape = shape };
        }

        /// Frees the allocated memory of the tensor.
        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            gpa.free(self.data);
        }

        /// Reshapes the N-dimensional tensor into a 2D matrix representation.
        pub fn asMatrix(self: @This()) Matrix(T) {
            var rows: usize = 1;
            var cols: usize = 1;

            for (0..(N + 1) / 2) |i| {
                rows *= self.shape[i];
            }

            for ((N + 1) / 2..N) |i| {
                cols *= self.shape[i];
            }

            return .{ .data = self.data, .shape = .{ rows, cols } };
        }

        /// Promotes the tensor to a tensor representation to satisfy the tensor interface.
        pub fn asTensor(self: @This()) @This() {
            return self;
        }

        /// Returns the tensor element at the specified multi-index coordinate.
        pub fn at(self: @This(), indx: [N]usize) T {
            var idx: usize = 0;
            var str: usize = 1;

            inline for (0..N) |k| {
                const j = N - 1 - k;

                std.debug.assert(indx[j] < self.shape[j]);

                idx += indx[j] * str;
                str *= self.shape[j];
            }

            return self.data[idx];
        }

        /// Allocates a tensor and initializes all elements to zero.
        pub fn initZero(shape: [N]usize, gpa: std.mem.Allocator) !@This() {
            var U = try @This().init(shape, gpa);

            U.zero();

            return U;
        }

        /// Returns a pointer to the tensor element at the specified multi-index coordinate.
        pub fn ptr(self: *@This(), indx: [N]usize) *T {
            var idx: usize = 0;
            var str: usize = 1;

            inline for (0..N) |k| {
                const j = N - 1 - k;

                std.debug.assert(indx[j] < self.shape[j]);

                idx += indx[j] * str;
                str *= self.shape[j];
            }

            return &self.data[idx];
        }

        /// Returns the rank (number of dimensions) of the tensor.
        pub fn rank(self: @This()) usize {
            return self.shape.len;
        }

        /// Sets all elements of the tensor to zero.
        pub fn zero(self: *@This()) void {
            for (0..self.data.len) |i| {
                self.data[i] = std.mem.zeroes(T);
            }
        }
    };
}

/// Returns a 1D vector type representing physical coordinates or state vectors.
pub fn Vector(comptime T: type) type {
    return struct {
        data: []T,
        shape: [1]usize,

        /// Allocates memory for a vector of the specified size.
        pub fn init(size: usize, gpa: std.mem.Allocator) !@This() {
            return .{ .data = try gpa.alloc(T, size), .shape = .{size} };
        }

        /// Frees the allocated memory of the vector.
        pub fn deinit(self: *@This(), gpa: std.mem.Allocator) void {
            gpa.free(self.data);
        }

        /// Promotes the 1D vector to a 2D column matrix (N x 1).
        pub fn asMatrix(self: @This()) Matrix(T) {
            return .{ .data = self.data, .shape = .{ self.shape[0], 1 } };
        }

        /// Promotes the 1D vector to a rank-1 tensor representation.
        pub fn asTensor(self: @This()) Tensor(T, 1) {
            return .{ .data = self.data, .shape = self.shape };
        }

        /// Returns the vector element at index i.
        pub fn at(self: @This(), i: usize) T {
            std.debug.assert(i < self.shape[0]);

            return self.data[i];
        }

        /// Divides all elements of the vector by a scalar in-place.
        pub fn divs(self: *@This(), scalar: T) void {
            for (0..self.data.len) |i| {
                self.data[i] = Value(T).init(self.data[i]).div(Value(T).init(scalar)).val;
            }
        }

        /// Wraps an existing slice into a 1D vector view.
        pub fn fromSlice(data: []T) @This() {
            return .{ .data = data, .shape = .{data.len} };
        }

        /// Allocates a vector and initializes all elements to zero.
        pub fn initZero(size: usize, gpa: std.mem.Allocator) !@This() {
            var v = try @This().init(size, gpa);

            v.zero();

            return v;
        }

        /// Returns the number of elements in the vector.
        pub fn length(self: @This()) usize {
            return self.shape[0];
        }

        /// Multiplies all elements of the vector by a scalar in-place.
        pub fn muls(self: *@This(), scalar: T) void {
            for (0..self.data.len) |i| {
                self.data[i] = Value(T).init(self.data[i]).mul(Value(T).init(scalar)).val;
            }
        }

        /// Returns a pointer to the vector element at index i.
        pub fn ptr(self: *@This(), i: usize) *T {
            std.debug.assert(i < self.shape[0]);

            return &self.data[i];
        }

        /// Returns the rank (number of dimensions) of the vector.
        pub fn rank(self: @This()) usize {
            return self.shape.len;
        }

        /// Returns a subvector view consisting of the first n elements.
        pub fn takeRows(self: @This(), n: usize) @This() {
            std.debug.assert(n <= self.shape[0]);

            return .{ .data = self.data[0..n], .shape = .{n} };
        }

        /// Sets all elements of the vector to zero.
        pub fn zero(self: *@This()) void {
            for (0..self.data.len) |i| {
                self.data[i] = std.mem.zeroes(T);
            }
        }
    };
}

/// Executes tensor or matrix operations specified by options and writes results to files.
pub fn run(comptime T: type, io: std.Io, opt: Options, log: bool, gpa: Allocator) !Result(T) {
    switch (opt.operation) {
        .eigh => |eigh_opt| return try runEigh(T, io, eigh_opt, log, gpa),
        .matmul => |matmul_opt| return try runMatmul(T, io, matmul_opt, log, gpa),
        .random => |rand_opt| return try runRandom(T, io, rand_opt, log, gpa),
    }
}

/// Computes eigenvalues and eigenvectors of a symmetric matrix from an input file.
pub fn runEigh(comptime T: type, io: std.Io, opt: EighOptions, log: bool, gpa: Allocator) !Result(T) {
    if (opt.nthreads == 0) {
        std.log.err("THREAD COUNT MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    cblas.openblas_set_num_threads(@intCast(opt.nthreads));

    if (log) {
        try printf(io, "\nREAD MATRIX: ", .{});
    }

    var timer = std.Io.Timestamp.now(io, .real);

    var A = try readMatrix(T, io, opt.matrix, gpa);
    defer A.deinit(gpa);

    if (A.nrow() != A.ncol()) {
        std.log.err("EIGENVALUE DECOMPOSITION REQUIRES A SQUARE MATRIX", .{});

        return error.InvalidInput;
    }

    if (log) {
        try printf(io, "{f}\n", .{timer.untilNow(io, .real)});
    }

    var W = try Vector(T).init(A.nrow(), gpa);
    errdefer W.deinit(gpa);

    var U = try Matrix(T).init(A.nrow(), A.ncol(), gpa);
    errdefer U.deinit(gpa);

    if (log) {
        try printf(io, "\nCOMPUTE EIGH: ", .{});
    }

    timer = std.Io.Timestamp.now(io, .real);

    try eigh(T, &W, &U, A);

    if (log) {
        try printf(io, "{f}\n", .{timer.untilNow(io, .real)});
    }

    if (opt.write.eigenvalues != null or opt.write.eigenvectors != null) {
        if (log) {
            try printf(io, "\nWRITE MATRIX: ", .{});
        }

        timer = std.Io.Timestamp.now(io, .real);

        if (opt.write.eigenvalues) |path| {
            try writeMatrix(T, io, path, W.asMatrix());
        }

        if (opt.write.eigenvectors) |path| {
            try writeMatrix(T, io, path, U);
        }

        if (log) {
            try printf(io, "{f}\n", .{timer.untilNow(io, .real)});
        }
    }

    if (opt.log.eigenvalues) {
        try printf(io, "\nEIGENVALUES:\n", .{});

        try printMatrix(T, io, W.asMatrix());
    }

    if (opt.log.eigenvectors) {
        try printf(io, "\nEIGENVECTORS:\n", .{});

        try printMatrix(T, io, U);
    }

    const tensors = try gpa.alloc(Matrix(T), 2);

    tensors[0], tensors[1] = .{ W.asMatrix(), U };

    return Result(T){ .tensors = tensors };
}

/// Executes matrix-matrix multiplication on input files using BLAS GEMM and exports the product.
pub fn runMatmul(comptime T: type, io: std.Io, opt: MatmulOptions, log: bool, gpa: Allocator) !Result(T) {
    if (opt.nthreads == 0) {
        std.log.err("THREAD COUNT MUST BE GREATER THAN 0", .{});

        return error.InvalidInput;
    }

    cblas.openblas_set_num_threads(@intCast(opt.nthreads));

    if (log) {
        try printf(io, "\nREAD MATRICES: ", .{});
    }

    var timer = std.Io.Timestamp.now(io, .real);

    var A = try readMatrix(T, io, opt.a, gpa);
    defer A.deinit(gpa);

    var B = try readMatrix(T, io, opt.b, gpa);
    defer B.deinit(gpa);

    if (log) {
        try printf(io, "{f}\n", .{timer.untilNow(io, .real)});
    }

    const m = if (opt.trans_a) A.ncol() else A.nrow();
    const n = if (opt.trans_b) B.nrow() else B.ncol();

    var C = try Matrix(T).initZero(m, n, gpa);
    errdefer C.deinit(gpa);

    if (log) {
        try printf(io, "\nCOMPUTE MATMUL: ", .{});
    }

    timer = std.Io.Timestamp.now(io, .real);

    mm(T, &C, A, B, opt.alpha, opt.beta, opt.trans_a, opt.trans_b);

    if (log) {
        try printf(io, "{f}\n", .{timer.untilNow(io, .real)});
    }

    if (opt.write.product) |path| {
        if (log) {
            try printf(io, "\nWRITE PRODUCT: ", .{});
        }

        timer = std.Io.Timestamp.now(io, .real);

        try writeMatrix(T, io, path, C);

        if (log) {
            try printf(io, "{f}\n", .{timer.untilNow(io, .real)});
        }
    }

    if (opt.log.product) {
        try printf(io, "\nPRODUCT:\n", .{});

        try printMatrix(T, io, C);
    }

    const tensors = try gpa.alloc(Matrix(T), 1);

    tensors[0] = C;

    return Result(T){ .tensors = tensors };
}

/// Generates a random matrix with specified dimensions and probability distribution.
pub fn runRandom(comptime T: type, io: std.Io, opt: RandomOptions, log: bool, gpa: Allocator) !Result(T) {
    if (opt.symmetric and opt.shape[0] != opt.shape[1]) {
        std.log.err("SYMMETRIC MATRIX GENERATION REQUIRES A SQUARE MATRIX", .{});

        return error.InvalidInput;
    }

    var A = try Matrix(T).init(opt.shape[0], opt.shape[1], gpa);
    errdefer A.deinit(gpa);

    if (log) {
        try printf(io, "\nINITIALIZE RANDOM MATRIX: ", .{});
    }

    var timer = std.Io.Timestamp.now(io, .real);

    const seed = if (opt.seed == 0) @as(u64, @truncate(@as(u96, @bitCast(timer.nanoseconds)))) else opt.seed;

    var split_mix = std.Random.SplitMix64.init(seed);

    var rng = std.Random.DefaultPrng.init(split_mix.next());

    const random = rng.random();

    switch (opt.distribution) {
        .normal => |norm_opt| {
            for (0..A.data.len) |i| {
                A.data[i] = norm_opt.mean + norm_opt.std * random.floatNorm(T);
            }
        },
        .uniform => |unif_opt| {
            const diff = unif_opt.max - unif_opt.min;

            for (0..A.data.len) |i| {
                A.data[i] = unif_opt.min + diff * random.float(T);
            }
        },
    }

    if (log) {
        try printf(io, "{f}\n", .{timer.untilNow(io, .real)});
    }

    if (opt.symmetric) {
        if (log) {
            try printf(io, "SYMMETRIZE RANDOM MATRIX: ", .{});
        }

        timer = std.Io.Timestamp.now(io, .real);

        A.symmetrize();

        if (log) {
            try printf(io, "{f}\n", .{timer.untilNow(io, .real)});
        }
    }

    if (opt.write.matrix) |path| {
        if (log) {
            try printf(io, "\nWRITE MATRIX: ", .{});
        }

        timer = std.Io.Timestamp.now(io, .real);

        try writeMatrix(T, io, path, A);

        if (log) {
            try printf(io, "{f}\n", .{timer.untilNow(io, .real)});
        }
    }

    if (opt.log.matrix) {
        try printf(io, "\nMATRIX:\n", .{});

        try printMatrix(T, io, A);
    }

    const tensors = try gpa.alloc(Matrix(T), 1);

    tensors[0] = A;

    return Result(T){ .tensors = tensors };
}
