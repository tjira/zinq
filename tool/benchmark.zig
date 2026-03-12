const builtin = @import("builtin");
const std = @import("std");
const zinq = @import("zinq");

const RealMatrix = zinq.real_matrix.RealMatrix;

const eigensystemJacobi = zinq.eigenproblem_solver.eigensystemJacobi;
const exportRealMatrix = zinq.device_write.exportRealMatrix;
const mean = zinq.array_functions.mean;
const mm = zinq.matrix_multiplication.mm;
const print = zinq.device_write.print;
const sd = zinq.array_functions.sd;
const timestamp = zinq.timestamp.timestamp;

const SAMPLES = 100;

pub fn EigensystemJacobiBenchStruct(comptime T: type) type {
    return struct {
        A: RealMatrix(T),
        B: RealMatrix(T),
        C: RealMatrix(T),

        pub const name = "JACOBI EIGENSOLVER";

        pub const min: usize =   2;
        pub const max: usize = 100;

        pub fn init(allocator: std.mem.Allocator, size: usize) !@This() {
            return .{
                .A = try RealMatrix(T).init(size, size, allocator),
                .B = try RealMatrix(T).init(size, size, allocator),
                .C = try RealMatrix(T).init(size, size, allocator)
            };
        }

        pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
            self.A.deinit(allocator);
            self.B.deinit(allocator);
            self.C.deinit(allocator);
        }

        pub fn run(self: *@This(), seed: usize) !u64 {
            self.A.randomize(seed); try self.A.symmetrize();

            var timer = try std.time.Timer.start();

            try eigensystemJacobi(RealMatrix, T, &self.B, &self.C, self.A);

            return timer.read();
        }
    };
}

pub fn MmBenchStruct(comptime T: type) type {
    return struct {
        A: RealMatrix(T),
        B: RealMatrix(T),
        C: RealMatrix(T),

        pub const name = "MATRIX MULTIPLICATION";

        pub const min: usize =   2;
        pub const max: usize = 300;

        pub fn init(allocator: std.mem.Allocator, size: usize) !@This() {
            return .{
                .A = try RealMatrix(T).init(size, size, allocator),
                .B = try RealMatrix(T).init(size, size, allocator),
                .C = try RealMatrix(T).init(size, size, allocator)
            };
        }

        pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
            self.A.deinit(allocator);
            self.B.deinit(allocator);
            self.C.deinit(allocator);
        }

        pub fn run(self: *@This(), seed: usize) !u64 {
            self.A.randomize(seed); self.B.randomize(seed);

            var timer = try std.time.Timer.start();

            try mm(T, &self.C, self.A, false, self.B, false);

            return timer.read();
        }
    };
}

pub fn main() !void {
    var gpa = std.heap.GeneralPurposeAllocator(.{}){}; const allocator = gpa.allocator(); const T = f64;

    defer {
        if (gpa.deinit() == .leak) std.log.err("MEMORY LEAK DETECTED IN THE ALLOCATOR\n", .{});
    }

    try print("ZIG VERSION: {d}.{d}.{d}, ZINQ VERSION: {s}", .{builtin.zig_version.major, builtin.zig_version.minor, builtin.zig_version.patch, zinq.config.zinq_version});

    {
        const ts = try timestamp(allocator); defer allocator.free(ts);

        try print(", TIMESTAMP: {s}\n", .{ts});
    }

    const benchmarks = .{
         EigensystemJacobiBenchStruct(T),
         MmBenchStruct(T)
    };

    var timer = try std.time.Timer.start();

    inline for (benchmarks) |bench_struct| {

        const results = try benchmark(T, bench_struct, SAMPLES, allocator); defer results.deinit(allocator);

        const fname = try std.fmt.allocPrint(allocator, "{s}.mat", .{bench_struct.name}); defer allocator.free(fname);

        std.mem.replaceScalar(u8, fname, ' ', '_');

        try exportRealMatrix(T, fname, results);
    }

    try print("\nBENCHMARKING COMPLETED IN {D}\n", .{timer.read()});
}

pub fn sortAndAnalyzeTimings(timings: []f64) !struct{min: u64, median: u64, mean: u64, max: u64, sd: u64} {
    std.mem.sort(f64, timings, {}, comptime std.sort.asc(f64));

    const timings_mean = mean(f64, timings);
    const timings_sd = sd(f64, timings);
    const timings_median = if (timings.len % 2 == 1) timings[timings.len / 2] else (timings[timings.len / 2 - 1] + timings[timings.len / 2]) / 2;
    const timings_min = timings[0];
    const timings_max = timings[timings.len - 1];

    return .{
        .min = @as(u64, @intFromFloat(@round(timings_min))),
        .median = @as(u64, @intFromFloat(@round(timings_median))),
        .mean = @as(u64, @intFromFloat(@round(timings_mean))),
        .max = @as(u64, @intFromFloat(@round(timings_max))),
        .sd = @as(u64, @intFromFloat(@round(timings_sd)))
    };
}

pub fn benchmark(comptime T: type, bench_type: anytype, samples: usize, allocator: std.mem.Allocator) !RealMatrix(T) {
    if (bench_type.min < 2 or bench_type.max < bench_type.min or samples == 0) return error.InvalidInput;

    var results = try RealMatrix(T).init(bench_type.max - bench_type.min + 1, 6, allocator); errdefer results.deinit(allocator);

    try print("\n{s} BENCHMARK FROM {d}x{d} TO {d}x{d}, {d} SAMPLES EACH\n", .{bench_type.name, bench_type.min, bench_type.min, bench_type.max, bench_type.max, samples});

    try print("|---------|----------|-----------|---------|----------|---------|\n", .{});
    try print("|  SHAPE  | MIN TIME |MEDIAN TIME|MEAN TIME| MAX TIME | STD DEV |\n", .{});
    try print("|---------|----------|-----------|---------|----------|---------|\n", .{});

    var timings = try allocator.alloc(f64, samples); defer allocator.free(timings);

    for (bench_type.min..bench_type.max + 1) |i| {

        var bench_struct = try bench_type.init(allocator, i); defer bench_struct.deinit(allocator);

        for (0..samples) |j| timings[j] = @as(f64, @floatFromInt(try bench_struct.run(j)));

        const analysis = try sortAndAnalyzeTimings(timings);

        results.ptr(i - bench_type.min, 0).* = @as(T, @floatFromInt(i              ));
        results.ptr(i - bench_type.min, 1).* = @as(T, @floatFromInt(analysis.min   ));
        results.ptr(i - bench_type.min, 2).* = @as(T, @floatFromInt(analysis.median));
        results.ptr(i - bench_type.min, 3).* = @as(T, @floatFromInt(analysis.mean  ));
        results.ptr(i - bench_type.min, 4).* = @as(T, @floatFromInt(analysis.max   ));
        results.ptr(i - bench_type.min, 5).* = @as(T, @floatFromInt(analysis.sd    ));

        for (1..6) |col| results.ptr(i - bench_type.min, col).* /= std.time.ns_per_ms;

        const print_payload = .{.value = i, .min = analysis.min, .median = analysis.median, .mean = analysis.mean, .max = analysis.max, .sd = analysis.sd};

        try print("|{[value]d:>4}x{[value]d:<4}|{[min]D:^10.3}|{[median]D:^11.3}|{[mean]D:^9.3}|{[max]D:^10.3}|{[sd]D:^9.3}|\n", print_payload);
    }

    try print("|---------|----------|-----------|---------|----------|---------|\n", .{});

    return results;
}
