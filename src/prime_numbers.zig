//! File with the prime number generation target as well as related functions.

const std = @import("std");

const device_write = @import("device_write.zig");
const real_vector = @import("real_vector.zig");

const RealVector = real_vector.RealVector;

const exportRealVector = device_write.exportRealVector;
const print = device_write.print;
const printJson = device_write.printJson;

/// Options for the prime number generation target.
pub fn Options(comptime _: type) type {
    return struct {
        pub const Filter = enum {
            all,
            mersenne
        };
        pub const Out = struct {
            interval: ?u32 = null,
            path: []const u8
        };

        count: u32 = 10,
        log_interval: u32 = 1,
        start: u32 = 2,

        filter: Filter = .all,
        output: ?Out = null
    };
}

/// Output structure for the prime number generation target.
pub fn Output(comptime T: type) type {
    return struct {
        prime_numbers: RealVector(T),

        allocator: std.mem.Allocator,

        /// Initialize the output structure.
        pub fn init(count: usize, allocator: std.mem.Allocator) !@This() {
            return @This(){
                .prime_numbers = try RealVector(T).init(count, allocator),
                .allocator = allocator
            };
        }

        /// Free the output structure.
        pub fn deinit(self: @This()) void {
            self.prime_numbers.deinit();
        }
    };
}

/// Run the prime number generation target.
pub fn run(comptime T: type, opt: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(opt);

    if (enable_printing) try print("\n{s:9} {s:18}\n", .{"INDEX", "PRIME NUMBER"});

    const max_output = if (opt.output != null and opt.output.?.interval != null) @min(opt.count, opt.output.?.interval.?) else opt.count;

    var output = try Output(T).init(max_output, allocator);

    output.prime_numbers.ptr(0).* = if (opt.start <= 2) 2 else switch (opt.filter) {
        .all => nextPrime(T, opt.start),
        .mersenne => try nextMersenne(T, opt.start, allocator),
    };

    for (0..opt.count) |i| {

        if (i > 0 or output.prime_numbers.len == 0) output.prime_numbers.ptr(i % max_output).* = switch (opt.filter) {
            .all => nextPrime(T, output.prime_numbers.at((i - 1) % max_output)),
            .mersenne => try nextMersenne(T, output.prime_numbers.at((i - 1) % max_output), allocator),
        };

        if (enable_printing and (i == 0 or (i + 1) % opt.log_interval == 0)) {
            try print("{d:9} {d:18}\n", .{i + 1, output.prime_numbers.at(i % max_output)});
        }

        if ((i + 1) % max_output == 0) {
            if (opt.output != null) {

                var buf_start: [32]u8 = undefined; var buf_end: [32]u8 = undefined;

                const start = try std.fmt.bufPrint(&buf_start, "{d}", .{i + 2 - max_output}); const end = try std.fmt.bufPrint(&buf_end, "{d}", .{i + 1});

                const path = if (opt.output.?.interval != null) try std.mem.concat(allocator, u8, &.{opt.output.?.path, "_", start, "-", end}) else opt.output.?.path;

                try exportRealVector(T, path, output.prime_numbers); if (opt.output.?.interval != null) allocator.free(path);
            }
        }
    }

    return output;
}

/// Check if the number is prime.
pub fn isPrime(comptime T: type, p: T) bool {
    return trialDivision(T, p);
}

/// Check if the number is Mersenne prime.
pub fn isMersenne(comptime T: type, p: T, allocator: std.mem.Allocator) !bool {
    return try lucasLehmer(T, p, allocator);
}

/// The Lucas-Lehmer test for Mersenne primes.
pub fn lucasLehmer(comptime T: type, p: T, allocator: std.mem.Allocator) !bool {
    if (!isPrime(T, p)) return false;

    var M = try std.math.big.int.Managed.initSet(allocator, 2); defer M.deinit();
    var s = try std.math.big.int.Managed.initSet(allocator, 4); defer s.deinit();

    var q = try std.math.big.int.Managed.init(allocator); defer q.deinit();

    try std.math.big.int.Managed.shiftLeft(&M, &M, @as(u32, @intCast(p - 1)));
    try std.math.big.int.Managed.addScalar(&M, &M, @as(i32, @intCast(0 - 1)));

    for (0..@as(usize, @intCast(p - 2))) |_| {

        try std.math.big.int.Managed.mul(&s, &s, &s);

        try std.math.big.int.Managed.addScalar(&s, &s, -2);

        try std.math.big.int.Managed.divFloor(&q, &s, &s, &M);
    }

    return std.math.big.int.Managed.eqlZero(s) or p == 2;
}

/// Get the next Mersenne prime number after a given number.
pub fn nextMersenne(comptime T: type, p: T, allocator: std.mem.Allocator) !T {
    if (p < 2) return 2;

    if (p == 2) return 3;

    var candidate: T = if (p % 2 == 0) p + 1 else p + 2;

    while (true) : (candidate += 2) {
        if (try isMersenne(T, candidate, allocator)) return candidate;
    }
}

/// Get the next prime number after a given number.
pub fn nextPrime(comptime T: type, p: T) T {
    if (p < 2) return 2;

    if (p == 2) return 3;

    var candidate: T = if (p % 2 == 0) p + 1 else p + 2;

    while (true) : (candidate += 2) {
        if (isPrime(T, candidate)) return candidate;
    }
}

/// Check if a number is prime using trial division.
pub fn trialDivision(comptime T: type, p: T) bool {
    if (p < 2) return false;

    if (p == 2) return true;

    if (p % 2 == 0) return false;

    var divisor: T = 3;

    while (divisor * divisor <= p) : (divisor += 2) {
        if (p % divisor == 0) return false;
    }

    return true;
}
