//! File with the prime number generation target as well as related functions.

const std = @import("std");

const device_write = @import("device_write.zig");
const error_handling = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");
const integer_arithmetic = @import("integer_arithmetic.zig");
const real_vector = @import("real_vector.zig");

const RealVector = real_vector.RealVector;

const addWithOverflow = integer_arithmetic.addWithOverflow;
const exportRealMatrix = device_write.exportRealMatrix;
const print = device_write.print;
const printJson = device_write.printJson;
const printRealMatrix = device_write.printRealMatrix;
const squareWithOverflow = integer_arithmetic.squareWithOverflow;
const throw = error_handling.throw;

const WRITE_BUFFER_SIZE = global_variables.WRITE_BUFFER_SIZE;

/// Options for the prime number generation target.
pub fn Options(comptime _: type) type {
    return struct {
        pub const Filter = enum {
            all,
            mersenne
        };
        pub const Mode = union(enum) {
            check: struct {
                number: []const u8,

                filter: Filter = .all
            },
            factorize: struct {
                number: []const u8
            },
            generate: struct {
                pub const OutputFiles = struct {
                    interval: ?u32 = null,
                    path: []const u8
                };

                count: u32,
                log_interval: u32 = 1,
                start: []const u8 = "2",

                filter: Filter = .all,
                output: ?OutputFiles = null
            }
        };

        bits: ?u32 = 64,

        mode: Mode,
    };
}

/// Output structure for the prime number generation target.
pub fn Output(comptime _: type) type {
    return struct {

        /// Free the output structure.
        pub fn deinit(_: @This(), _: std.mem.Allocator) void {}
    };
}

/// Run the prime number generation target.
pub fn run(comptime T: type, opt: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(opt);

    const bits: i32 = if (opt.bits != null) @intCast(opt.bits.?) else -1;

    return switch (bits) {
        8 => try dispatch(u8, opt, enable_printing, allocator),
        16 => try dispatch(u16, opt, enable_printing, allocator),
        32 => try dispatch(u32, opt, enable_printing, allocator),
        64 => try dispatch(u64, opt, enable_printing, allocator),
        128 => try dispatch(u128, opt, enable_printing, allocator),
        256 => try dispatch(u256, opt, enable_printing, allocator),
        512 => try dispatch(u512, opt, enable_printing, allocator),
        1024 => try dispatch(u1024, opt, enable_printing, allocator),
        -1 => try dispatch(std.math.big.int.Managed, opt, enable_printing, allocator),
        else => return throw(Output(T), "BIT SIZE HAS TO BE POWER OF TWO BETWEEN 8 AND 1024, VALUE '{d}' IS INVALID", .{bits})
    };
}

/// Helper function to avoid repeating the mode switch logic
pub fn dispatch(comptime T: type, opt: anytype, enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) {

        const type_name = if (comptime T == std.math.big.int.Managed) "BIG" else @typeName(T);

        const max_value = if (comptime T == std.math.big.int.Managed) std.math.inf(f128) else @as(f128, @floatFromInt(std.math.maxInt(T)));

        try print("\nUSING TYPE: {s}, MAX VALUE: {e:.3}\n", .{type_name, max_value});
    }

    return switch (opt.mode) {
        .check => try checkMode(T, opt.mode.check, enable_printing, allocator),
        .factorize => try factorizeMode(T, opt.mode.factorize, enable_printing, allocator),
        .generate => try generateMode(T, opt.mode.generate, enable_printing, allocator),
    };
}

/// Check mode.
pub fn checkMode(comptime T: type, opt: anytype, enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    var number_big = try std.math.big.int.Managed.init(allocator); defer number_big.deinit();

    try number_big.setString(10, opt.number); const number = if (comptime T == std.math.big.int.Managed) number_big else try number_big.toInt(T);

    const is_prime = switch (opt.filter) {
        .all => try isPrime(number, allocator),
        .mersenne => try isMersenne(number, allocator),
    };

    if (enable_printing) {

        const s = if (opt.filter == .mersenne) " MERSENNE" else "";

        try print("\n{s} IS{s} A{s} PRIME NUMBER\n", .{opt.number, if (is_prime) "" else " NOT", s});
    }

    return Output(T){};
}

/// Factorize mode.
pub fn factorizeMode(comptime T: type, opt: anytype, enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    var number_big = try std.math.big.int.Managed.init(allocator); defer number_big.deinit();

    try number_big.setString(10, opt.number); const number = if (comptime T == std.math.big.int.Managed) number_big else try number_big.toInt(T);

    const factors = try factorize(number, allocator); defer allocator.free(factors);

    var prime_numbers = try allocator.alloc(std.math.big.int.Managed, factors.len); defer allocator.free(prime_numbers);

    for (factors, 0..) |factor, i| {
        prime_numbers[i] = if (comptime T == std.math.big.int.Managed) factor else try std.math.big.int.Managed.initSet(allocator, factor);
    }

    if (enable_printing) {

        if (factors.len == 0) {
            try print("\n{s} HAS NO PRIME FACTORS\n", .{opt.number}); return Output(T){};
        }

        try print("\n", .{});

        for (prime_numbers, 1..) |factor, i| {

            const factor_str = try factor.toString(allocator, 10, std.fmt.Case.lower); defer allocator.free(factor_str);

            try print("FACTOR {[index]d:[width]}: {[string]s}\n", .{.index = i, .width = std.math.log10_int(factors.len) + 1, .string = factor_str});
        }
    }

    return Output(T){};
}

/// Generate mode.
pub fn generateMode(comptime T: type, opt: anytype, enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try print("\n{s:9} {s:20} {s:4}\n", .{"INDEX", "PRIME NUMBER", "TIME"});

    var start_big = try std.math.big.int.Managed.init(allocator); defer start_big.deinit();

    try start_big.setString(10, opt.start); try start_big.addScalar(&start_big, -1);

    const start = if (comptime T == std.math.big.int.Managed) start_big else try start_big.toInt(T);

    const max_output = if (opt.output != null and opt.output.?.interval != null) @min(opt.count, opt.output.?.interval.?) else opt.count;

    var prime_numbers: []T = try allocator.alloc(T, max_output); defer allocator.free(prime_numbers);

    prime_numbers[0] = switch (opt.filter) {
        .all => try nextPrime(start, allocator),
        .mersenne => try nextMersenne(start, allocator)
    };

    var timer = try std.time.Timer.start();

    for (0..opt.count) |i| {

        if (comptime T == std.math.big.int.Managed) if (i >= max_output) prime_numbers[i % max_output].deinit();

        if (i > 0) prime_numbers[i % max_output] = switch (opt.filter) {
            .all => try nextPrime(prime_numbers[(i - 1) % max_output], allocator),
            .mersenne => try nextMersenne(prime_numbers[(i - 1) % max_output], allocator)
        };

        if (enable_printing and (i == 0 or (i + 1) % opt.log_interval == 0)) {

            const prime_string = if (comptime T == std.math.big.int.Managed) try prime_numbers[i % max_output].toString(allocator, 10, std.fmt.Case.lower) else "";

            const format_string = "{d:9} " ++ (if (comptime T == std.math.big.int.Managed) "{s:20}" else "{d:20}") ++ " {D}\n";

            try print(format_string, .{i + 1, if (comptime T == std.math.big.int.Managed) prime_string else prime_numbers[i % max_output], timer.read()});

            if (comptime T == std.math.big.int.Managed) allocator.free(prime_string);
        }

        if ((i + 1) % max_output == 0) if (opt.output != null) {

            var buf_start: [32]u8 = undefined; var buf_end: [32]u8 = undefined;

            const start_i_str = try std.fmt.bufPrint(&buf_start, "{d}", .{i + 2 - max_output}); const end_i_str = try std.fmt.bufPrint(&buf_end, "{d}", .{i + 1});

            const path = if (opt.output.?.interval != null) try std.mem.concat(allocator, u8, &.{opt.output.?.path, "_", start_i_str, "-", end_i_str}) else opt.output.?.path;

            try exportPrimeNumbersAsRealMatrix(path, prime_numbers, allocator);

            if (opt.output.?.interval != null) allocator.free(path);
        };
    }

    return Output(T){};
}

/// Factorize a number into its prime factors.
pub fn factorize(p: anytype, allocator: std.mem.Allocator) ![]@TypeOf(p) {
    const T = @TypeOf(p); const ap = T == std.math.big.int.Managed;

    var one = if (comptime ap) try std.math.big.int.Managed.initSet(allocator, 1) else @as(T, 1); defer if (comptime ap) one.deinit();
    var two = if (comptime ap) try std.math.big.int.Managed.initSet(allocator, 2) else @as(T, 2); defer if (comptime ap) two.deinit();
    var rem = if (comptime ap) try std.math.big.int.Managed.initSet(allocator, 0) else @as(T, 0); defer if (comptime ap) rem.deinit();
    var pfd = if (comptime ap) try std.math.big.int.Managed.initSet(allocator, 0) else @as(T, 0); defer if (comptime ap) pfd.deinit();

    if (if (comptime ap) p.order(two) == std.math.Order.lt else p < 2) return &[_]T{};

    var factors: std.ArrayList(T) = .empty;

    var n = if (comptime ap) try p.clone() else p; var pd = if (comptime ap) try std.math.big.int.Managed.initSet(allocator, 0) else @as(T, 0);

    while (if (comptime ap) !n.eql(one) else n != 1) {

        var pd_next = try nextPrime(pd, allocator); defer if (comptime ap) pd_next.deinit();

        if (comptime ap) try pd.copy(pd_next.toConst()) else pd = pd_next;

        if (comptime ap) try pfd.divFloor(&rem, &n, &pd);

        while (if (comptime ap) rem.eqlZero() else n % pd == 0) {

            try factors.append(allocator, if (comptime ap) try pd.clone() else pd);

            if (comptime ap) try n.copy(pfd.toConst()) else n /= pd;

            if (comptime ap) try pfd.divFloor(&rem, &n, &pd);

        }
    }

    if (comptime ap) n.deinit(); if (comptime ap) pd.deinit();

    return factors.toOwnedSlice(allocator);
}

/// Check if the number is prime.
pub fn isPrime(p: anytype, allocator: std.mem.Allocator) !bool {
    return trialDivision(p, allocator);
}

/// Check if the number is Mersenne prime.
pub fn isMersenne(p: anytype, allocator: std.mem.Allocator) !bool {
    if (!try isPrime(p, allocator)) return false;

    const T = @TypeOf(p); const ap = T == std.math.big.int.Managed; const p32 = if (comptime ap) try p.toInt(u32) else @as(u32, @intCast(p));

    var M = try std.math.big.int.Managed.initSet(allocator, 2); defer M.deinit();

    try M.shiftLeft(&M, p32 - 1); try M.addScalar(&M, -1);

    return try lucasLehmer(&M, p32, allocator);
}

/// The Lucas-Lehmer test for Mersenne primes.
pub fn lucasLehmer(M: *const std.math.big.int.Managed, p: u32, allocator: std.mem.Allocator) !bool {
    if (p == 2) return true;

    var s = try std.math.big.int.Managed.initSet(allocator, 4); defer s.deinit();
    var q = try std.math.big.int.Managed.initSet(allocator, 0); defer q.deinit();

    try s.ensureCapacity((2 * p + (@bitSizeOf(usize) - 1)) / @bitSizeOf(usize));
    try q.ensureCapacity((2 * p + (@bitSizeOf(usize) - 1)) / @bitSizeOf(usize));

    for (0..p - 2) |_| {

        try s.sqr(&s); try q.shiftRight(&s, p); try s.truncate(&s, .unsigned, p); try s.add(&s, &q);

        if (s.order(M.*) != .lt) try s.sub(&s, M);

        try s.addScalar(&s, -2);

        if (!s.isPositive()) try s.add(&s, M);
    }

    return std.math.big.int.Managed.eqlZero(s);
}

/// Get the next Mersenne prime number after a given number.
pub fn nextMersenne(p: anytype, allocator: std.mem.Allocator) !@TypeOf(p) {
    const T = @TypeOf(p); const ap = T == std.math.big.int.Managed;

    var two = if (comptime ap) try std.math.big.int.Managed.initSet(allocator, 2) else @as(T, 2); defer if (comptime ap) two.deinit();

    if (if (comptime ap) p.order(two) == std.math.Order.lt else p < 2) {
        return if (comptime ap) try std.math.big.int.Managed.initSet(allocator, 2) else 2;
    }

    var candidate = if (comptime ap) try p.clone() else if (p % 2 == 0) p - 1 else p;

    if (comptime ap) if (p.isEven()) try candidate.addScalar(&p, -1);

    while (true) {

        if (comptime ap) try candidate.addScalar(&candidate, 2) else candidate = try addWithOverflow(candidate, 2);

        if (try isMersenne(candidate, allocator)) return candidate;
    }
}

/// Get the next prime number after a given number.
pub fn nextPrime(p: anytype, allocator: std.mem.Allocator) !@TypeOf(p) {
    const T = @TypeOf(p); const ap = T == std.math.big.int.Managed;

    var two = if (comptime ap) try std.math.big.int.Managed.initSet(allocator, 2) else @as(T, 2); defer if (comptime ap) two.deinit();

    if (if (comptime ap) p.order(two) == std.math.Order.lt else p < 2) {
        return if (comptime ap) try std.math.big.int.Managed.initSet(allocator, 2) else 2;
    }

    var candidate = if (comptime ap) try p.clone() else if (p % 2 == 0) p - 1 else p;

    if (comptime ap) if (p.isEven()) try candidate.addScalar(&p, -1);

    while (true) {

        if (comptime ap) try candidate.addScalar(&candidate, 2) else candidate = try addWithOverflow(candidate, 2);

        if (try isPrime(candidate, allocator)) return candidate;
    }
}

/// Check if a number is prime using trial division.
pub fn trialDivision(p: anytype, allocator: std.mem.Allocator) !bool {
    const T = @TypeOf(p); const ap = T == std.math.big.int.Managed;

    var div = if (comptime ap) try std.math.big.int.Managed.initSet(allocator, 3) else @as(T, 3); defer if (comptime ap) div.deinit();
    var dsq = if (comptime ap) try std.math.big.int.Managed.initSet(allocator, 9) else @as(T, 9); defer if (comptime ap) dsq.deinit();
    var pfd = if (comptime ap) try std.math.big.int.Managed.initSet(allocator, 0) else @as(T, 0); defer if (comptime ap) pfd.deinit();
    var rem = if (comptime ap) try std.math.big.int.Managed.initSet(allocator, 0) else @as(T, 0); defer if (comptime ap) rem.deinit();
    var two = if (comptime ap) try std.math.big.int.Managed.initSet(allocator, 2) else @as(T, 2); defer if (comptime ap) two.deinit();

    if (if (comptime ap) p.order(two) == std.math.Order.lt else p < 2) return false;

    if (if (comptime ap) p.isEven() else p % 2 == 0) return if (comptime ap) p.eql(two) else p == 2;

    while (if (comptime ap) dsq.order(p) != std.math.Order.gt else dsq <= p) {

        if (comptime ap) try pfd.divFloor(&rem, &p, &div);

        if (if (comptime ap) rem.eqlZero() else p % div == 0) return false;

        if (comptime ap) try div.addScalar(&div, 2) else div += 2;

        if (comptime ap) try dsq.sqr(&div) else dsq = try squareWithOverflow(div);
    }

    return true;
}

/// Writes the array of prime numbers as a real matrix to a file path.
pub fn exportPrimeNumbersAsRealMatrix(path: []const u8, prime_numbers: anytype, allocator: std.mem.Allocator) !void {
    const ap = @TypeOf(prime_numbers[0]) == std.math.big.int.Managed;

    var file = try std.fs.cwd().createFile(path, .{}); defer file.close();

    var buffer: [WRITE_BUFFER_SIZE]u8 = undefined;

    var writer = file.writer(&buffer); var writer_interface = &writer.interface;

    try writer_interface.print("{d} {d}\n", .{prime_numbers.len, 2});

    for (0..prime_numbers.len) |i| {

        const prime_str = if (comptime ap) try prime_numbers[i].toString(allocator, 10, std.fmt.Case.lower) else "";

        try writer_interface.print("{d:20} " ++ if (comptime ap) "{s:20}" else "{d:20}" ++ "\n", .{i + 1, if (comptime ap) prime_str else prime_numbers[i]});

        if (comptime ap) allocator.free(prime_str);
    }

    try writer_interface.flush();
}
