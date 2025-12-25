//! File with the prime number generation target as well as related functions.

const std = @import("std");

const device_write = @import("device_write.zig");
const global_variables = @import("global_variables.zig");
const real_vector = @import("real_vector.zig");

const RealVector = real_vector.RealVector;

const exportRealMatrix = device_write.exportRealMatrix;
const print = device_write.print;
const printJson = device_write.printJson;
const printRealMatrix = device_write.printRealMatrix;

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

        mode: Mode
    };
}

/// Output structure for the prime number generation target.
pub fn Output(comptime _: type) type {
    return struct {
        prime_numbers: []std.math.big.int.Managed,

        /// Initialize the output structure.
        pub fn init(count: usize, allocator: std.mem.Allocator) !@This() {
            return @This(){
                .prime_numbers = try allocator.alloc(std.math.big.int.Managed, count),
            };
        }

        /// Free the output structure.
        pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
            for (self.prime_numbers) |*p| p.deinit();

            allocator.free(self.prime_numbers);
        }
    };
}

/// Run the prime number generation target.
pub fn run(comptime T: type, opt: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(opt);

    switch (opt.mode) {
        .check => return try checkMode(T, opt.mode.check, enable_printing, allocator),
        .factorize => return try factorizeMode(T, opt.mode.factorize, enable_printing, allocator),
        .generate => return try generateMode(T, opt.mode.generate, enable_printing, allocator),
    }
}

/// Check mode.
pub fn checkMode(comptime T: type, opt: anytype, enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    var number = try std.math.big.int.Managed.init(allocator); defer number.deinit();

    try number.setString(10, opt.number);

    const is_prime = switch (opt.filter) {
        .all => try isPrime(number, allocator),
        .mersenne => try isMersenne(try number.toInt(u64), allocator),
    };

    const output = try Output(T).init(if (is_prime) 1 else 0, allocator);

    if (is_prime) {
        output.prime_numbers[0] = try number.clone();
    }

    if (enable_printing) {

        const s = if (opt.filter == .mersenne) " MERSENNE" else "";

        try print("\n{s} IS{s} A{s} PRIME NUMBER\n", .{opt.number, if (is_prime) "" else " NOT", s});
    }

    return output;
}

/// Factorize mode.
pub fn factorizeMode(comptime T: type, opt: anytype, enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    var number = try std.math.big.int.Managed.init(allocator); defer number.deinit();

    try number.setString(10, opt.number);

    const factors = try factorize(number, allocator);

    const output = Output(T){.prime_numbers = factors};

    if (enable_printing) {

        if (factors.len == 0) {
            try print("\n{s} HAS NO PRIME FACTORS\n", .{opt.number}); return output;
        }

        try print("\n", .{});

        for (factors, 1..) |factor, i| {

            const factor_str = try factor.toString(allocator, 10, std.fmt.Case.lower); defer allocator.free(factor_str);

            try print("FACTOR {d:4}: {s}\n", .{i, factor_str});
        }
    }

    return output;
}

/// Generate mode.
pub fn generateMode(comptime T: type, opt: anytype, enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try print("\n{s:9} {s:18}\n", .{"INDEX", "PRIME NUMBER"});

    var start = try std.math.big.int.Managed.init(allocator); defer start.deinit();

    try start.setString(10, opt.start);

    var two = try std.math.big.int.Managed.initSet(allocator, 2); defer two.deinit();

    const max_output = if (opt.output != null and opt.output.?.interval != null) @min(opt.count, opt.output.?.interval.?) else opt.count;

    var output = try Output(T).init(max_output, allocator);

    output.prime_numbers[0] = if (start.order(two) != std.math.Order.gt) try two.clone() else switch (opt.filter) {
        .all => try nextPrime(start, allocator),
        .mersenne => try std.math.big.int.Managed.initSet(allocator, try nextMersenne(try start.toInt(u64), allocator))
    };

    for (0..opt.count) |i| {

        if (i >= max_output) output.prime_numbers[i % max_output].deinit();

        if (i > 0 or output.prime_numbers.len == 0) output.prime_numbers[i % max_output] = switch (opt.filter) {
            .all => try nextPrime(output.prime_numbers[(i - 1) % max_output], allocator),
            .mersenne => try std.math.big.int.Managed.initSet(allocator, try nextMersenne(try output.prime_numbers[(i - 1) % max_output].toInt(u64), allocator))
        };

        if (enable_printing and (i == 0 or (i + 1) % opt.log_interval == 0)) {

            const prime_string = try output.prime_numbers[i % max_output].toString(allocator, 10, std.fmt.Case.lower); defer allocator.free(prime_string);

            try print("{d:9} {s:18}\n", .{i + 1, prime_string});
        }

        if ((i + 1) % max_output == 0) if (opt.output != null) {

            var buf_start: [32]u8 = undefined; var buf_end: [32]u8 = undefined;

            const start_i_str = try std.fmt.bufPrint(&buf_start, "{d}", .{i + 2 - max_output}); const end_i_str = try std.fmt.bufPrint(&buf_end, "{d}", .{i + 1});

            const path = if (opt.output.?.interval != null) try std.mem.concat(allocator, u8, &.{opt.output.?.path, "_", start_i_str, "-", end_i_str}) else opt.output.?.path;

            try exportPrimeNumbersAsRealMatrix(path, output.prime_numbers, allocator);

            if (opt.output.?.interval != null) allocator.free(path);
        };
    }

    return output;
}

/// Factorize a number into its prime factors.
pub fn factorize(p: std.math.big.int.Managed, allocator: std.mem.Allocator) ![]std.math.big.int.Managed {
    var one = try std.math.big.int.Managed.initSet(allocator, 1); defer one.deinit();
    var two = try std.math.big.int.Managed.initSet(allocator, 2); defer two.deinit();
    var rem = try std.math.big.int.Managed.initSet(allocator, 0); defer rem.deinit();
    var pfd = try std.math.big.int.Managed.initSet(allocator, 0); defer pfd.deinit();

    if (p.order(two) == std.math.Order.lt) return &[_]std.math.big.int.Managed{};

    var factors: std.ArrayList(std.math.big.int.Managed) = .empty;

    var n = try p.clone(); var pd = try std.math.big.int.Managed.initSet(allocator, 0);

    while (!n.eql(one)) {

        var pd_next = try nextPrime(pd, allocator); defer pd_next.deinit(); try pd.copy(pd_next.toConst());

        try std.math.big.int.Managed.divFloor(&pfd, &rem, &n, &pd);

        while (rem.eqlZero()) {

            try factors.append(allocator, try pd.clone()); try n.copy(pfd.toConst());

            try std.math.big.int.Managed.divFloor(&pfd, &rem, &n, &pd);

        }
    }

    n.deinit(); pd.deinit();

    return factors.toOwnedSlice(allocator);
}

/// Check if the number is prime.
pub fn isPrime(p: std.math.big.int.Managed, allocator: std.mem.Allocator) !bool {
    return trialDivision(p, allocator);
}

/// Check if the number is Mersenne prime.
pub fn isMersenne(p: anytype, allocator: std.mem.Allocator) !bool {
    return try lucasLehmer(p, allocator);
}

/// The Lucas-Lehmer test for Mersenne primes.
pub fn lucasLehmer(p: anytype, allocator: std.mem.Allocator) !bool {
    var p_big = try std.math.big.int.Managed.initSet(allocator, p); defer p_big.deinit();

    if (!try isPrime(p_big, allocator)) return false;

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
pub fn nextMersenne(p: anytype, allocator: std.mem.Allocator) !@TypeOf(p) {
    if (p < 2) return 2;

    if (p == 2) return 3;

    var candidate = if (p % 2 == 0) p + 1 else p + 2;

    while (true) : (candidate += 2) {
        if (try isMersenne(candidate, allocator)) return candidate;
    }
}

/// Get the next prime number after a given number.
pub fn nextPrime(p: std.math.big.int.Managed, allocator: std.mem.Allocator) !std.math.big.int.Managed {
    var pfd = try std.math.big.int.Managed.initSet(allocator, 0); defer pfd.deinit();
    var rem = try std.math.big.int.Managed.initSet(allocator, 0); defer rem.deinit();
    var two = try std.math.big.int.Managed.initSet(allocator, 2); defer two.deinit();

    if (p.order(two) == std.math.Order.lt) {
        return try std.math.big.int.Managed.initSet(allocator, 2);
    }

    if (p.eql(two)) {
        return try std.math.big.int.Managed.initSet(allocator, 3);
    }

    var candidate = try p.clone();

    try std.math.big.int.Managed.addScalar(&candidate, &candidate, @as(u64, if (p.isEven()) 1 else 2));

    while (true) : (try std.math.big.int.Managed.addScalar(&candidate, &candidate, 2)) {
        if (try isPrime(candidate, allocator)) return candidate;
    }
}

/// Check if a number is prime using trial division.
pub fn trialDivision(p: std.math.big.int.Managed, allocator: std.mem.Allocator) !bool {
    var div = try std.math.big.int.Managed.initSet(allocator, 3); defer div.deinit();
    var dsq = try std.math.big.int.Managed.initSet(allocator, 9); defer dsq.deinit();
    var pfd = try std.math.big.int.Managed.initSet(allocator, 0); defer pfd.deinit();
    var rem = try std.math.big.int.Managed.initSet(allocator, 0); defer rem.deinit();
    var two = try std.math.big.int.Managed.initSet(allocator, 2); defer two.deinit();

    if (p.order(two) == std.math.Order.lt) return false;

    if (p.eql(two)) return true;

    if (p.isEven()) return false;

    while (dsq.order(p) != std.math.Order.gt) {

        try std.math.big.int.Managed.divFloor(&pfd, &rem, &p, &div);

        if (rem.eqlZero()) return false;

        try std.math.big.int.Managed.addScalar(&div, &div, 2);

        try std.math.big.int.Managed.sqr(&dsq, &div);
    }

    return true;
}

/// Writes the array of prime numbers as a real matrix to a file path.
pub fn exportPrimeNumbersAsRealMatrix(path: []const u8, prime_numbers: []std.math.big.int.Managed, allocator: std.mem.Allocator) !void {
    var file = try std.fs.cwd().createFile(path, .{}); defer file.close();

    var buffer: [WRITE_BUFFER_SIZE]u8 = undefined;

    var writer = file.writer(&buffer); var writer_interface = &writer.interface;

    try writer_interface.print("{d} {d}\n", .{prime_numbers.len, 2});

    for (0..prime_numbers.len) |i| {

        const prime_str = try prime_numbers[i].toString(allocator, 10, std.fmt.Case.lower); defer allocator.free(prime_str);

        try writer_interface.print("{d:20} {s:20}\n", .{i + 1, prime_str});
    }

    try writer_interface.flush();
}
