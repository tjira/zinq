//! File with the prime number generation target as well as related functions.

const std = @import("std");

const device_write = @import("device_write.zig");

const print = device_write.print;

/// Options for the prime number generation target.
pub fn Options(comptime _: type) type {
    return struct {
        pub const Filter = enum {
            all
        };

        count: u32 = 100,
        log_interval: u32 = 10,
        start: u32 = 2,

        filter: Filter = .all,
    };
}

/// Output structure for the prime number generation target.
pub fn Output(comptime T: type) type {
    return struct {
        prime_numbers: []T,

        allocator: std.mem.Allocator,

        /// Initialize the output structure.
        pub fn init(count: T, allocator: std.mem.Allocator) !@This() {
            return @This(){
                .prime_numbers = try allocator.alloc(T, count),
                .allocator = allocator
            };
        }

        /// Free the output structure.
        pub fn deinit(self: @This()) void {
            self.allocator.free(self.prime_numbers);
        }
    };
}

/// Run the prime number generation target.
pub fn run(comptime T: type, options: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try print("\nGENERATING {d} PRIME NUMBERS STARTING FROM {d}\n\n", .{options.count, options.start});

    if (enable_printing) try print("{s:9} {s:18}\n", .{"INDEX", "PRIME NUMBER"});

    const output = try Output(T).init(@as(T, options.count), allocator);

    var p = @as(T, if (options.start < 2) 2 else options.start); var i: usize = 0;

    while (i < options.count) : (p = nextPrime(T, p)) {

        output.prime_numbers[i] = p; i += 1;

        if (enable_printing and (i == 1 or i % options.log_interval == 0)) {
            try print("{d:9} {d:18}\n", .{i, p});
        }
    }

    return output;
}

/// Get the next prime number after a given number.
pub fn nextPrime(comptime T: type, p: T) T {
    if (p < 2) return 2;

    if (p == 2) return 3;

    var candidate: T = if (p % 2 == 0) p + 1 else p + 2;

    while (true) : (candidate += 2) {
        if (trialDivision(T, candidate)) return candidate;
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
