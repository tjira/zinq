//! Some specific math functions.

const std = @import("std");

/// Binomial confidence interval given the probability and number of trials.
pub fn binomialConfInt(p: anytype, n: usize) @TypeOf(p) {
    return 1.96 * std.math.sqrt(p * (1 - p) / @as(@TypeOf(p), @floatFromInt(n)));
}

/// Calculate the combination number.
pub fn comb(n: anytype, k: @TypeOf(n)) @TypeOf(n, k) {
    var nck: @TypeOf(n) = 1;

    for (k + 1..n + 1) |i| nck *= i;
    for (2..n - k + 1) |i| nck /= i;

    return nck;
}

/// Generate all combinations of n elements from an array.
pub fn combinations(comptime T: type, array: []const T, n: usize, allocator: std.mem.Allocator) !std.ArrayList([]T) {
    var result = std.ArrayList([]T){}; var current = std.ArrayList(usize){}; defer current.deinit(allocator);

    const Backtrack = struct { fn get(res: *std.ArrayList([]T), arr: []const T, m: usize, curr: *std.ArrayList(usize), start: usize, alloc: std.mem.Allocator) !void {
        if (curr.items.len == m) {

            var last = try alloc.alloc(T, m);

            for (curr.items, 0..) |index, i| {
                last[i] = arr[index];
            }

            try res.append(alloc, last); return;
        }

        for (start..arr.len) |i| {
            try curr.append(alloc, i); try get(res, arr, m, curr, i + 1, alloc); _ = curr.pop();
        }
    }};

    try Backtrack.get(&result, array, n, &current, 0, allocator); return result;
}

/// Calculate the double factorial of a number.
pub fn dfact(n: anytype) @TypeOf(n) {
    if (n == -1 or n == 0 or n == 1) return 1;

    if (n == 2) return 2;

    return n * dfact(n - 2);
}

/// Heaviside step function.
pub fn h(x: anytype) @TypeOf(x) {
    if (x < 0) return 0; return 1;
}

/// Multiplication modulo function to prevent overflow.
pub fn mulMod(a: anytype, b: anytype, m: anytype) @TypeOf(a, b, m) {
    const T = @typeInfo(@TypeOf(a, b, m));

    return @intCast(@as(std.meta.Int(T.int.signedness, 2 * T.int.bits), a) * b % m);
}

/// Power function floats raised to integer exponents.
pub fn powi(base: anytype, exp: usize) @TypeOf(base) {
    var result: @TypeOf(base) = 1;

    for (0..exp) |_| {
        result *= base;
    }

    return result;
}

/// Exponentiation modulo function.
pub fn powMod(base: anytype, exp: anytype, m: anytype) @TypeOf(base, exp, m) {
    var result: @TypeOf(base, exp, m) = 1; var b = base % m; var e = exp;

    while (e > 0) {

        if (e % 2 == 1) result = mulMod(result, b, m);

        b = mulMod(b, b, m); e /= 2;
    }

    return result;
}

/// Reverse the bits of a number.
pub fn revk(value: anytype, k: u6) @TypeOf(value) {
    var result: @TypeOf(value) = 0; var i: u6 = 0;

    while (i < k) : (i += 1) {
        result |= ((value >> @intCast(i)) & 1) << @intCast(k - 1 - i);
    }

    return result;
}
