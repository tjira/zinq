//! Some specific math functions.

const std = @import("std");

/// Binomial confidence interval given the probability and number of trials.
pub fn binomialConfInt(p: anytype, n: usize) @TypeOf(p) {
    return 1.96 * std.math.sqrt(p * (1 - p) / @as(@TypeOf(p), @floatFromInt(n)));
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
