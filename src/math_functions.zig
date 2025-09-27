//! Some specific math functions.

const std = @import("std");

/// Calculate the double factorial of a number.
pub fn dfact(n: anytype) @TypeOf(n) {
    if (n == -1 or n == 0 or n == 1) return 1;

    if (n == 2) return 2;

    return n * dfact(n - 2);
}

/// Power function floats raised to integer exponents.
pub fn powi(base: anytype, exp: usize) @TypeOf(base) {
    var result: @TypeOf(base) = 1;

    for (0..exp) |_| {
        result *= base;
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
