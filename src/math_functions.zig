//! Some specific math functions.

const std = @import("std");

/// Calculate the double factorial of a number.
pub fn dfact(n: anytype) @TypeOf(n) {
    if (n == -1 or n == 0 or n == 1) return 1;

    if (n == 2) return 2;

    return n * dfact(n - 2);
}

/// Reverse the bits of a number.
pub fn revk(value: anytype, k: u6) @TypeOf(value) {
    var result: @TypeOf(value) = 0; var i: u6 = 0;

    while (i < k) : (i += 1) {
        result |= ((value >> @intCast(i)) & 1) << @intCast(k - 1 - i);
    }

    return result;
}
