//! Some specific math functions.

const std = @import("std");

/// Reverse the bits of a number.
pub fn reverse_bits(value: anytype, count: u6) @TypeOf(value) {
    var result: @TypeOf(value) = 0; var i: u6 = 0;

    while (i < count) : (i += 1) {
        result |= ((value >> @intCast(i)) & 1) << @intCast(count - 1 - i);
    }

    return result;
}

/// Product of elements in an array.
pub fn prod(comptime T: type, arr: []const T) T {
    var product: T = 1;

    for (arr) |value| {
        product *= value;
    }

    return product;
}

/// Sign function.
pub fn sgn(x: anytype) @TypeOf(x) {
    return if (x > 0) 1 else if (x < 0) -1 else 0;
}
