//! File with functions specific to array of numbers.

const std = @import("std");

/// Sum of elements in an array.
pub fn sum(comptime T: type, arr: []const T) T {
    var total: T = 0;

    for (arr) |value| {
        total += value;
    }

    return total;
}

pub fn mean(comptime T: type, arr: []const T) T {
    var result: T = 0;
    const len = @as(T, @floatFromInt(arr.len));

    for (arr) |value| result += value / len;

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

/// Calculate sample standard deviation of an array.
pub fn sd(comptime T: type, arr: []const T) T {
    const len = @as(T, @floatFromInt(arr.len));
    const mean_value = mean(T, arr);

    var variance: T = 0;

    for (arr) |value| {
        const diff = value - mean_value;
        variance += diff * diff;
    }

    return std.math.sqrt(variance / (len - 1));
}
