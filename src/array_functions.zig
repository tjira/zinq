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
    if (arr.len == 0) return 0;
    return sum(T, arr) / @as(T, @floatFromInt(arr.len));
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

test "array functions" {
    const testing = std.testing;

    const data = [_]f64{ 1.0, 2.0, 3.0, 4.0, 5.0 };

    try testing.expectEqual(@as(f64, 15.0), sum(f64, &data));
    try testing.expectEqual(@as(f64, 3.0), mean(f64, &data));
    try testing.expectEqual(@as(f64, 120.0), prod(f64, &data));

    try testing.expectApproxEqAbs(@as(f64, 1.5811388300841898), sd(f64, &data), 1e-10);
}
