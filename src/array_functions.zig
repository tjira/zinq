//! File with functions specific to array of numbers.

/// Sum of elements in an array.
pub fn sum(comptime T: type, arr: []const T) T {
    var total: T = 0;

    for (arr) |value| {
        total += value;
    }

    return total;
}

/// Product of elements in an array.
pub fn prod(comptime T: type, arr: []const T) T {
    var product: T = 1;

    for (arr) |value| {
        product *= value;
    }

    return product;
}
