//! This module provides integer arithmetic operations with overflow checking as well as some other utilities.

const std = @import("std");

const error_handling = @import("error_handling.zig");
const fourier_transform = @import("fourier_transform.zig");
const math_functions = @import("math_functions.zig");

const mulMod = math_functions.mulMod;
const ntt = fourier_transform.ntt;
const throw = error_handling.throw;

/// Adds two numbers and checks for overflow.
pub fn addWithOverflow(x: anytype, y: anytype) !@TypeOf(x, y) {
    const ov = @addWithOverflow(x, y);

    if (ov[1] != 0) return throw(@TypeOf(x), "OVERFLOW WHEN ADDING {d} AND {d} FOR {s} TYPE", .{x, y, @typeName(@TypeOf(x))});

    return ov[0];
}

/// Strassen squaring of a number.
pub fn strassenSquare(result: *std.math.big.int.Managed, x: *const std.math.big.int.Managed, allocator: std.mem.Allocator) !void {
    var size: usize = 1; while (size < 2 * x.len() * @bitSizeOf(usize)) size <<= 1;

    var bits = try allocator.alloc(usize, size); defer allocator.free(bits); @memset(bits, 0);

    for (0..x.len()) |i| inline for (0..@bitSizeOf(usize)) |j| {
        bits[i * @bitSizeOf(usize) + j] = @truncate(x.limbs[i] >> @intCast(j) & 1);
    };

    const mod: @TypeOf(bits[0]) = 998244353; const g: @TypeOf(bits[0]) = 3;

    ntt(bits, mod, g, false);

    for (0..bits.len) |i| {
        bits[i] = mulMod(bits[i], bits[i], mod);
    }

    ntt(bits, mod, g, true);

    try result.set(0);

    for (0..bits.len) |i| if (bits[i] != 0) {

        var carry = try std.math.big.int.Managed.initSet(allocator, bits[i]); defer carry.deinit();

        var shifted = try std.math.big.int.Managed.init(allocator); defer shifted.deinit();

        try shifted.shiftLeft(&carry, i);

        try result.add(result, &shifted);
    };
}

/// Squares a number and checks for overflow.
pub fn squareWithOverflow(x: anytype) !@TypeOf(x) {
    const ov = @mulWithOverflow(x, x);

    if (ov[1] != 0) return throw(@TypeOf(x), "OVERFLOW WHEN SQUARING {d} FOR {s} TYPE", .{x, @typeName(@TypeOf(x))});

    return ov[0];
}
