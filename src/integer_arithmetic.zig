//! This module provides integer arithmetic operations with overflow checking as well as some other utilities.

const std = @import("std");

const error_handling = @import("error_handling.zig");

const throw = error_handling.throw;

/// Adds two numbers and checks for overflow.
pub fn addWithOverflow(x: anytype, y: anytype) !@TypeOf(x, y) {
    const ov = @addWithOverflow(x, y);

    if (ov[1] != 0) return throw(@TypeOf(x), "OVERFLOW WHEN ADDING {d} AND {d} FOR {s} TYPE", .{x, y, @typeName(@TypeOf(x))});

    return ov[0];
}

/// Squares a number and checks for overflow.
pub fn squareWithOverflow(x: anytype) !@TypeOf(x) {
    const ov = @mulWithOverflow(x, x);

    if (ov[1] != 0) return throw(@TypeOf(x), "OVERFLOW WHEN SQUARING {d} FOR {s} TYPE", .{x, @typeName(@TypeOf(x))});

    return ov[0];
}
