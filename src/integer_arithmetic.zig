//! This module provides integer arithmetic operations with overflow checking as well as some other utilities.

/// Adds two numbers and checks for overflow.
pub fn addWithOverflow(x: anytype, y: anytype) !@TypeOf(x, y) {
    const ov = @addWithOverflow(x, y);

    if (ov[1] != 0) return error.Overflow;

    return ov[0];
}

/// Squares a number and checks for overflow.
pub fn squareWithOverflow(x: anytype) !@TypeOf(x) {
    const ov = @mulWithOverflow(x, x);

    if (ov[1] != 0) return error.Overflow;

    return ov[0];
}
