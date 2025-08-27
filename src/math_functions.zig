//! Some specific math functions.

const std = @import("std");

/// Sign function.
pub fn sgn(x: anytype) @TypeOf(x) {
    return if (x > 0) 1 else if (x < 0) -1 else 0;
}
