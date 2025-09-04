//! File with integral mathematical functions.

const std = @import("std");

const global_variables = @import("global_variables.zig");

const GAMMAINC_CUTOFF = global_variables.GAMMAINC_CUTOFF;

/// Boys function with zero n.
pub fn boys(comptime T: type, x: T, a: T) T {
    return if (x > 0) 0.5 * std.math.gamma(T, a + 0.5) * gammainc(T, x, a + 0.5) / std.math.pow(T, x, a + 0.5) else 1 / (2 * a + 1);
}

/// Regularized lower incomplete gamma function.
pub fn gammainc(comptime T: type, x: T, a: T) T {
    if (x < a + 1) {

        var b = 1 / a; var c = 1 / a; var i: T = 1;

        while (true) : (i += 1) {
            c *= x / (a + i); b += c; if (c < GAMMAINC_CUTOFF) break;
        }

        return b * std.math.exp(-x + a * std.math.log(T, std.math.e, x)) / std.math.gamma(T, a);
    }

    else {

        var b = x + 1 - a; var c: T = 1e10; var d = 1 / b; var e = d; var i: T = 1;

        while (true) : (i += 1) {
            const f = i * (a - i); b += 2; d = f * d + b; c = b + f / c; d = 1 / d; const g = d * c; e *= g; if (@abs(g - 1) < GAMMAINC_CUTOFF) break;
        }

        return 1 - std.math.exp(-x + a * std.math.log(T, std.math.e, x) - std.math.log(T, std.math.e, std.math.gamma(T, a))) * e;
    }
}
