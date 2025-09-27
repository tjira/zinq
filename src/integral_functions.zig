//! File with integral mathematical functions.

const std = @import("std");

const global_variables = @import("global_variables.zig");

const BOYS_CUTOFF = global_variables.BOYS_CUTOFF;

/// Boys function.
pub fn boys(n: usize, t: anytype) @TypeOf(t) {
    if (t < BOYS_CUTOFF) {
        return 1 / @as(@TypeOf(t), @floatFromInt(2 * n + 1)) - t / @as(@TypeOf(t), @floatFromInt(2 * n + 3));
    } else if (n == 0) {
        return std.math.sqrt(std.math.pi / (4 * t)) * erf(std.math.sqrt(t));
    } else {
        return (@as(@TypeOf(t), @floatFromInt(n)) - 0.5) * boys(n - 1, t) / t - std.math.exp(-t) / (2 * t);
    }
}

/// Error function. Using the Cody's rational Chebyshev approximations.
pub fn erf(x: anytype) @TypeOf(x) {
    if (x < 0){
        return -erf(-x);
    }

    else if (x < 0.46875) {

        const p = [_]@TypeOf(x){
            3.209377589138469472562e+3,
            3.774852376853020208137e+2,
            1.138641541510501556495e+2,
            3.161123743870565596947e+0,
            1.857777061846031526730e-1
        };

        const q = [_]@TypeOf(x){
            2.844236833439170622273e+3,
            1.282616526077372275645e+3,
            2.440246379344441733056e+2,
            2.360129095234412093499e+1,
            1.000000000000000000000e+0
        };

        const t = x * x;

        const P = p[0] + t * (p[1] + t * (p[2] + t * (p[3] + t * p[4])));
        const Q = q[0] + t * (q[1] + t * (q[2] + t * (q[3] + t * q[4])));

        return x * P / Q;
    }

    else if (x < 4) {

        const p = [_]@TypeOf(x){
            1.23033935479799725272e+3,
            2.05107837782607146532e+3,
            1.71204761263407058314e+3,
            8.81952221241769090411e+2,
            2.98635138197400131132e+2,
            6.61191906371416294775e+1,
            8.88314979438837594118e+0,
            5.64188496988670089180e-1,
            2.15311535474403846343e-8
        };

        const q = [_]@TypeOf(x){
            1.23033935480374942043e+3,
            3.43936767414372163696e+3,
            4.36261909014324715820e+3,
            3.29079923573345962678e+3,
            1.62138957456669018874e+3,
            5.37181101862009857509e+2,
            1.17693950891312499305e+2,
            1.57449261107098347253e+1,
            1.00000000000000000000e+0
        };

        const t = x;

        const P = p[0] + t * (p[1] + t * (p[2] + t * (p[3] + t * (p[4] + t * (p[5] + t * (p[6] + t * (p[7] + t * p[8])))))));
        const Q = q[0] + t * (q[1] + t * (q[2] + t * (q[3] + t * (q[4] + t * (q[5] + t * (q[6] + t * (q[7] + t * q[8])))))));

        return 1 - std.math.exp(-x * x) * P / Q;
    }

    else {

        const p = [_]@TypeOf(x){
            -6.58749161529837803157e-4,
            -1.60837851487422766278e-2,
            -1.25781726111229246204e-1,
            -3.60344899949804439429e-1,
            -3.05326634961232344035e-1,
            -1.63153871373020978498e-2
        };

        const q = [_]@TypeOf(x){
            2.33520497626869185443e-3,
            6.05183413124413191178e-2,
            5.27905102951428412248e-1,
            1.87295284992346047209e+0,
            2.56852019228982242072e+0,
            1.00000000000000000000e+0
        };

        const t = 1.0 / (x * x);

        const P = p[0] + t * (p[1] + t * (p[2] + t * (p[3] + t * (p[4] + t * p[5]))));
        const Q = q[0] + t * (q[1] + t * (q[2] + t * (q[3] + t * (q[4] + t * q[5]))));

        return 1 - (std.math.exp(-x * x) / x * (1.0 / std.math.sqrt(std.math.pi) + P / Q / (x * x)));
    }
}
