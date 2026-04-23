//! Functions for linear interpolation.

const std = @import("std");

const global_variables = @import("global_variables.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

/// Function to get a function value of a point in a grid using multilinear interpolation.
pub fn lerp(comptime T: type, grid: RealMatrix(T), column: usize, r: RealVector(T)) !T {
    const max_dim = 16;

    if (r.len > max_dim) {
        std.log.err("LERP DIMENSION EXCEEDS MAXIMUM, MAXIMUM IS {d} BUT GOT {d}, THIS CAN BE CHANGED IN THE SOURCE CODE", .{ max_dim, r.len });

        return error.InvalidInput;
    }

    const size = @as(usize, @intFromFloat(@round(std.math.pow(T, @as(T, @floatFromInt(grid.rows)), 1 / @as(T, @floatFromInt(r.len))))));

    if (std.math.pow(usize, size, r.len) != grid.rows) {
        std.log.err("GRID SIZE INCONSISTENT WITH LERP DIMENSION, EXPECTED {d} BUT GOT {d}", .{ size, grid.rows });

        return error.InvalidInput;
    }

    var indices: [max_dim]usize = undefined;
    var weights: [max_dim]T = undefined;
    var strides: [max_dim]usize = undefined;

    for (0..r.len) |d| {
        strides[d] = std.math.pow(usize, size, r.len - d - 1);

        var low: usize = 0;
        var high: usize = size;
        var mid = (low + high) / 2;
        const target = r.at(d);

        while (low < high) : (mid = (low + high) / 2) {
            if (grid.at(mid * strides[d], d) <= target) low = mid + 1 else high = mid;
        }

        indices[d] = @min(@max(low, 1), size - 1);

        const x0 = grid.at((indices[d] - 1) * strides[d], d);
        const x1 = grid.at(indices[d] * strides[d], d);

        weights[d] = (target - x0) / (x1 - x0);
    }

    var result: T = 0;

    for (0..@as(usize, 1) << @as(u5, @intCast(r.len))) |k| {
        var weight: T = 1;
        var index: usize = 0;

        for (0..r.len) |d| {
            const use_upper = ((k >> @as(u5, @intCast(r.len - d - 1))) & 1) == 1;

            weight *= if (use_upper) weights[d] else 1.0 - weights[d];
            index += if (use_upper) indices[d] * strides[d] else (indices[d] - 1) * strides[d];
        }

        result += grid.at(index, column) * weight;
    }

    return result;
}
