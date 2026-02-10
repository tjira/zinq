//! Functions for linear interpolation.

const std = @import("std");

const error_handling = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const throw = error_handling.throw;

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const MAX_LERP_DIM = global_variables.MAX_LERP_DIM;

/// Function to get a function value of a point in a grid using multilinear interpolation.
pub fn lerp(comptime T: type, grid: RealMatrix(T), column: usize, r: RealVector(T)) !T {
    if (r.len > MAX_LERP_DIM) return throw(T, "LINEAR INTERPOLATION ONLY IMPLEMENTED FOR UP TO {d} DIMENSIONS", .{MAX_LERP_DIM});

    const size = @as(usize, @intFromFloat(@round(std.math.pow(T, @as(T, @floatFromInt(grid.rows)), 1 / @as(T, @floatFromInt(r.len))))));

    if (std.math.pow(usize, size, r.len) != grid.rows) return throw(T, "DATA PASSED TO LERP FUNCTION DOES NOT HAVE THE SAME NUMBER OF POINTS IN EACH DIMENSION", .{});

    var indices: [MAX_LERP_DIM]usize = undefined;
    var weights: [MAX_LERP_DIM]T = undefined;
    var strides: [MAX_LERP_DIM]usize = undefined;

    for (0..r.len) |d| {

        strides[d] = std.math.pow(usize, size, r.len - d - 1);

        var low: usize = 0; var high: usize = size; var mid = (low + high) / 2; const target = r.at(d);

        while (low < high) : (mid = (low + high) / 2) {
            if (grid.at(mid * strides[d], d) <= target) low = mid + 1 else high = mid;
        }
        
        indices[d] = @min(@max(low, 1), size - 1);

        const x0 = grid.at((indices[d] - 1) * strides[d], d);
        const x1 = grid.at(indices[d] * strides[d],       d);

        weights[d] = (target - x0) / (x1 - x0);
    }

    var result: T = 0;

    for (0..@as(usize, 1) << @as(u5, @intCast(r.len))) |k| {

        var weight: T = 1; var index: usize = 0;

        for (0..r.len) |d| {

            const use_upper = ((k >> @as(u5, @intCast(r.len - d - 1))) & 1) == 1;

            weight *= if (use_upper) weights[d] else 1.0 - weights[d];
            index += if (use_upper) indices[d] * strides[d] else (indices[d] - 1) * strides[d];
        }

        result += grid.at(index, column) * weight;
    }

    return result;
}
