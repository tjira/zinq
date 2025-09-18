//! Functions for linear interpolation.

const std = @import("std");

const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

/// Function to get a function value of a point in a grid using linear interpolation.
pub fn lerp(comptime T: type, grid: RealMatrix(T), column: usize, r: RealVector(T)) !T {
    if (r.len != 1) return error.InvalidDimension;

    var i: usize = 1;

    while (i < grid.rows and r.at(0) > grid.at(i, r.len - 1)) : (i += 1) {}

    if (i >= grid.rows) i = grid.rows - 1;

    const x0 = grid.at(i - 1, 0);
    const x1 = grid.at(i,     0);

    const y0 = grid.at(i - 1, column);
    const y1 = grid.at(i,     column);

    const t = (r.at(0) - x0) / (x1 - x0);

    return std.math.lerp(y0, y1, t);
}
