//! Functions for linear interpolation.

const std = @import("std");

const error_handling = @import("error_handling.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const throw = error_handling.throw;

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

/// Function to get a function value of a point in a grid using multilinear interpolation.
pub fn lerp(comptime T: type, grid: RealMatrix(T), column: usize, r: RealVector(T)) !T {
    switch (r.len) {
        1 => return lerp1D(T, grid, column, r),
        2 => return lerp2D(T, grid, column, r),
        3 => return lerp3D(T, grid, column, r),
        else => return throw(T, "LINEAR INTERPOLATION ONLY IMPLEMENTED FOR 1D, 2D, AND 3D DATA", .{})
    }
}

/// Function to get a function value of a point in a grid using 1D linear interpolation.
pub fn lerp1D(comptime T: type, grid: RealMatrix(T), column: usize, r: RealVector(T)) T {
    var i = grid.column(r.len - 1).bisectRight(r.at(0));

    if (i >= grid.rows) i = grid.rows - 1;

    const x0 = grid.at(i - 1, r.len - 1);
    const x1 = grid.at(i,     r.len - 1);

    const t0 = (r.at(0) - x0) / (x1 - x0);

    const f0 = grid.at(i - 1, column);
    const f1 = grid.at(i,     column);

    return std.math.lerp(f0, f1, t0);
}

/// Function to get a function value of a point in a grid using 2D linear interpolation.
pub fn lerp2D(comptime T: type, grid: RealMatrix(T), column: usize, r: RealVector(T)) T {
    const size: usize = @as(usize, @intFromFloat(@round(std.math.pow(T, @as(T, @floatFromInt(grid.rows)), 1.0 / 2.0))));

    var i = grid.column(r.len - 1).slice(0, size).bisectRight(r.at(0));
    var j = grid.column(r.len - 1).slice(0, size).bisectRight(r.at(1));

    if (i >= grid.rows) i = grid.rows - 1;
    if (j >= grid.rows) j = grid.rows - 1;

    const x0 = grid.at(i - 1, r.len - 1);
    const x1 = grid.at(i,     r.len - 1);

    const y0 = grid.at(j - 1, r.len - 1);
    const y1 = grid.at(j,     r.len - 1);

    const t0 = (r.at(0) - x0) / (x1 - x0);
    const t1 = (r.at(1) - y0) / (y1 - y0);

    const f00 = grid.at((i - 1) * size + j - 1, column);
    const f01 = grid.at((i - 1) * size + j,     column);
    const f10 = grid.at(i       * size + j - 1, column);
    const f11 = grid.at(i       * size + j,     column);

    const f0 = std.math.lerp(f00, f10, t0);
    const f1 = std.math.lerp(f01, f11, t0);

    return std.math.lerp(f0, f1, t1);
}

/// Function to get a function value of a point in a grid using 3D linear interpolation.
pub fn lerp3D(comptime T: type, grid: RealMatrix(T), column: usize, r: RealVector(T)) T {
    const size: usize = @as(usize, @intFromFloat(@round(std.math.pow(T, @as(T, @floatFromInt(grid.rows)), 1.0 / 3.0))));

    var i = grid.column(r.len - 1).slice(0, size).bisectRight(r.at(0));
    var j = grid.column(r.len - 1).slice(0, size).bisectRight(r.at(1));
    var k = grid.column(r.len - 1).slice(0, size).bisectRight(r.at(2));

    if (i >= grid.rows) i = grid.rows - 1;
    if (j >= grid.rows) j = grid.rows - 1;
    if (k >= grid.rows) k = grid.rows - 1;

    const x0 = grid.at(i - 1, r.len - 1);
    const x1 = grid.at(i,     r.len - 1);

    const y0 = grid.at(j - 1, r.len - 1);
    const y1 = grid.at(j,     r.len - 1);

    const z0 = grid.at(k - 1, r.len - 1);
    const z1 = grid.at(k,     r.len - 1);

    const t0 = (r.at(0) - x0) / (x1 - x0);
    const t1 = (r.at(1) - y0) / (y1 - y0);
    const t2 = (r.at(2) - z0) / (z1 - z0);

    const f000 = grid.at((i - 1) * size * size + (j - 1) * size + k - 1, column);
    const f001 = grid.at((i - 1) * size * size + (j - 1) * size + k,     column);
    const f010 = grid.at((i - 1) * size * size + j       * size + k - 1, column);
    const f011 = grid.at((i - 1) * size * size + j       * size + k,     column);
    const f100 = grid.at(i       * size * size + (j - 1) * size + k - 1, column);
    const f101 = grid.at(i       * size * size + (j - 1) * size + k,     column);
    const f110 = grid.at(i       * size * size + j       * size + k - 1, column);
    const f111 = grid.at(i       * size * size + j       * size + k,     column);

    const f00 = std.math.lerp(f000, f100, t0);
    const f01 = std.math.lerp(f001, f101, t0);
    const f10 = std.math.lerp(f010, f110, t0);
    const f11 = std.math.lerp(f011, f111, t0);

    const f0 = std.math.lerp(f00, f10, t1);
    const f1 = std.math.lerp(f01, f11, t1);

    return std.math.lerp(f0, f1, t2);
}
