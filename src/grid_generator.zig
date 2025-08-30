//! Functions related to the generation of grids.

const std = @import("std");

const real_matrix = @import("real_matrix.zig");

const RealMatrix = real_matrix.RealMatrix;

/// Generate a grid in the k-space. The result is allocated and returned. Both the start and end values are included.
pub fn momentumGridAlloc(comptime T: type, limits: []const []const T, npoint: usize, allocator: std.mem.Allocator) !RealMatrix(T) {
    var momentum_grid = try RealMatrix(T).init(std.math.pow(usize, npoint, limits.len), limits.len, allocator);

    for (0..momentum_grid.rows) |i| for (0..momentum_grid.cols) |j| {

        const axis = momentum_grid.cols - j - 1;

        const grid_index = @as(T, @floatFromInt(i / std.math.pow(usize, npoint, axis) % npoint));

        const shift = if (grid_index < @as(T, @floatFromInt(npoint / 2))) grid_index else grid_index - @as(T, @floatFromInt(npoint));

        momentum_grid.ptr(i, j).* = 2 * std.math.pi * shift / @as(T, @floatFromInt(npoint)) / (limits[j][1] - limits[j][0]) * @as(T, @floatFromInt(npoint - 1));
    };

    return momentum_grid;
}

/// Generate a grid in the position space. The array is allocated inside the function.
pub fn positionGridAlloc(comptime T: type, limits: []const []const T, npoint: usize, allocator: std.mem.Allocator) !RealMatrix(T) {
    var position_grid = try RealMatrix(T).init(std.math.pow(usize, npoint, limits.len), limits.len, allocator);

    for (0..position_grid.rows) |i| for (0..position_grid.cols) |j| {

        const axis = position_grid.cols - j - 1;

        const grid_index = @as(T, @floatFromInt(i / std.math.pow(usize, npoint, axis) % npoint));

        position_grid.ptr(i, j).* = limits[j][0] + grid_index * (limits[j][1] - limits[j][0]) / @as(T, @floatFromInt(npoint - 1));
    };

    return position_grid;
}
