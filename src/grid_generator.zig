//! Functions related to the generation of grids.

const std = @import("std");

const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

/// Get the momentum point given the row index.
pub fn momentumAtRow(comptime T: type, point: *RealVector(T), index: usize, ndim: usize, npoint: usize, limits: []const []const T) void {
    for (0..ndim) |i| {

        const axis = ndim - i - 1;

        const grid_index = @as(T, @floatFromInt(index / std.math.pow(usize, npoint, axis) % npoint));

        const shift = if (grid_index < @as(T, @floatFromInt(npoint / 2))) grid_index else grid_index - @as(T, @floatFromInt(npoint));

        point.ptr(i).* = 2 * std.math.pi * shift / @as(T, @floatFromInt(npoint)) / (limits[i][1] - limits[i][0]) * @as(T, @floatFromInt(npoint - 1));
    }
}

/// Generate a grid in the k-space. The result is allocated and returned. Both the start and end values are included.
pub fn momentumGridAlloc(comptime T: type, limits: []const []const T, npoint: usize, allocator: std.mem.Allocator) !RealMatrix(T) {
    var momentum_grid = try RealMatrix(T).init(std.math.pow(usize, npoint, limits.len), limits.len, allocator);

    for (0..momentum_grid.rows) |i| {

        var row = momentum_grid.row(i);

        momentumAtRow(T, &row, i, limits.len, npoint, limits);
    }

    return momentum_grid;
}

/// Get the position point given the row index.
pub fn positionAtRow(comptime T: type, point: *RealVector(T), index: usize, ndim: usize, npoint: usize, limits: []const []const T) void {
    for (0..ndim) |i| {

        const axis = ndim - i - 1;

        const grid_index = @as(T, @floatFromInt(index / std.math.pow(usize, npoint, axis) % npoint));

        point.ptr(i).* = limits[i][0] + grid_index * (limits[i][1] - limits[i][0]) / @as(T, @floatFromInt(npoint - 1));
    }
}

/// Generate a grid in the position space. The array is allocated inside the function.
pub fn positionGridAlloc(comptime T: type, limits: []const []const T, npoint: usize, allocator: std.mem.Allocator) !RealMatrix(T) {
    var position_grid = try RealMatrix(T).init(std.math.pow(usize, npoint, limits.len), limits.len, allocator);

    for (0..position_grid.rows) |i| {

        var row = position_grid.row(i);

        positionAtRow(T, &row, i, limits.len, npoint, limits);
    }

    return position_grid;
}
