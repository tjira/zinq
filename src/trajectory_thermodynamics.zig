//! File that contains various thermodynamic calculations for trajectories.

const std = @import("std");

const determinant = @import("determinant.zig");
const device_write = @import("device_write.zig");
const global_variables = @import("global_variables.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const determinantHermitian = determinant.determinantHermitian;
const mmAlloc = matrix_multiplication.mmAlloc;

const AU2K = global_variables.AU2K;
const kB = global_variables.kB;
const Eh = global_variables.Eh;

pub fn schlitterEntropy(comptime T: type, positions: RealMatrix(T), masses: RealVector(T), temp: T, allocator: std.mem.Allocator) !T {
    var M = try RealMatrix(T).initZero(masses.len, masses.len, allocator); defer M.deinit(allocator);

    for (0..masses.len) |i| M.ptr(i, i).* = masses.at(i);

    const sigma = try positions.cov(allocator); defer sigma.deinit(allocator);

    var Ms = try mmAlloc(T, M, false, sigma, false, allocator); defer Ms.deinit(allocator);

    for (0..Ms.rows) |i| for (0..Ms.cols) |j| {
        Ms.ptr(i, j).* = if (i == j) 1 + temp / AU2K * std.math.e * std.math.e * Ms.at(i, j) else temp / AU2K * std.math.e * std.math.e * Ms.at(i, j);
    };

    const det = try determinantHermitian(T, Ms, allocator);

    return 0.5 * kB * std.math.log(T, std.math.e, det) / Eh;
}
