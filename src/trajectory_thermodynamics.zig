//! File that contains various thermodynamic calculations for trajectories.

const std = @import("std");

const complex_matrix = @import("complex_matrix.zig");
const determinant = @import("determinant.zig");
const device_write = @import("device_write.zig");
const fourier_transform = @import("fourier_transform.zig");
const global_variables = @import("global_variables.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");
const system_alignment = @import("system_alignment.zig");

const ComplexMatrix = complex_matrix.ComplexMatrix;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const cfft1 = fourier_transform.cfft1;
const determinantHermitian = determinant.determinantHermitian;
const determinantAbs = determinant.determinantAbs;
const determinantLog = determinant.determinantLog;
const mmAlloc = matrix_multiplication.mmAlloc;
const alignTrajectory = system_alignment.alignTrajectory;

const AU2K = global_variables.AU2K;
const kB = global_variables.kB;
const Eh = global_variables.Eh;

/// Function to calculate the Schlitter entropy of a trajectory given the positions, masses, and temperature.
pub fn schlitterEntropy(comptime T: type, positions: RealMatrix(T), masses: RealVector(T), temp: T, ab_initio: bool, allocator: std.mem.Allocator) !T {
    const aligned = if (ab_initio) try alignTrajectory(T, positions, allocator) else positions;
    defer if (ab_initio) aligned.deinit(allocator);

    var M = try RealMatrix(T).initZero(masses.len, masses.len, allocator);
    defer M.deinit(allocator);

    for (0..masses.len) |i| M.ptr(i, i).* = masses.at(i);

    const sigma = try aligned.cov(allocator);
    defer sigma.deinit(allocator);

    var Ms = try mmAlloc(T, M, false, sigma, false, allocator);
    defer Ms.deinit(allocator);

    for (0..Ms.rows) |i| for (0..Ms.cols) |j| {
        Ms.ptr(i, j).* = if (i == j) 1 + temp / AU2K * std.math.e * std.math.e * Ms.at(i, j) else temp / AU2K * std.math.e * std.math.e * Ms.at(i, j);
    };

    const detlog = try determinantLog(T, Ms, allocator);

    return 0.5 * kB * detlog / Eh;
}

/// Function to calculate the QHA entropy of a trajectory given the positions, masses, and temperature.
pub fn qhaEntropy(comptime T: type, positions: RealMatrix(T), masses: RealVector(T), temp: T, ab_initio: bool, allocator: std.mem.Allocator) !T {
    const aligned = if (ab_initio) try alignTrajectory(T, positions, allocator) else positions;
    defer if (ab_initio) aligned.deinit(allocator);

    var M = try RealMatrix(T).initZero(masses.len, masses.len, allocator);
    defer M.deinit(allocator);

    for (0..masses.len) |i| M.ptr(i, i).* = masses.at(i);

    const sigma = try aligned.cov(allocator);
    defer sigma.deinit(allocator);

    var Ms = try mmAlloc(T, M, false, sigma, false, allocator);
    defer Ms.deinit(allocator);

    for (0..Ms.rows) |i| for (0..Ms.cols) |j| {
        Ms.ptr(i, j).* = if (i == j) temp / AU2K * std.math.e * std.math.e * Ms.at(i, j) else temp / AU2K * std.math.e * std.math.e * Ms.at(i, j);
    };

    const detlog = try determinantLog(T, Ms, allocator);

    return 0.5 * kB * detlog / Eh;
}
