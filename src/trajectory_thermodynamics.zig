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
const mmAlloc = matrix_multiplication.mmAlloc;
const alignTrajectory = system_alignment.alignTrajectory;

const AU2K = global_variables.AU2K;
const kB = global_variables.kB;
const Eh = global_variables.Eh;

/// Function to calculate the Schlitter entropy of a trajectory given the positions, masses, and temperature.
pub fn schlitterEntropy(comptime T: type, positions: RealMatrix(T), masses: RealVector(T), temp: T, ab_initio: bool, allocator: std.mem.Allocator) !T {
    const aligned = if (ab_initio) try alignTrajectory(T, positions, allocator) else positions; defer if (ab_initio) aligned.deinit(allocator);

    var M = try RealMatrix(T).initZero(masses.len, masses.len, allocator); defer M.deinit(allocator);

    for (0..masses.len) |i| M.ptr(i, i).* = masses.at(i);

    const sigma = try aligned.cov(allocator); defer sigma.deinit(allocator);

    var Ms = try mmAlloc(T, M, false, sigma, false, allocator); defer Ms.deinit(allocator);

    for (0..Ms.rows) |i| for (0..Ms.cols) |j| {
        Ms.ptr(i, j).* = if (i == j) 1 + temp / AU2K * std.math.e * std.math.e * Ms.at(i, j) else temp / AU2K * std.math.e * std.math.e * Ms.at(i, j);
    };

    const det = try determinantAbs(T, Ms, allocator);

    return 0.5 * kB * std.math.log(T, std.math.e, det) / Eh;
}

/// Function to calculate the Spectral Resolved Estimate of the entropy of a trajectory given the momenta, masses, and temperature.
pub fn sreEntropy(comptime T: type, positions: RealMatrix(T), masses: RealVector(T), temp: T, time_step: T, dof: usize, ab_initio: bool, allocator: std.mem.Allocator) !T {
    const aligned = if (ab_initio) try alignTrajectory(T, positions, allocator) else positions; defer if (ab_initio) aligned.deinit(allocator);

    var q = try ComplexMatrix(T).initZero(try std.math.powi(usize, 2, std.math.log2_int_ceil(usize, aligned.rows)), aligned.cols, allocator); defer q.deinit(allocator);

    for (0..aligned.rows) |i| for (0..aligned.cols) |j| {q.ptr(i, j).*.re = aligned.at(i, j);};

    for (0..aligned.cols) |j| {

        const mean = aligned.column(j).mean();

        for (0..aligned.rows) |i| q.ptr(i, j).*.re = (aligned.at(i, j) - mean) * std.math.sqrt(masses.at(j));
    }

    var omega = try RealVector(T).initZero(q.rows, allocator); defer omega.deinit(allocator);

    for (0..omega.len) |i| {

        const k = @as(T, @floatFromInt(i)) - if (i < (omega.len + 1) / 2) 0 else @as(T, @floatFromInt(omega.len));

        omega.ptr(i).* = k * (2.0 * std.math.pi) / (@as(T, @floatFromInt(omega.len)) * time_step);
    }

    for (0..q.cols) |j| {
        var column = q.column(j); try cfft1(T, &column, -1);
    }

    var D = try RealVector(T).initZero(q.rows, allocator); defer D.deinit(allocator);

    for (0..q.cols) |j| for (0..q.rows) |i| {D.ptr(i).* += q.at(i, j).squaredMagnitude() * omega.at(i) * omega.at(i);};

    var S: T = 0; var norm: T = 0;

    for (0..D.len) |i| {

        const x = omega.at(i) / (kB * temp / Eh);

        if (x < 1e-6 or x > 1e3) continue;

        S += D.at(i) * kB / Eh * (x / (std.math.exp(x) - 1) - std.math.log(T, std.math.e, 1 - std.math.exp(-x)));

        norm += D.at(i);
    }

    return if (norm > 0) @as(T, @floatFromInt(dof)) * S / norm else 0;
}
