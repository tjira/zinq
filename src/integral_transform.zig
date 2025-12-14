//! Molecular integral transform file.

const std = @import("std");

const real_matrix = @import("real_matrix.zig");
const real_tensor_four = @import("real_tensor_four.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealTensor4 = real_tensor_four.RealTensor4;

/// Extend one-electron integral matrices to spinorbital basis.
pub fn oneAO2AS(comptime T: type, A: *RealMatrix(T), allocator: std.mem.Allocator) !void {
    try A.expand(2 * A.rows, 2 * A.cols, allocator);

    for (0..A.rows / 2) |i| for (0..A.cols / 2) |j| {
        A.ptr(A.rows / 2 + i, A.cols / 2 + j).* = A.at(i, j);
    };
}

pub fn oneAO2MO(comptime T: type, A_MO: *RealMatrix(T), A_AO: RealMatrix(T), C_MO: RealMatrix(T)) void {
    for (0..C_MO.cols) |p| for (0..p + 1) |q| {

        var sum: T = 0;

        for (0..C_MO.rows) |mu| {

            sum += C_MO.at(mu, p) * A_AO.at(mu, mu) * C_MO.at(mu, q);

            for (0..mu) |nu| {
                sum += A_AO.at(mu, nu) * (C_MO.at(mu, p) * C_MO.at(nu, q) + C_MO.at(nu, p) * C_MO.at(mu, q));
            }
        }

        A_MO.ptr(q, p).* = sum;

        if (p != q) A_MO.ptr(p, q).* = sum;
    };
}

pub fn oneAO2MS(comptime T: type, A_MS: *RealMatrix(T), A_AO: RealMatrix(T), C_MO: RealMatrix(T)) void {
    A_MS.zero();

    for (0..C_MO.cols) |p| for (0..p + 1) |q| {
            
        var sum: T = 0;

        for (0..C_MO.rows) |mu| {

            sum += C_MO.at(mu, p) * A_AO.at(mu, mu) * C_MO.at(mu, q);

            for (0..mu) |nu| {
                sum += A_AO.at(mu, nu) * (C_MO.at(mu, p) * C_MO.at(nu, q) + C_MO.at(nu, p) * C_MO.at(mu, q));

            }
        }

        A_MS.ptr(2 * p, 2 * q).* = sum; A_MS.ptr(2 * p + 1, 2 * q + 1).* = sum;

        if (p != q) {
            A_MS.ptr(2 * q, 2 * p).* = sum; A_MS.ptr(2 * q + 1, 2 * p + 1).* = sum;
        }
    };
}

/// Extend two-electron integral tensor to spinorbital basis.
pub fn twoAO2AS(comptime T: type, A: *RealTensor4(T), allocator: std.mem.Allocator) !void {
    try A.expand([4]usize{2 * A.shape[0], 2 * A.shape[1], 2 * A.shape[2], 2 * A.shape[3]}, allocator);

    for (0..A.shape[0] / 2) |i| for (0..A.shape[1] / 2) |j| for (0..A.shape[2] / 2) |k| for (0..A.shape[3] / 2) |l| {
        A.ptr(i,                  j,                  k + A.shape[2] / 2, l + A.shape[3] / 2).* = A.at(i, j, k, l);
        A.ptr(i + A.shape[0] / 2, j + A.shape[1] / 2, k,                  l                 ).* = A.at(i, j, k, l);
        A.ptr(i + A.shape[0] / 2, j + A.shape[1] / 2, k + A.shape[2] / 2, l + A.shape[3] / 2).* = A.at(i, j, k, l);
    };
}

/// Transforms the two-electron integrals from the AO basis to the MO basis.
pub fn twoAO2MO(comptime T: type, A_MO: *RealTensor4(T), A_AO: RealTensor4(T), C_MO: RealMatrix(T)) void {
    A_MO.zero();

    for (0..C_MO.rows) |mu| for (0..C_MO.rows) |nu| for (0..C_MO.cols) |r| for (0..C_MO.cols) |s| {

        var value: T = 0;

        for (0..C_MO.rows) |rho| for (0..C_MO.rows) |sig| {
            value += C_MO.at(rho, r) * C_MO.at(sig, s) * A_AO.at(mu, nu, rho, sig);
        };

        if (value == 0) continue;

        for (0..C_MO.cols) |p| for (0..C_MO.cols) |q| {
             A_MO.ptr(p, q, r, s).* += C_MO.at(mu, p) * C_MO.at(nu, q) * value;
        };
    };
}

/// Transforms the two-electron integrals from the AO basis to the MS basis.
pub fn twoAO2MS(comptime T: type, A_MS: *RealTensor4(T), A_AO: RealTensor4(T), C_MO: RealMatrix(T)) void {
    A_MS.zero();

    for (0..C_MO.rows) |mu| for (0..C_MO.rows) |nu| for (0..C_MO.cols) |r| for (0..C_MO.cols) |s| {

        var value: T = 0;

        for (0..C_MO.rows) |rho| for (0..C_MO.rows) |sig| {
            value += C_MO.at(rho, r) * C_MO.at(sig, s) * A_AO.at(mu, nu, rho, sig);
        };

        if (value == 0) continue;

        for (0..C_MO.cols) |p| for (0..C_MO.cols) |q| {

            const contribution = C_MO.at(mu, p) * C_MO.at(nu, q) * value;
            
            A_MS.ptr(2 * p + 0, 2 * q + 0, 2 * r + 0, 2 * s + 0).* += contribution;
            A_MS.ptr(2 * p + 0, 2 * q + 0, 2 * r + 1, 2 * s + 1).* += contribution;
            A_MS.ptr(2 * p + 1, 2 * q + 1, 2 * r + 0, 2 * s + 0).* += contribution;
            A_MS.ptr(2 * p + 1, 2 * q + 1, 2 * r + 1, 2 * s + 1).* += contribution;
        };
    };
}
