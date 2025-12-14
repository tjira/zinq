//! Molecular integral transform file.

const std = @import("std");

const matrix_multiplication = @import("matrix_multiplication.zig");
const real_matrix = @import("real_matrix.zig");
const real_tensor_four = @import("real_tensor_four.zig");

const RealMatrix = real_matrix.RealMatrix;
const RealTensor4 = real_tensor_four.RealTensor4;

const mm = matrix_multiplication.mm;

/// Transforms coefficients obtained from the HF calculation from the MO basis to the MS basis.
pub fn cfsAO2AS(comptime T: type, C_MS: *RealMatrix(T), C_MO: RealMatrix(T)) void {
    C_MS.zero();

    for (0..C_MO.rows) |i| for (0..C_MO.cols) |j| {
        C_MS.ptr(i            , 2 * j    ).* = C_MO.at(i, j);
        C_MS.ptr(i + C_MO.rows, 2 * j + 1).* = C_MO.at(i, j);
    };
}

/// Extend one-electron integral matrices to spinorbital basis.
pub fn oneAO2AS(comptime T: type, A_AS: *RealMatrix(T), A_AO: RealMatrix(T)) void {
    A_AS.zero();

    for (0..A_AO.rows) |i| for (0..A_AO.cols) |j| {
        A_AS.ptr(i            , j            ).* = A_AO.at(i, j);
        A_AS.ptr(i + A_AO.rows, j + A_AO.cols).* = A_AO.at(i, j);
    };
}

/// Transforms the one-electron integrals from the AO basis to the MO basis or from the AS basis to the MS basis.
pub fn oneAO2MO(comptime T: type, A_MO: *RealMatrix(T), A_AO: RealMatrix(T), C_MO: RealMatrix(T), allocator: std.mem.Allocator) !void {
    var AC = try RealMatrix(T).init(A_AO.rows, A_AO.cols, allocator); defer AC.deinit(allocator);

    try mm(T, &AC,  A_AO, false, C_MO, false);
    try mm(T, A_MO, C_MO, true,  AC,   false);
}

/// Transforms the one-electron integrals from the AO basis to the MS basis.
pub fn oneAO2MS(comptime T: type, A_MS: *RealMatrix(T), A_AO: RealMatrix(T), C_MO: RealMatrix(T), allocator: std.mem.Allocator) !void {
    var C_MS = try RealMatrix(T).init(2 * C_MO.rows, 2 * C_MO.cols, allocator); defer C_MS.deinit(allocator);
    var A_AS = try RealMatrix(T).init(    A_MS.rows,     A_MS.cols, allocator); defer A_AS.deinit(allocator);

    cfsAO2AS(T, &C_MS, C_MO); oneAO2AS(T, A_MS, A_AO);

    try mm(T, &A_AS, A_MS.*, false, C_MS, false);
    try mm(T, A_MS,  C_MS,   true,  A_AS, false);
}

/// Extend two-electron integral tensor to spinorbital basis.
pub fn twoAO2AS(comptime T: type, A_AS: *RealTensor4(T), A_AO: RealTensor4(T)) void {
    A_AS.zero();

    for (0..A_AO.shape[0]) |i| for (0..A_AO.shape[1]) |j| for (0..A_AO.shape[2]) |k| for (0..A_AO.shape[3]) |l| {
        A_AS.ptr(i,                 j,                 k                , l                ).* = A_AO.at(i, j, k, l);
        A_AS.ptr(i,                 j                , k + A_AO.shape[2], l + A_AO.shape[3]).* = A_AO.at(i, j, k, l);
        A_AS.ptr(i + A_AO.shape[0], j + A_AO.shape[1], k,                 l                ).* = A_AO.at(i, j, k, l);
        A_AS.ptr(i + A_AO.shape[0], j + A_AO.shape[1], k + A_AO.shape[2], l + A_AO.shape[3]).* = A_AO.at(i, j, k, l);
    };
}

/// Transforms the two-electron integrals from the AO basis to the MO basis or from the AS basis to the MS basis.
pub fn twoAO2MO(comptime T: type, A_MO: *RealTensor4(T), A_AO: RealTensor4(T), C_MO: RealMatrix(T)) void {
    for (0..A_MO.shape[0]) |i| for (0..A_MO.shape[1]) |j| for (0..A_MO.shape[2]) |k| for (0..A_MO.shape[3]) |l| {

        var sum: T = 0;

        for (0..C_MO.rows) |p| for (0..C_MO.rows) |q| for (0..C_MO.rows) |r| for (0..C_MO.rows) |s| {
            sum += C_MO.at(p, i) * C_MO.at(q, j) * C_MO.at(r, k) * C_MO.at(s, l) * A_AO.at(p, q, r, s);
        };

        A_MO.ptr(i, j, k, l).* = sum;
    };
}

/// Transforms the two-electron integrals from the AO basis to the MS basis.
pub fn twoAO2MS(comptime T: type, A_MS: *RealTensor4(T), A_AO: RealTensor4(T), C_MO: RealMatrix(T), allocator: std.mem.Allocator) !void {
    var C_MS = try RealMatrix (T).init(2 * C_MO.rows, 2 * C_MO.cols, allocator); defer C_MS.deinit(allocator);
    var A_AS = try RealTensor4(T).init(A_MS.shape,                   allocator); defer A_AS.deinit(allocator);

    cfsAO2AS(T, &C_MS, C_MO); twoAO2AS(T, &A_AS, A_AO);

    for (0..A_MS.shape[0]) |i| for (0..A_MS.shape[1]) |j| for (0..A_MS.shape[2]) |k| for (0..A_MS.shape[3]) |l| {

        var sum: T = 0;

        for (0..C_MS.rows) |p| for (0..C_MS.rows) |q| for (0..C_MS.rows) |r| for (0..C_MS.rows) |s| {
            sum += C_MS.at(p, i) * C_MS.at(q, j) * C_MS.at(r, k) * C_MS.at(s, l) * A_AS.at(p, q, r, s);
        };

        A_MS.ptr(i, j, k, l).* = sum;
    };
}
