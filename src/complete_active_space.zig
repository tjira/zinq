//! This module contains functions to generate determinants for complete active space (CAS) calculations in quantum chemistry. It includes functions to generate all possible determinants, as well as determinants with only singlet excitations.

const std = @import("std");

const math_functions = @import("math_functions.zig");
const real_matrix = @import("real_matrix.zig");

const comb = math_functions.comb;
const combinations = math_functions.combinations;

const RealMatrix = real_matrix.RealMatrix;

/// Generates the determinants for the CI calculations including singlet and triplet excitations.
pub fn generateAllDeterminants(nbf: usize, nel: usize, allocator: std.mem.Allocator) !RealMatrix(usize) {
    var D = try RealMatrix(usize).init(comb(nbf, nel), nel, allocator);

    for (0..D.cols) |i| D.ptr(0, i).* = i;

    for (1..D.rows) |i| {

        var row_i = D.row(i); try D.row(i - 1).copyTo(&row_i); var index: usize = 0;

        for (0..D.cols) |j| if (D.at(i, D.cols - j - 1) != nbf - j - 1) {
            D.ptr(i, D.cols - j - 1).* += 1; index = D.cols - j - 1; break;
        };

        for (index + 1..D.cols) |j| D.ptr(i, j).* = D.at(i, j - 1) + 1;
    }

    return D;
}

/// Generates CAS determinants for the CI calculations.
pub fn generateCasDeterminants(nel: usize, nbf: usize, nocc: usize, spin: bool, allocator: std.mem.Allocator) !RealMatrix(usize) {
    const CASD = if (spin) try generateAllDeterminants(nbf, nel, allocator) else try generateAllSingletDeterminants(nbf, nel, allocator); defer CASD.deinit(allocator);

    var D = try RealMatrix(usize).init(CASD.rows, nocc, allocator);

    for (0..CASD.rows) |i| {
        for (0..nocc - nel) |j| {D.ptr(i, j).* = j;} for (0..CASD.cols) |j| D.ptr(i, j + nocc - nel).* = CASD.at(i, j) + nocc - nel;
    }

    return D;
}

/// Generates the determinants for the CI calculations including only singlet excitations.
pub fn generateAllSingletDeterminants(nbf: usize, nel: usize, allocator: std.mem.Allocator) !RealMatrix(usize) {
    const data_alpha = try allocator.alloc(usize, nbf / 2); defer allocator.free(data_alpha);
    const data_beta  = try allocator.alloc(usize, nbf / 2); defer allocator.free(data_beta );
    
    for (0..data_alpha.len) |i| data_alpha[i] = 2 * i + 0;
    for (0..data_beta.len ) |i| data_beta[i]  = 2 * i + 1;

    var dets_alpha = try combinations(usize, data_alpha, nel / 2, allocator); defer dets_alpha.deinit(allocator);
    var dets_beta  = try combinations(usize, data_beta,  nel / 2, allocator); defer  dets_beta.deinit(allocator);

    var D = try RealMatrix(usize).init(dets_alpha.items.len * dets_beta.items.len, nel, allocator);

    for (0..dets_alpha.items.len) |i| for (0..dets_beta.items.len) |j| for (0..nel / 2) |k| {
        D.ptr(i * dets_beta.items.len + j, k          ).* = dets_alpha.items[i][k];
        D.ptr(i * dets_beta.items.len + j, k + nel / 2).* =  dets_beta.items[j][k];
    };

    for (dets_alpha.items) |*e| allocator.free(e.*);
    for (dets_beta.items ) |*e| allocator.free(e.*);

    return D;
}
