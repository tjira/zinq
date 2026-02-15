//! File with determinant calculation functions.

const std = @import("std");

const eigenproblem_solver = @import("eigenproblem_solver.zig");
const real_matrix = @import("real_matrix.zig");

const RealMatrix = real_matrix.RealMatrix;

const eigensystemHermitianAlloc = eigenproblem_solver.eigensystemHermitianAlloc;

/// Calculate the determinant of a symmetric matrix using its eigenvalue decomposition. The matrix must be hermitian (symmetric in the real case).
pub fn determinantHermitian(comptime T: type, A: RealMatrix(T), allocator: std.mem.Allocator) !T {
    const AJC = try eigensystemHermitianAlloc(T, A, allocator); defer AJC.J.deinit(allocator); defer AJC.C.deinit(allocator);

    var det: T = 1;

    for (0..AJC.J.rows) |i| det *= AJC.J.at(i, i);

    return det;
}
