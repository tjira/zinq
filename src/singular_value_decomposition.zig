//! Finds singular value decomposition of a matrix using the eigenvalue decomposition of the covariance matrix.

const std = @import("std");

const eigenproblem_solver = @import("eigenproblem_solver.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const real_matrix = @import("real_matrix.zig");

const eigensystemHermitianAlloc = eigenproblem_solver.eigensystemHermitianAlloc;
const mmAlloc = matrix_multiplication.mmAlloc;

const RealMatrix = real_matrix.RealMatrix;

pub fn svdAlloc(comptime T: type, A: RealMatrix(T), allocator: std.mem.Allocator) !struct {U: RealMatrix(T), S: RealMatrix(T), VT: RealMatrix(T)} {
    if (!A.isSquare()) {

        std.log.err("SVD IS NOT IMPLEMENTED FOR NON-SQUARE MATRICES", .{});

        return error.NotImplemented;
    }

    const AAT = try mmAlloc(T, A, false, A, true, allocator);
    const ATA = try mmAlloc(T, A, true, A, false, allocator);

    const AATJC = try eigensystemHermitianAlloc(T, AAT, allocator);
    const ATAJC = try eigensystemHermitianAlloc(T, ATA, allocator);

    defer {
        AAT.deinit(allocator); ATA.deinit(allocator); ATAJC.J.deinit(allocator);
    }

    var U = AATJC.C; var S = AATJC.J; var VT = ATAJC.C; try VT.transpose();

    for (0..S.rows) |i| {
        S.ptr(i, i).* = if (S.at(i, i) > 0) std.math.sqrt(S.at(i, i)) else 0;
    }

    for (0..A.cols) |k| {

        var dot: T = 0;

        for (0..A.rows) |i| {

            var sum: T = 0;

            for (0..A.cols) |j| {
                sum += A.at(i, j) * VT.at(k, j); 
            }

            dot += sum * U.at(i, k); 
        }

        if (dot < 0) for (0..U.rows) |i| {
            U.ptr(i, k).* = -U.at(i, k);
        };
    }

    return .{.U = U, .S = S, .VT = VT};
}
