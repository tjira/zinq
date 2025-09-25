//! Solvers for eigenvalue problems.

const std = @import("std");

const error_handling = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");
const real_matrix = @import("real_matrix.zig");

const RealMatrix = real_matrix.RealMatrix;

const throw = error_handling.throw;

const MAX_JACOBI_ITERATIONS = global_variables.MAX_JACOBI_ITERATIONS;

/// Diagonalize a real symmetric matrix A. The provided matrix is overwritten by the diagonal form.
pub fn diagonalizeSymmetric(comptime T: type, A: *RealMatrix(T)) !void {
    try eigensystemJacobi(T, A, null, A.*);
}

/// Solve the eigenproblem for a real symmetric system using the Jacobi method.
pub fn eigensystemJacobi(comptime T: type, J: *RealMatrix(T), C: ?*RealMatrix(T), A: RealMatrix(T)) !void {
    if (!A.isSquare()) return throw(void, "MATRIX MUST BE SQUARE TO SOLVE THE EIGENPROBLEM", .{});
    if (!A.isSymmetric(0)) return throw(void, "THEM MATRIX YOU ARE PASSING TO THE JACOBI EIGENSOLVER IS NOT SYMMETRIC", .{});
    if (A.rows != J.rows or A.cols != J.cols) return throw(void, "EIGENVALUE MATRIX MUST HAVE THE SAME DIMENSIONS AS THE INPUT MATRIX", .{});
    if (C != null and (C.?.rows != A.rows or C.?.cols != A.cols)) return throw(void, "EIGENVECTOR MATRIX MUST HAVE THE SAME DIMENSIONS AS THE INPUT MATRIX", .{});

    const n = A.rows; var iter: usize = 0; if (J.data.ptr != A.data.ptr) try A.copyTo(J);

    if (C != null) C.?.identity(); if (n == 1) return;

    const tol = 16 * std.math.floatEps(T) * A.frobeniusNorm();

    while (J.offDiagonalFrobeniusNorm() > tol) : (iter += 1) for (0..n - 1) |i| for (i + 1..n) |j| {

        if (iter == MAX_JACOBI_ITERATIONS) {
            return throw(void, "JACOBI EIGENSOLVER DID NOT CONVERGE IN {d} ITERATIONS WITH {e:5.3} OFF-DIAGONAL NORM", .{MAX_JACOBI_ITERATIONS, J.offDiagonalFrobeniusNorm()});
        }

        const a = J.at(i, i); const b = J.at(j, j); const g = J.at(i, j);

        if (@abs(g) < std.math.floatEps(T) * std.math.sqrt(a * a + b * b)) continue;

        const tau = 0.5 * (b - a) / g;

        const t = std.math.copysign(@as(T, 1), tau) / (@abs(tau) + std.math.hypot(1, tau));

        const c = 1 / std.math.sqrt(1 + t * t); const s = t * c;

        for (0..n) |k| {

            if (k == i or k == j) continue;

            const u = J.at(k, i);
            const v = J.at(k, j);

            const uk = c * u - s * v;
            const vk = s * u + c * v;

            J.ptr(k, i).* = uk; J.ptr(k, j).* = vk;
            J.ptr(i, k).* = uk; J.ptr(j, k).* = vk;
        }

        J.ptr(i, i).* = a - t * g; J.ptr(i, j).* = 0;
        J.ptr(j, j).* = b + t * g; J.ptr(j, i).* = 0;

        if (C != null) for (0..n) |k| {

            const u = C.?.at(k, i);
            const v = C.?.at(k, j);

            C.?.ptr(k, i).* = c * u - s * v;
            C.?.ptr(k, j).* = s * u + c * v;
        };
    };

    for (0..n) |i| for (0..n) |j| if (i != j) {J.ptr(i, j).* = 0;};

    for (0..n - 1) |_| for (0..n - 1) |j| if (J.at(j, j) > J.at(j + 1, j + 1)) {

        std.mem.swap(T, J.ptr(j, j), J.ptr(j + 1, j + 1));

        if (C != null) for (0..n) |k| {
            std.mem.swap(T, C.?.ptr(k, j), C.?.ptr(k, j + 1));
        };
    };
}

/// Solve the eigenproblem for a real symmetric system.
pub fn eigensystemSymmetric(comptime T: type, J: *RealMatrix(T), C: *RealMatrix(T), A: RealMatrix(T)) !void {
    try eigensystemJacobi(T, J, C, A);
}

/// Solve the eigenproblem for a real symmetric system, returning the eigenvalues and eigenvectors.
pub fn eigensystemSymmetricAlloc(comptime T: type, A: RealMatrix(T), allocator: std.mem.Allocator) !struct {J: RealMatrix(T), C: RealMatrix(T)} {
    var J = try RealMatrix(T).init(A.rows, A.cols, allocator);
    var C = try RealMatrix(T).init(A.rows, A.cols, allocator);

    try eigensystemSymmetric(T, &J, &C, A);

    return .{.J = J, .C = C};
}

/// Function to fix the gauge of eigenvectors.
pub fn fixGauge(comptime T: type, eigenvectors: *RealMatrix(T), reference: RealMatrix(T)) !void {
    if (eigenvectors.rows != reference.rows or eigenvectors.cols != reference.cols) {
        return throw(void, "EIGENVECTOR AND REFERENCE MATRICES MUST HAVE THE SAME DIMENSIONS TO FIX THE GAUGE", .{});
    }

    for (0..eigenvectors.cols) |i| {

        var dot_product: T = 0;

        for (0..eigenvectors.rows) |j| {
            dot_product += eigenvectors.at(j, i) * reference.at(j, i);
        }

        if (dot_product < 0) for (0..eigenvectors.rows) |j| {
            eigenvectors.ptr(j, i).* *= -1;
        };
    }
}
