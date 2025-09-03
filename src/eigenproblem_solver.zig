//! Solvers for eigenvalue problems.

const std = @import("std");

const real_matrix = @import("real_matrix.zig");

const RealMatrix = real_matrix.RealMatrix;

/// Solve the eigenproblem for a real symmetric system using the Jacobi method.
pub fn diagonalizeJacobi(comptime T: type, A: *RealMatrix(T)) void {
    std.debug.assert(A.rows == A.cols);

    const n = A.rows;

    if (n == 1) return;

    const tol = 16 * std.math.floatEps(T) * A.frobeniusNorm();

    while (A.offDiagonalFrobeniusNorm() > tol) for (0..n - 1) |i| for (i + 1..n) |j| {

        const a = A.at(i, i); const b = A.at(j, j); const g = A.at(i, j);

        if (@abs(g) < std.math.floatEps(T) * std.math.sqrt(a * a + b * b)) continue;

        const tau = 0.5 * (b - a) / g;

        const t = std.math.copysign(@as(T, 1), tau) / (@abs(tau) + std.math.hypot(1, tau));

        const c = 1 / std.math.sqrt(1 + t * t); const s = t * c;

        for (0..n) |k| {

            if (k == i or k == j) continue;

            const u = A.at(k, i);
            const v = A.at(k, j);

            const uk = c * u - s * v;
            const vk = s * u + c * v;

            A.ptr(k, i).* = uk; A.ptr(k, j).* = vk;
            A.ptr(i, k).* = uk; A.ptr(j, k).* = vk;
        }

        A.ptr(i, i).* = a - t * g; A.ptr(i, j).* = 0;
        A.ptr(j, j).* = b + t * g; A.ptr(j, i).* = 0;
    };

    for (0..n - 1) |_| for (0..n - 1) |j| if (A.at(j, j) > A.at(j + 1, j + 1)) {
        std.mem.swap(T, A.ptr(j, j), A.ptr(j + 1, j + 1));
    };
}

/// Diagonalize a real symmetric matrix A. The provided matrix is overwritten by the diagonal form.
pub fn diagonalizeSymmetric(comptime T: type, A: *RealMatrix(T)) !void {
    diagonalizeJacobi(T, A);
}

/// Solve the eigenproblem for a real symmetric system using the Jacobi method.
pub fn eigensystemJacobi(comptime T: type, J: *RealMatrix(T), C: *RealMatrix(T), A: RealMatrix(T)) void {
    std.debug.assert(A.rows == A.cols);
    std.debug.assert(J.rows == A.rows);
    std.debug.assert(J.cols == A.cols);
    std.debug.assert(C.rows == A.rows);
    std.debug.assert(C.cols == A.cols);

    const n = A.rows; @memcpy(J.data, A.data); C.identity();

    if (n == 1) return;

    const tol = 16 * std.math.floatEps(T) * A.frobeniusNorm();

    while (J.offDiagonalFrobeniusNorm() > tol) for (0..n - 1) |i| for (i + 1..n) |j| {

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

        for (0..n) |k| {

            const u = C.at(k, i);
            const v = C.at(k, j);

            C.ptr(k, i).* = c * u - s * v;
            C.ptr(k, j).* = s * u + c * v;
        }
    };

    for (0..n - 1) |_| for (0..n - 1) |j| if (J.at(j, j) > J.at(j + 1, j + 1)) {

        std.mem.swap(T, J.ptr(j, j), J.ptr(j + 1, j + 1));

        for (0..n) |k| std.mem.swap(T, C.ptr(k, j), C.ptr(k, j + 1));
    };
}

/// Solve the eigenproblem for a real symmetric system.
pub fn eigensystemSymmetric(comptime T: type, J: *RealMatrix(T), C: *RealMatrix(T), A: RealMatrix(T)) !void {
    eigensystemJacobi(T, J, C, A);
}

/// Function to fix the gauge of eigenvectors.
pub fn fixGauge(comptime T: type, eigenvectors: *RealMatrix(T), reference: RealMatrix(T)) void {
    std.debug.assert(eigenvectors.rows == reference.rows);
    std.debug.assert(eigenvectors.cols == reference.cols);

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
