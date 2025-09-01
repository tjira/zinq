//! Solvers for eigenvalue problems.

const std = @import("std");

const global_variables = @import("global_variables.zig");
const real_matrix = @import("real_matrix.zig");
const math_functions = @import("math_functions.zig");

const RealMatrix = real_matrix.RealMatrix;

const sgn = math_functions.sgn;

const JACOBI_EIGENPROBLEM_TOLERANCE = global_variables.JACOBI_EIGENPROBLEM_TOLERANCE;

/// Diagonalize a real symmetric matrix A. The provided matrix is overwritten by the diagonal form.
pub fn diagonalizeSymmetric(comptime T: type, A: *RealMatrix(T)) !void {
    if (A.rows == 1) return;

    if (A.rows == 2) {
        diagonalizeSymmetric2x2(T, A);
    }

    else return error.NotImplemented;
}

/// Diagonalize a real symmetric 2x2 matrix A using an analytical formula. The provided matrix is overwritten by the diagonal form.
pub fn diagonalizeSymmetric2x2(comptime T: type, A: *RealMatrix(T)) void {
    const a00 = A.at(0, 0); const a01 = A.at(0, 1); const a11 = A.at(1, 1);

    const tol = 16 * @max(@max(@abs(a00), @abs(a11)), @abs(a01)) * std.math.floatEps(T);

    if (@abs(a01) < tol) {

        A.ptr(0, 0).* = if (a00 < a11) a00 else a11;
        A.ptr(1, 1).* = if (a00 < a11) a11 else a00;

        A.ptr(0, 1).* = 0;
        A.ptr(1, 0).* = 0;

        return;
    }

    const tau = 0.5 * (a11 - a00) / a01;

    const t = std.math.copysign(@as(T, 1), tau) / (@abs(tau) + std.math.hypot(1, tau));

    A.ptr(0, 0).* = a00 - t * a01;
    A.ptr(1, 1).* = a11 + t * a01;

    A.ptr(0, 1).* = 0;
    A.ptr(1, 0).* = 0;
}

/// Solve the eigenproblem for a real symmetric system.
pub fn eigensystemSymmetric(comptime T: type, A_eigenvalues: *RealMatrix(T), A_eigenvectors: *RealMatrix(T), A: RealMatrix(T)) !void {
    if (A.rows == 1) {
        @memcpy(A_eigenvalues.data, A.data); A_eigenvectors.fill(1);
    }

    else if (A.rows == 2) {
        eigensystemSymmetric2x2(T, A_eigenvalues, A_eigenvectors, A);
    }

    else return error.NotImplemented;
}

/// Eigenproblem for a real symmetric 2x2 system using an analytical formula.
pub fn eigensystemSymmetric2x2(comptime T: type, A_eigenvalues: *RealMatrix(T), A_eigenvectors: *RealMatrix(T), A: RealMatrix(T)) void {
    const a00 = A.at(0, 0); const a01 = A.at(0, 1); const a11 = A.at(1, 1);

    const tol = 16 * @max(@max(@abs(a00), @abs(a11)), @abs(a01)) * std.math.floatEps(T);

    if (@abs(a01) < tol) {

        A_eigenvalues.ptr(0, 0).* = if (a00 < a11) a00 else a11;
        A_eigenvalues.ptr(1, 1).* = if (a00 < a11) a11 else a00;

        A_eigenvectors.ptr(0, 0).* = if (a00 < a11) 1 else 0;
        A_eigenvectors.ptr(0, 1).* = if (a00 < a11) 0 else 1;
        A_eigenvectors.ptr(1, 0).* = if (a00 < a11) 0 else 1;
        A_eigenvectors.ptr(1, 1).* = if (a00 < a11) 1 else 0;

        return;
    }

    const tau = 0.5 * (a11 - a00) / a01;

    const t = std.math.copysign(@as(T, 1), tau) / (@abs(tau) + std.math.hypot(1, tau));

    const c = 1 / std.math.sqrt(1 + t * t); const s = t * c;

    A_eigenvalues.ptr(0, 0).* = a00 - t * a01;
    A_eigenvalues.ptr(1, 1).* = a11 + t * a01;

    A_eigenvalues.ptr(0, 1).* = 0;
    A_eigenvalues.ptr(1, 0).* = 0;

    A_eigenvectors.ptr(0, 0).* =  c;
    A_eigenvectors.ptr(0, 1).* =  s;
    A_eigenvectors.ptr(1, 0).* = -s;
    A_eigenvectors.ptr(1, 1).* =  c;

    if (A_eigenvalues.at(0, 0) > A_eigenvalues.at(1, 1)) {

        std.mem.swap(T, A_eigenvalues.ptr(0, 0), A_eigenvalues.ptr(1, 1));

        std.mem.swap(T, A_eigenvectors.ptr(0, 0), A_eigenvectors.ptr(0, 1));
        std.mem.swap(T, A_eigenvectors.ptr(1, 0), A_eigenvectors.ptr(1, 1));
    }
}
