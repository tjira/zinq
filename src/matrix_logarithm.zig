//! Matrix logarithm evaluations for general and orthogonal matrices using real Schur decomposition.

const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const Vector = @import("tensor.zig").Vector;

const ensureF64 = @import("linear_algebra.zig").ensureF64;
const gees = @import("linear_algebra.zig").gees;
const gemm = @import("linear_algebra.zig").gemm;

/// Computes the skew-symmetric logarithm of an orthogonal matrix using caller-provided workspace memory.
pub fn logmOrthogonalWithWorkspace(comptime T: type, L: *Matrix(T), A: Matrix(T), work: []T) !void {
    std.debug.assert(A.nrow() == A.ncol());

    std.debug.assert(L.nrow() == A.nrow() and L.ncol() == A.ncol());

    std.debug.assert(work.len >= logmWorkspaceSize(A.nrow()));

    return logmSliceWithWorkspace(T, L.data, A.data, work, true);
}

/// Computes the required workspace slice length for real Schur matrix logarithm evaluations.
pub fn logmWorkspaceSize(n: usize) usize {
    return if (n <= 2) 0 else 2 * n * n + 2 * n;
}

/// Computes the principal matrix logarithm of a square matrix using the real Schur decomposition.
fn logm(comptime T: type, L: *Matrix(T), A: Matrix(T), gpa: Allocator) !void {
    std.debug.assert(A.nrow() == A.ncol());

    std.debug.assert(L.nrow() == A.nrow() and L.ncol() == A.ncol());

    return logmSlice(T, L.data, A.data, false, gpa);
}

/// Computes the principal matrix logarithm of a square matrix using caller-provided workspace memory.
fn logmWithWorkspace(comptime T: type, L: *Matrix(T), A: Matrix(T), work: []T) !void {
    std.debug.assert(A.nrow() == A.ncol());

    std.debug.assert(L.nrow() == A.nrow() and L.ncol() == A.ncol());

    std.debug.assert(work.len >= logmWorkspaceSize(A.nrow()));

    return logmSliceWithWorkspace(T, L.data, A.data, work, false);
}

/// Computes the skew-symmetric matrix logarithm of an orthogonal matrix using real Schur decomposition.
fn logmOrthogonal(comptime T: type, L: *Matrix(T), A: Matrix(T), gpa: Allocator) !void {
    std.debug.assert(A.nrow() == A.ncol());

    std.debug.assert(L.nrow() == A.nrow() and L.ncol() == A.ncol());

    return logmSlice(T, L.data, A.data, true, gpa);
}

/// Computes matrix logarithm on contiguous slices by dynamically allocating workspace.
fn logmSlice(comptime T: type, L: []T, A: []const T, comptime orthogonal: bool, gpa: Allocator) !void {
    const n = std.math.sqrt(A.len);

    if (n <= 2) return logmSliceWithWorkspace(T, L, A, &.{}, orthogonal);

    const work = try gpa.alloc(T, logmWorkspaceSize(n));
    defer gpa.free(work);

    return logmSliceWithWorkspace(T, L, A, work, orthogonal);
}

/// Computes matrix logarithm on contiguous slices using real Schur decomposition and workspace memory.
fn logmSliceWithWorkspace(comptime T: type, L: []T, A: []const T, work: []T, comptime orthogonal: bool) !void {
    ensureF64(T, if (orthogonal) "logmOrthogonal" else "logm");

    const n = std.math.sqrt(A.len);

    std.debug.assert(A.len == n * n and L.len == n * n);

    if (n == 0) return;

    if (n == 1) {
        L[0] = try scalarLog(T, A[0], orthogonal);

        return;
    }

    if (n == 2) {
        const ln_r, const theta = try blockLog(T, A[0], A[1], A[2], A[3], orthogonal);

        L[0], L[1], L[2], L[3] = .{ ln_r, -theta, theta, ln_r };

        return;
    }

    std.debug.assert(work.len >= 2 * n * n + 2 * n);

    const nn = n * n;

    const T_slice, const Z_slice = .{ work[0..nn], work[nn .. 2 * nn] };

    const wr = work[2 * nn + 0 .. 2 * nn + 1 * n];
    const wi = work[2 * nn + n .. 2 * nn + 2 * n];

    var T_mat = Matrix(T){ .data = T_slice, .shape = .{ n, n } };
    var Z_mat = Matrix(T){ .data = Z_slice, .shape = .{ n, n } };

    var wr_vec = Vector(T){ .data = wr, .shape = .{n} };
    var wi_vec = Vector(T){ .data = wi, .shape = .{n} };

    const A_mat = Matrix(T){ .data = @constCast(A), .shape = .{ n, n } };

    try gees(T, &T_mat, &Z_mat, &wr_vec, &wi_vec, A_mat);

    var k: usize = 0;

    while (k < n) {
        const wik = wi[k];

        if (wik == 0) {
            const ln_val = try scalarLog(T, T_slice[k * n + k], orthogonal);

            @memset(T_slice[k * n .. (k + 1) * n], 0);

            T_slice[k * n + k] = ln_val;

            k += 1;
        }

        if (wik != 0) {
            const idx0, const idx1 = .{ k * n + k, (k + 1) * n + k };

            const ln_r, const theta = try blockLog(T, T_slice[idx0], T_slice[idx0 + 1], T_slice[idx1], T_slice[idx1 + 1], orthogonal);

            @memset(T_slice[k * n .. (k + 2) * n], 0);

            T_slice[idx0 + 0] = ln_r;

            T_slice[idx0 + 1] = 0 - theta;
            T_slice[idx1 + 0] = 0 + theta;

            T_slice[idx1 + 1] = ln_r;

            k += 2;
        }
    }

    var L_mat = Matrix(T){ .data = L, .shape = .{ n, n } };

    gemm(T, &L_mat, Z_mat, T_mat, 1, 0, false, false);

    gemm(T, &T_mat, L_mat, Z_mat, 1, 0, false, true);

    if (orthogonal) {
        for (0..n) |i| {
            L[i * n + i] = 0;

            for (i + 1..n) |j| {
                const val = 0.5 * (T_slice[i * n + j] - T_slice[j * n + i]);

                L[i * n + j] = 0 + val;
                L[j * n + i] = 0 - val;
            }
        }
    }

    if (!orthogonal) {
        @memcpy(L, T_slice);
    }
}

/// Computes the principal scalar logarithm or returns zero for orthogonal matrices.
fn scalarLog(comptime T: type, val: T, comptime orthogonal: bool) !T {
    if (orthogonal) return 0;

    if (val <= 0) return error.LapackError;

    return std.math.ln(val);
}

/// Computes radial log and rotation angle for 2x2 real quasi-triangular Schur blocks.
fn blockLog(comptime T: type, a: T, b: T, c: T, d: T, comptime orthogonal: bool) !struct { T, T } {
    const theta = std.math.atan2(c - b, a + d);

    if (orthogonal) return .{ 0, theta };

    const det_val = a * d - b * c;

    if (det_val <= 0) return error.LapackError;

    return .{ 0.5 * std.math.ln(det_val), theta };
}
