//! Solvers for eigenvalue problems.

const std = @import("std");

const complex_matrix = @import("complex_matrix.zig");
const device_read = @import("device_read.zig");
const device_write = @import("device_write.zig");
const error_handling = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");
const real_matrix = @import("real_matrix.zig");

const Complex = std.math.complex.Complex;
const ComplexMatrix = complex_matrix.ComplexMatrix;
const RealMatrix = real_matrix.RealMatrix;

const exportComplexMatrix = device_write.exportComplexMatrix;
const exportRealMatrix = device_write.exportRealMatrix;
const printJson = device_write.printJson;
const printComplexMatrix = device_write.printComplexMatrix;
const printRealMatrix = device_write.printRealMatrix;
const readComplexMatrix = device_read.readComplexMatrix;
const readRealMatrix = device_read.readRealMatrix;
const throw = error_handling.throw;

const MAX_JACOBI_ITERATIONS = global_variables.MAX_JACOBI_ITERATIONS;

/// Options for the eigenvalue solver program.
pub fn Options(comptime _: type) type {
    return struct {
        pub const Print = struct {
            input_matrix: bool = true,
            eigenvalues: bool = false,
            eigenvectors: bool = false,
        };
        pub const Write = struct{
            eigenvalues: ?[]const u8 = null,
            eigenvectors: ?[]const u8 = null
        };

        matrix_file: []const u8,
        hermitian: bool = false,
        real: bool = true,

        print: Print = .{},
        write: Write = .{},
    };
}

/// Output of the eigenvalue solver program.
pub fn Output(comptime _: type) type {
    return struct {

        // Deinitialize the output matrices.
        pub fn deinit(_: @This(), _: std.mem.Allocator) void {}
    };
}

/// Run the eigenvalue solver with the provided options.
pub fn run(comptime T: type, opt: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(opt);

    if (!opt.hermitian) return throw(Output(T), "EIGENPROBLEM SOLVER ONLY SUPPORTS HERMITIAN MATRICES CURRENTLY", .{});

    if (opt.real) {

        const A = try readRealMatrix(T, opt.matrix_file, allocator); defer A.deinit(allocator);

        if (enable_printing and opt.print.input_matrix) {
            try device_write.print("\nINPUT MATRIX:\n", .{}); try printRealMatrix(T, A);
        }

        const AJC = try eigensystemHermitianAlloc(T, A, allocator); defer AJC.J.deinit(allocator); defer AJC.C.deinit(allocator);

        if (enable_printing) {
            if (opt.print.eigenvalues) {
                try device_write.print("\nEIGENVALUES:\n", .{}); try printRealMatrix(T, AJC.J);
            }
            if (opt.print.eigenvectors) {
                try device_write.print("\nEIGENVECTORS:\n", .{}); try printRealMatrix(T, AJC.C);
            }
        }

        if (opt.write.eigenvalues) |path| try exportRealMatrix(T, path, AJC.J);
        if (opt.write.eigenvectors) |path| try exportRealMatrix(T, path, AJC.C);
    }

    else {

        const A = try readComplexMatrix(T, opt.matrix_file, allocator); defer A.deinit(allocator);

        if (enable_printing and opt.print.input_matrix) {
            try device_write.print("\nINPUT MATRIX:\n", .{}); try printComplexMatrix(T, A);
        }

        const AJC = try eigensystemHermitianAlloc(T, A, allocator); defer AJC.J.deinit(allocator); defer AJC.C.deinit(allocator);

        if (enable_printing) {
            if (opt.print.eigenvalues) {
                try device_write.print("\nEIGENVALUES:\n", .{}); try printComplexMatrix(T, AJC.J);
            }
            if (opt.print.eigenvectors) {
                try device_write.print("\nEIGENVECTORS:\n", .{}); try printComplexMatrix(T, AJC.C);
            }
        }

        if (opt.write.eigenvalues) |path| try exportComplexMatrix(T, path, AJC.J);
        if (opt.write.eigenvectors) |path| try exportComplexMatrix(T, path, AJC.C);
    }

    return Output(T){};
}

/// Diagonalize a complex hermitian matrix A. The provided matrix is overwritten by the diagonal form.
pub fn diagonalizeHermitian(comptime T: type, A: anytype) !void {
    if (comptime @TypeOf(A.*) == RealMatrix(T)) {try eigensystemJacobi(RealMatrix, T, A, null, A.*);}
    else if (comptime @TypeOf(A.*) == ComplexMatrix(T)) {try eigensystemJacobi(ComplexMatrix, T, A, null, A.*);}
    else @compileError("UNSUPPORTED INPUT MATRIX TYPE IN HERMITIAN DIAGONALIZATION");
}

/// Solve the eigenproblem for a complex hermitian system.
pub fn eigensystemHermitian(comptime T: type, J: anytype, C: anytype, A: anytype) !void {
    if (comptime @TypeOf(A) == RealMatrix(T)) {try eigensystemJacobi(RealMatrix, T, J, C, A);}
    else if (comptime @TypeOf(A) == ComplexMatrix(T)) {try eigensystemJacobi(ComplexMatrix, T, J, C, A);}
    else @compileError("UNSUPPORTED INPUT MATRIX TYPE IN HERMITIAN EIGENSOLVER");
}

/// Solve the eigenproblem for a complex hermitian system, returning the eigenvalues and eigenvectors.
pub fn eigensystemHermitianAlloc(comptime T: type, A: anytype, allocator: std.mem.Allocator) !struct {J: @TypeOf(A), C: @TypeOf(A)} {
    var J = try @TypeOf(A).init(A.rows, A.cols, allocator);
    var C = try @TypeOf(A).init(A.rows, A.cols, allocator);

    try eigensystemHermitian(T, &J, &C, A);

    return .{.J = J, .C = C};
}

/// Solve the eigenproblem for a symmetric or Hermitian system using the Jacobi method.
pub fn eigensystemJacobi(comptime M: fn (comptime type) type, comptime T: type, J: *M(T), C: ?*M(T), A: M(T)) !void {
    if (!A.isSquare()) return throw(void, "MATRIX MUST BE SQUARE TO SOLVE THE EIGENPROBLEM", .{});
    if (A.rows != J.rows or A.cols != J.cols) return throw(void, "EIGENVALUE MATRIX MUST HAVE THE SAME DIMENSIONS AS THE INPUT MATRIX", .{});
    if (C != null and (C.?.rows != A.rows or C.?.cols != A.cols)) return throw(void, "EIGENVECTOR MATRIX MUST HAVE THE SAME DIMENSIONS AS THE INPUT MATRIX", .{});

    if (comptime M == RealMatrix) {
        if (!A.isSymmetric(0)) return throw(void, "THE MATRIX YOU ARE PASSING TO THE JACOBI EIGENSOLVER IS NOT SYMMETRIC", .{});
    } else {
        if (!A.isHermitian(0)) return throw(void, "THE MATRIX YOU ARE PASSING TO THE JACOBI EIGENSOLVER IS NOT HERMITIAN", .{});
    }

    const n = A.rows; var iter: usize = 0; if (J.data.ptr != A.data.ptr) try A.copyTo(J);

    if (C != null) C.?.identity(); if (n == 1) return;

    const tol = 16 * std.math.floatEps(T) * A.frobeniusNorm();

    while (J.offDiagonalFrobeniusNorm() > tol) : (iter += 1) for (0..n - 1) |i| for (i + 1..n) |j| {

        if (iter == MAX_JACOBI_ITERATIONS) {
            return throw(void, "JACOBI EIGENSOLVER DID NOT CONVERGE IN {d} ITERATIONS WITH {e:5.3} OFF-DIAGONAL NORM", .{MAX_JACOBI_ITERATIONS, J.offDiagonalFrobeniusNorm()});
        }

        const a = if (comptime M == RealMatrix) J.at(i, i) else J.at(i, i).re;
        const b = if (comptime M == RealMatrix) J.at(j, j) else J.at(j, j).re;

        const g = if (comptime M == RealMatrix) J.at(i, j) else J.at(i, j).magnitude();

        if (@abs(g) <= std.math.floatEps(T) * std.math.sqrt(a * a + b * b)) continue;

        const tau = 0.5 * (b - a) / g;

        const t = std.math.copysign(@as(T, 1), tau) / (@abs(tau) + std.math.hypot(1, tau));

        const c = if (comptime M == RealMatrix) 1 / std.math.sqrt(1 + t * t) else Complex(T).init(1 / std.math.sqrt(1 + t * t), 0);

        const s = if (comptime M == RealMatrix) t * c else std.math.complex.exp(Complex(T).init(0, std.math.complex.arg(J.at(i, j)))).mul(c).mul(Complex(T).init(t, 0));

        for (0..n) |k| {

            if (k == i or k == j) continue;

            const u = J.at(k, i);
            const v = J.at(k, j);

            const uk = if (comptime M == RealMatrix) c * u - s * v else c.mul(u).sub(s.conjugate().mul(v));
            const vk = if (comptime M == RealMatrix) s * u + c * v else s.mul(u).add(c            .mul(v));

            J.ptr(k, i).* = uk; J.ptr(i, k).* = if (comptime M == RealMatrix) uk else uk.conjugate();
            J.ptr(k, j).* = vk; J.ptr(j, k).* = if (comptime M == RealMatrix) vk else vk.conjugate();
        }

        J.ptr(i, i).* = if (comptime M == RealMatrix) a - t * g else Complex(T).init(a - t * g, 0); J.ptr(i, j).* = if (comptime M == RealMatrix) 0 else Complex(T).init(0, 0);
        J.ptr(j, j).* = if (comptime M == RealMatrix) b + t * g else Complex(T).init(b + t * g, 0); J.ptr(j, i).* = if (comptime M == RealMatrix) 0 else Complex(T).init(0, 0);

        if (C != null) for (0..n) |k| {

            const u = C.?.at(k, i);
            const v = C.?.at(k, j);

            C.?.ptr(k, i).* = if (comptime M == RealMatrix) c * u - s * v else c.mul(u).sub(s.conjugate().mul(v));
            C.?.ptr(k, j).* = if (comptime M == RealMatrix) s * u + c * v else s.mul(u).add(c            .mul(v));
        };
    };

    for (0..n) |i| for (0..n) |j| if (i != j) {J.ptr(i, j).* = if (comptime M == RealMatrix) 0 else Complex(T).init(0, 0);};

    for (0..n - 1) |_| for (0..n - 1) |j| {

        const Jii = if (comptime M == RealMatrix) J.at(j,     j    ) else J.at(j,     j    ).re;
        const Jjj = if (comptime M == RealMatrix) J.at(j + 1, j + 1) else J.at(j + 1, j + 1).re;

        if (Jii > Jjj) {

            if (comptime M == RealMatrix) {
                std.mem.swap(T, J.ptr(j, j), J.ptr(j + 1, j + 1));
            } else {
                std.mem.swap(T, &J.ptr(j, j).re, &J.ptr(j + 1, j + 1).re);
                std.mem.swap(T, &J.ptr(j, j).im, &J.ptr(j + 1, j + 1).im);
            }

            if (C != null) for (0..n) |k| {

                if (comptime M == RealMatrix) {
                    std.mem.swap(T, C.?.ptr(k, j), C.?.ptr(k, j + 1));
                } else {
                    std.mem.swap(T, &C.?.ptr(k, j).re, &C.?.ptr(k, j + 1).re);
                    std.mem.swap(T, &C.?.ptr(k, j).im, &C.?.ptr(k, j + 1).im);
                }
            };
        }
    };
}

/// Solve the eigenproblem for a complex hermitian system.
pub fn eigensystemJacobiAlloc(comptime M: fn (comptime type) type, comptime T: type, A: M(T), allocator: std.mem.Allocator) !struct {J: M(T), C: M(T)} {
    var J = try M(T).init(A.rows, A.cols, allocator);
    var C = try M(T).init(A.rows, A.cols, allocator);

    try eigensystemJacobi(M, T, &J, &C, A);

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
