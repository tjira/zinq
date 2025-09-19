//! Integral engine.

const std = @import("std");

const basis_set = @import("basis_set.zig");
const classical_particle = @import("classical_particle.zig");
const contracted_gaussian = @import("contracted_gaussian.zig");
const device_write = @import("device_write.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const matrix_multiplication = @import("matrix_multiplication.zig");
const molecular_integrals = @import("molecular_integrals.zig");
const real_matrix = @import("real_matrix.zig");
const real_tensor_four = @import("real_tensor_four.zig");

const BasisSet = basis_set.BasisSet;
const ClassicalParticle = classical_particle.ClassicalParticle;
const ContractedGaussian = contracted_gaussian.ContractedGaussian;
const RealMatrix = real_matrix.RealMatrix;
const RealTensor4 = real_tensor_four.RealTensor4;

const coulomb = molecular_integrals.coulomb;
const eigensystemSymmetric = eigenproblem_solver.eigensystemSymmetric;
const exportRealMatrix = device_write.exportRealMatrix;
const exportRealTensorFour = device_write.exportRealTensorFour;
const kinetic = molecular_integrals.kinetic;
const mmRealReal = matrix_multiplication.mmRealReal;
const mmRealRealTrans = matrix_multiplication.mmRealRealTrans;
const nuclear = molecular_integrals.nuclear;
const overlap = molecular_integrals.overlap;
const print = device_write.print;
const printJson = device_write.printJson;
const printRealMatrix = device_write.printRealMatrix;

/// Hartree-Fock target options.
pub fn Options(comptime T: type) type {
    return struct {
        const Write = struct {
            coefficient: ?[]const u8 = null,
            density: ?[]const u8 = null,
            fock: ?[]const u8 = null,
        };

        system: []const u8,
        basis: []const u8,

        maxiter: u32 = 100,
        threshold: T = 1e-8,

        write: Write = .{}
    };
}

/// Hartree-Fock target output.
pub fn Output(comptime T: type) type {
    return struct {
        C: RealMatrix(T),
        F: RealMatrix(T),
        K: RealMatrix(T),
        P: RealMatrix(T),
        S: RealMatrix(T),
        V: RealMatrix(T),
        J: RealTensor4(T),

        allocator: std.mem.Allocator,

        /// Deinitialize the output struct.
        pub fn deinit(self: @This()) void {
            self.C.deinit();
            self.F.deinit();
            self.K.deinit();
            self.P.deinit();
            self.S.deinit();
            self.V.deinit();
            self.J.deinit();
        }
    };
}

/// Run the Hartree-Fock target.
pub fn run(comptime T: type, options: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) {try print("\n", .{}); try printJson(options);}

    const system = try classical_particle.read(T, options.system, 0, allocator); defer system.deinit();
    var basis = try BasisSet(T).init(system, options.basis, allocator); defer basis.deinit();

    const S = try overlap(T, basis, allocator);
    const K = try kinetic(T, basis, allocator);
    const V = try nuclear(T, system, basis, allocator);
    const J = try coulomb(T, basis, allocator);

    var C = try RealMatrix(T).init(basis.nbf(), basis.nbf(), allocator);
    var F = try RealMatrix(T).init(basis.nbf(), basis.nbf(), allocator);
    var P = try RealMatrix(T).init(basis.nbf(), basis.nbf(), allocator);

    var FX = try RealMatrix(T).init(basis.nbf(), basis.nbf(), allocator); defer FX.deinit();

    var X = try getXMatrix(T, S, allocator); defer X.deinit();

    const VNN = system.nuclearRepulsionEnergy();

    var energy: T = 0; var energy_prev: T = 1; var iter: usize = 0;

    if (enable_printing) try print("\n{s:4} {s:20} {s:8} {s:4}\n", .{"ITER", "ENERGY", "|DELTA E|", "TIME"});

    while (@abs(energy - energy_prev) > options.threshold) : (iter += 1) {

        if (iter >= options.maxiter) return error.HartreeFockNotConverged;

        var timer = try std.time.Timer.start();

        F.zero();

        for (0..basis.nbf()) |i| for (0..basis.nbf()) |j| {

            F.ptr(i, j).* += V.at(i, j) + K.at(i, j);

            for (0..basis.nbf()) |k| for (0..basis.nbf()) |l| {
                F.ptr(k, l).* += P.at(i, j) * (2.0 * J.at(i, j, k, l) - J.at(i, l, k, j));
            };
        };

        mmRealReal(T, &C, F, X); mmRealReal(T, &FX, X, C);

        try eigensystemSymmetric(T, &P, &C, FX);

        mmRealReal(T, &P, X, C); @memcpy(C.data, P.data);

        for (0..basis.nbf()) |i| for (0..basis.nbf()) |j| {

            P.ptr(i, j).* = 0;

            for (0..system.nocc()) |m| {
                P.ptr(i, j).* += C.at(i, m) * C.at(j, m);
            }
        };

        energy_prev = energy; energy = 0;

        for (0..basis.nbf()) |i| for (0..basis.nbf()) |j| {
            energy += P.at(i, j) * (K.at(i, j) + V.at(i, j) + F.at(i, j));
        };

        if (enable_printing) try print("{d:4} {d:20.14} {e:9.3} {D}\n", .{iter + 1, energy + VNN, @abs(energy - energy_prev), timer.read()});
    }

    if (enable_printing) try print("\nHF ENERGY: {d:.14}\n", .{energy + VNN});

    return .{
        .C = C,
        .F = F,
        .K = K,
        .P = P,
        .S = S,
        .V = V,
        .J = J,
        .allocator = allocator,
    };
}

/// Function to get the X matrix.
pub fn getXMatrix(comptime T: type, S: RealMatrix(T), allocator: std.mem.Allocator) !RealMatrix(T) {
    var X = try RealMatrix(T).init(S.rows, S.cols, allocator);

    var S_C = try RealMatrix(T).init(S.rows, S.cols, allocator); defer S_C.deinit();

    var mm_temp = try RealMatrix(T).init(S.rows, S.cols, allocator); defer mm_temp.deinit();

    try eigensystemSymmetric(T, &X, &S_C, S);

    for (0..S.rows) |i| X.ptr(i, i).* = 1.0 / std.math.sqrt(X.at(i, i));

    mmRealRealTrans(T, &mm_temp, X, S_C);
    mmRealReal(T, &X, S_C, mm_temp);

    return X;
}
