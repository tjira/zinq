//! File with functions related to the frequency analysis of a system.

const std = @import("std");

const classical_particle = @import("classical_particle.zig");
const eigenproblem_solver = @import("eigenproblem_solver.zig");
const error_handling = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");
const real_matrix = @import("real_matrix.zig");
const real_vector = @import("real_vector.zig");

const ClassicalParticle = classical_particle.ClassicalParticle;
const RealMatrix = real_matrix.RealMatrix;
const RealVector = real_vector.RealVector;

const throw = error_handling.throw;
const eigensystemSymmetricAlloc = eigenproblem_solver.eigensystemSymmetricAlloc;

const AN2M = global_variables.AN2M;
const Eh = global_variables.Eh;
const a0 = global_variables.a0;
const amu = global_variables.amu;
const c = global_variables.c;

/// Returns the harmonic frequencies of the system.
pub fn particleHarmonicFrequencies(comptime T: type, system: ClassicalParticle(T), hessian: RealMatrix(T), allocator: std.mem.Allocator) !RealVector(T) {
    var freqs = try RealVector(T).init(hessian.rows, allocator); var HM = try hessian.clone(); defer HM.deinit();

    for (system.atoms.?) |atom| {if (atom - 1 > AN2M.len) return throw(RealVector(T), "ATOMIC NUMBER OUT OF RANGE IN FREQUENCY ANALYSIS", .{});}

    for (0..hessian.rows) |i| for (0..hessian.rows) |j| {
        HM.ptr(i, j).* /= std.math.sqrt(AN2M[system.atoms.?[i / 3]] * AN2M[system.atoms.?[j / 3]]);
    };

    const HMJC = try eigensystemSymmetricAlloc(T, HM, allocator); defer HMJC.J.deinit(); defer HMJC.C.deinit();

    const factor = 5e-3 / std.math.pi / c * std.math.sqrt(Eh / amu / a0 / a0);

    for (0..freqs.len) |i| freqs.ptr(i).* = factor * std.math.sign(HMJC.J.at(i, i)) * std.math.sqrt(@abs(HMJC.J.at(i, i)));

    return freqs;
}
