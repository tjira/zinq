//! Contracted Gaussian type definitions.

const std = @import("std");

const array_functions = @import("array_functions.zig");
const classical_particle = @import("classical_particle.zig");
const error_handling = @import("error_handling.zig");
const primitive_gaussian = @import("primitive_gaussian.zig");

const ClassicalParticle = classical_particle.ClassicalParticle;
const PrimitiveGaussian = primitive_gaussian.PrimitiveGaussian;

const sum = array_functions.sum;
const throw = error_handling.throw;

/// Contracted Gaussian type.
pub fn ContractedGaussian(comptime T: type) type {
    return struct {
        A: [3]T,
        a: [3]T,
        alpha: []T,
        c: []T,
        l: i8,

        allocator: std.mem.Allocator,

        /// Initialize a contracted Gaussian.
        pub fn init(A: [3]T, a: [3]T, c: []const T, alpha: []const T, allocator: std.mem.Allocator) !ContractedGaussian(T) {
            if (c.len != alpha.len) {
                return throw(ContractedGaussian(T), "THE LENGTH OF THE COEFFICIENTS AND EXPONENTS MUST BE THE SAME", .{});
            }

            const cg = ContractedGaussian(T) {
                .A = A,
                .a = a,
                .alpha = try allocator.alloc(T, alpha.len),
                .c = try allocator.alloc(T, c.len),
                .l = @as(i8, @intFromFloat(sum(T, &a))),
                .allocator = allocator
            };

            var N: T = 0;

            for (c, 0..) |ci, i| cg.c[i] = ci; for (alpha, 0..) |ai, i| cg.alpha[i] = ai;

            for (cg.c, 0..) |ci, i| for (cg.c, 0..) |cj, j| {

                const pgi = PrimitiveGaussian(T){.A = cg.A, .a = cg.a, .alpha = cg.alpha[i], .l = cg.l};
                const pgj = PrimitiveGaussian(T){.A = cg.A, .a = cg.a, .alpha = cg.alpha[j], .l = cg.l};

                N += ci * cj * pgi.overlap(pgj) / (pgi.norm() * pgj.norm());
            };

            for(cg.c) |*ci| ci.* /= std.math.sqrt(N);

            return cg;
        }

        /// Deinitialize a contracted Gaussian.
        pub fn deinit(self: ContractedGaussian(T)) void {
            self.allocator.free(self.c); self.allocator.free(self.alpha);
        }

        /// Compute the Coulomb integral between four contracted Gaussians.
        pub fn coulomb(self: ContractedGaussian(T), other1: ContractedGaussian(T), other2: ContractedGaussian(T), other3: ContractedGaussian(T)) T {
            var e: T = 0;

            for (self.c, 0..) |ci, i| for (other1.c, 0..) |cj, j| for (other2.c, 0..) |ck, k| for (other3.c, 0..) |cl, l| {

                const pgi = PrimitiveGaussian(T){.A =   self.A, .a =   self.a, .alpha =   self.alpha[i], .l =   self.l};
                const pgj = PrimitiveGaussian(T){.A = other1.A, .a = other1.a, .alpha = other1.alpha[j], .l = other1.l};
                const pgk = PrimitiveGaussian(T){.A = other2.A, .a = other2.a, .alpha = other2.alpha[k], .l = other2.l};
                const pgl = PrimitiveGaussian(T){.A = other3.A, .a = other3.a, .alpha = other3.alpha[l], .l = other3.l};

                e += ci * cj * ck * cl * pgi.coulomb(pgj, pgk, pgl) / (pgi.norm() * pgj.norm() * pgk.norm() * pgl.norm());
            };

            return e;
        }

        /// Compute the kinetic integral between two contracted Gaussians.
        pub fn kinetic(self: ContractedGaussian(T), other: ContractedGaussian(T)) T {
            var t: T = 0;

            for (self.c, 0..) |ci, i| for (other.c, 0..) |cj, j| {

                const pgi = PrimitiveGaussian(T){.A =  self.A, .a =  self.a, .alpha =  self.alpha[i], .l =  self.l};
                const pgj = PrimitiveGaussian(T){.A = other.A, .a = other.a, .alpha = other.alpha[j], .l = other.l};

                t += ci * cj * pgi.kinetic(pgj) / (pgi.norm() * pgj.norm());
            };

            return t;
        }

        /// Compute the nuclear integral between two contracted Gaussians.
        pub fn nuclear(self: ContractedGaussian(T), other: ContractedGaussian(T), system: ClassicalParticle(T)) T {
            var v: T = 0;

            for (self.c, 0..) |ci, i| for (other.c, 0..) |cj, j| {

                const pgi = PrimitiveGaussian(T){.A =  self.A, .a =  self.a, .alpha =  self.alpha[i], .l =  self.l};
                const pgj = PrimitiveGaussian(T){.A = other.A, .a = other.a, .alpha = other.alpha[j], .l = other.l};

                v += ci * cj * pgi.nuclear(pgj, system) / (pgi.norm() * pgj.norm());
            };

            return v;
        }

        /// Compute the overlap integral between two contracted Gaussians.
        pub fn overlap(self: ContractedGaussian(T), other: ContractedGaussian(T)) T {
            var s: T = 0;

            for (self.c, 0..) |ci, i| for (other.c, 0..) |cj, j| {

                const pgi = PrimitiveGaussian(T){.A =  self.A, .a =  self.a, .alpha =  self.alpha[i], .l =  self.l};
                const pgj = PrimitiveGaussian(T){.A = other.A, .a = other.a, .alpha = other.alpha[j], .l = other.l};

                s += ci * cj * pgi.overlap(pgj) / (pgi.norm() * pgj.norm());
            };

            return s;
        }
    };
}
