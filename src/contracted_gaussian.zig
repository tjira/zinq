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
        angular: [3]T,
        center: [3]T,
        alpha: []T,
        coef: []T,

        allocator: std.mem.Allocator,

        /// Initialize a contracted Gaussian.
        pub fn init(center: [3]T, angular: [3]T, coef: []const T, alpha: []const T, allocator: std.mem.Allocator) !ContractedGaussian(T) {
            if (coef.len != alpha.len) {
                return throw(ContractedGaussian(T), "THE LENGTH OF THE COEFFICIENTS AND EXPONENTS MUST BE THE SAME", .{});
            }

            const cg = ContractedGaussian(T) {
                .center = center,
                .angular = angular,
                .alpha = try allocator.alloc(T, alpha.len),
                .coef = try allocator.alloc(T, coef.len),
                .allocator = allocator
            };

            for (coef, 0..) |ci, i| cg.coef[i] = ci; for (alpha, 0..) |ai, i| cg.alpha[i] = ai;

            var N: T = 0;

            for (cg.coef, 0..) |ci, i| {

                const pgi = PrimitiveGaussian(T){.center = cg.center, .angular = cg.angular, .alpha = cg.alpha[i]};

                for (cg.coef, 0..) |cj, j| {

                    const pgj = PrimitiveGaussian(T){.center = cg.center, .angular = cg.angular, .alpha = cg.alpha[j]};

                    N += ci * cj * pgi.overlap(pgj) / (pgi.norm() * pgj.norm());
                }
            }

            for(cg.coef) |*ci| ci.* /= std.math.sqrt(N);

            return cg;
        }

        /// Deinitialize a contracted Gaussian.
        pub fn deinit(self: ContractedGaussian(T)) void {
            self.allocator.free(self.coef); self.allocator.free(self.alpha);
        }

        /// Compute the Coulomb integral between four contracted Gaussians.
        pub fn coulomb(self: ContractedGaussian(T), other1: ContractedGaussian(T), other2: ContractedGaussian(T), other3: ContractedGaussian(T)) T {
            var e: T = 0;

            for (self.coef, 0..) |ci, i| {

                const pgi = PrimitiveGaussian(T){.center = self.center, .angular = self.angular, .alpha = self.alpha[i]};

                for (other1.coef, 0..) |cj, j| {

                    const pgj = PrimitiveGaussian(T){.center = other1.center, .angular = other1.angular, .alpha = other1.alpha[j]};

                    for (other2.coef, 0..) |ck, k| {

                        const pgk = PrimitiveGaussian(T){.center = other2.center, .angular = other2.angular, .alpha = other2.alpha[k]};

                        for (other3.coef, 0..) |cl, l| {

                            const pgl = PrimitiveGaussian(T){.center = other3.center, .angular = other3.angular, .alpha = other3.alpha[l]};

                            e += ci * cj * ck * cl * pgi.coulomb(pgj, pgk, pgl) / (pgi.norm() * pgj.norm() * pgk.norm() * pgl.norm());
                        }
                    }
                }
            }

            return e;
        }

        /// Compute the kinetic integral between two contracted Gaussians.
        pub fn kinetic(self: ContractedGaussian(T), other: ContractedGaussian(T)) T {
            var t: T = 0;

            for (self.coef, 0..) |ci, i| {

                const pgi = PrimitiveGaussian(T){.center = self.center, .angular = self.angular, .alpha = self.alpha[i]};

                for (other.coef, 0..) |cj, j| {

                    const pgj = PrimitiveGaussian(T){.center = other.center, .angular = other.angular, .alpha = other.alpha[j]};

                    t += ci * cj * pgi.kinetic(pgj) / (pgi.norm() * pgj.norm());
                }
            }

            return t;
        }

        /// Compute the nuclear integral between two contracted Gaussians.
        pub fn nuclear(self: ContractedGaussian(T), other: ContractedGaussian(T), system: ClassicalParticle(T)) T {
            var v: T = 0;

            for (self.coef, 0..) |ci, i| {

                const pgi = PrimitiveGaussian(T){.center = self.center, .angular = self.angular, .alpha = self.alpha[i]};

                for (other.coef, 0..) |cj, j| {

                    const pgj = PrimitiveGaussian(T){.center = other.center, .angular = other.angular, .alpha = other.alpha[j]};

                    v += ci * cj * pgi.nuclear(pgj, system) / (pgi.norm() * pgj.norm());
                }
            }

            return v;
        }

        /// Compute the overlap integral between two contracted Gaussians.
        pub fn overlap(self: ContractedGaussian(T), other: ContractedGaussian(T)) T {
            var s: T = 0;

            for (self.coef, 0..) |ci, i| {

                const pgi = PrimitiveGaussian(T){.center = self.center, .angular = self.angular, .alpha = self.alpha[i]};

                for (other.coef, 0..) |cj, j| {

                    const pgj = PrimitiveGaussian(T){.center = other.center, .angular = other.angular, .alpha = other.alpha[j]};

                    s += ci * cj * pgi.overlap(pgj) / (pgi.norm() * pgj.norm());
                }
            }

            return s;
        }
    };
}
