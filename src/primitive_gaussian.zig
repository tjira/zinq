//! Primitive Gaussian type definitions.

const std = @import("std");

const array_functions = @import("array_functions.zig");
const classical_particle = @import("classical_particle.zig");
const integral_functions = @import("integral_functions.zig");
const math_functions = @import("math_functions.zig");

const ClassicalParticle = classical_particle.ClassicalParticle;

const dfact = math_functions.dfact;
const boys = integral_functions.boys;
const powi = math_functions.powi;
const sum = array_functions.sum;

/// Primitive Gaussian type.
pub fn PrimitiveGaussian(comptime T: type) type {
    return struct {
        center: [3]T, angular: [3]T, alpha: T,

        /// Compute the coulomb integral between four primitive Gaussians.
        pub fn coulomb(self: PrimitiveGaussian(T), other1: PrimitiveGaussian(T), other2: PrimitiveGaussian(T), other3: PrimitiveGaussian(T)) T {
            const p = self.alpha + other1.alpha; const q = other2.alpha + other3.alpha; const beta = p * q / (p + q); var e: T = 0;

            const RP: [3]T = .{
                (self.alpha * self.center[0] + other1.alpha * other1.center[0]) / (self.alpha + other1.alpha),
                (self.alpha * self.center[1] + other1.alpha * other1.center[1]) / (self.alpha + other1.alpha),
                (self.alpha * self.center[2] + other1.alpha * other1.center[2]) / (self.alpha + other1.alpha)
            };

            const RQ: [3]T = .{
                (other2.alpha * other2.center[0] + other3.alpha * other3.center[0]) / (other2.alpha + other3.alpha),
                (other2.alpha * other2.center[1] + other3.alpha * other3.center[1]) / (other2.alpha + other3.alpha),
                (other2.alpha * other2.center[2] + other3.alpha * other3.center[2]) / (other2.alpha + other3.alpha)
            };

            const RPQ: [3]T = .{RP[0] - RQ[0], RP[1] - RQ[1], RP[2] - RQ[2]};

            for (0..@as(usize, @intFromFloat(self.angular[0] + other1.angular[0])) + 1) |t| {

                const Eij = hermite_coef(.{self.angular[0], other1.angular[0]}, .{self.alpha, other1.alpha}, self.center[0] - other1.center[0], @floatFromInt(t));

                for (0..@as(usize, @intFromFloat(self.angular[1] + other1.angular[1])) + 1) |u| {

                    const Ekl = hermite_coef(.{self.angular[1], other1.angular[1]}, .{self.alpha, other1.alpha}, self.center[1] - other1.center[1], @floatFromInt(u));

                    for (0..@as(usize, @intFromFloat(self.angular[2] + other1.angular[2])) + 1) |v| {

                        const Emn = hermite_coef(.{self.angular[2], other1.angular[2]}, .{self.alpha, other1.alpha}, self.center[2] - other1.center[2], @floatFromInt(v));

                        for (0..@as(usize, @intFromFloat(other2.angular[0] + other3.angular[0])) + 1) |tau| {

                            const Eop = hermite_coef(.{other2.angular[0], other3.angular[0]}, .{other2.alpha, other3.alpha}, other2.center[0] - other3.center[0], @floatFromInt(tau));

                            for (0..@as(usize, @intFromFloat(other2.angular[1] + other3.angular[1])) + 1) |nu| {

                                const Eqr = hermite_coef(.{other2.angular[1], other3.angular[1]}, .{other2.alpha, other3.alpha}, other2.center[1] - other3.center[1], @floatFromInt(nu));

                                for (0..@as(usize, @intFromFloat(other2.angular[2] + other3.angular[2])) + 1) |phi| {

                                    const Est = hermite_coef(.{other2.angular[2], other3.angular[2]}, .{other2.alpha, other3.alpha}, other2.center[2] - other3.center[2], @floatFromInt(phi));

                                    const sign: T = if ((tau + nu + phi) % 2 == 0) 1 else -1;

                                    e += sign * Eij * Ekl * Emn * Eop * Eqr * Est * hermite_integral([3]T{@floatFromInt(t + tau), @floatFromInt(u + nu), @floatFromInt(v + phi)}, RPQ, beta, 0);
                                }
                            }
                        }
                    }
                }
            }

            return 2 * std.math.pow(T, std.math.pi, 2.5) / (p * q * std.math.sqrt(p + q)) * e;
        }

        /// Compute the Hermite Gaussian coefficients in one dimension.
        pub fn hermite_coef(ij: [2]T, ab: [2]T, Q: T, t: T) T {
            const p = ab[0] + ab[1]; const q = ab[0] * ab[1] / p;

            if (ij[0] == 0 and ij[1] == 0 and t == 0) {return std.math.exp(-q * Q * Q);}

            else if (ij[0] > 0) {

                const E1 = hermite_coef(.{ij[0] - 1, ij[1]}, ab, Q, t - 1);
                const E2 = hermite_coef(.{ij[0] - 1, ij[1]}, ab, Q, t    );
                const E3 = hermite_coef(.{ij[0] - 1, ij[1]}, ab, Q, t + 1);

                return (1 / (2 * p)) * E1 - (q * Q / ab[0]) * E2 + (t + 1) * E3;
            }

            else if (ij[1] > 0) {

                const E1 = hermite_coef(.{ij[0], ij[1] - 1}, ab, Q, t - 1);
                const E2 = hermite_coef(.{ij[0], ij[1] - 1}, ab, Q, t    );
                const E3 = hermite_coef(.{ij[0], ij[1] - 1}, ab, Q, t + 1);

                return (1 / (2 * p)) * E1 + (q * Q / ab[1]) * E2 + (t + 1) * E3;
            }

            return 0;
        }

        /// Compute the Hermite integrals.
        pub fn hermite_integral(tuv: [3]T, RPC: [3]T, p: T, n: usize) T {
            if (tuv[0] == 0 and tuv[1] == 0 and tuv[2] == 0) {
                return powi(-2 * p, n) * boys(n, p * (RPC[0] * RPC[0] + RPC[1] * RPC[1] + RPC[2] * RPC[2]));
            }

            else if (tuv[0] > 0){
                return (tuv[0] - 1) * hermite_integral(.{tuv[0] - 2, tuv[1], tuv[2]}, RPC, p, n + 1) + RPC[0] * hermite_integral(.{tuv[0] - 1, tuv[1], tuv[2]}, RPC, p, n + 1);
            }

            else if (tuv[1] > 0) {
                return (tuv[1] - 1) * hermite_integral(.{tuv[0], tuv[1] - 2, tuv[2]}, RPC, p, n + 1) + RPC[1] * hermite_integral(.{tuv[0], tuv[1] - 1, tuv[2]}, RPC, p, n + 1);
            }

            else if (tuv[2] > 0) {
                return (tuv[2] - 1) * hermite_integral(.{tuv[0], tuv[1], tuv[2] - 2}, RPC, p, n + 1) + RPC[2] * hermite_integral(.{tuv[0], tuv[1], tuv[2] - 1}, RPC, p, n + 1);
            }

            return 0;
        }

        /// Compute the kinetic integral between two primitive Gaussians.
        pub fn kinetic(self: PrimitiveGaussian(T), other: PrimitiveGaussian(T)) T {
            const pgpl = PrimitiveGaussian(T){.center = other.center, .angular = .{other.angular[0] + 2, other.angular[1], other.angular[2]}, .alpha = other.alpha};
            const pgpm = PrimitiveGaussian(T){.center = other.center, .angular = .{other.angular[0], other.angular[1] + 2, other.angular[2]}, .alpha = other.alpha};
            const pgpn = PrimitiveGaussian(T){.center = other.center, .angular = .{other.angular[0], other.angular[1], other.angular[2] + 2}, .alpha = other.alpha};
            const pgml = PrimitiveGaussian(T){.center = other.center, .angular = .{other.angular[0] - 2, other.angular[1], other.angular[2]}, .alpha = other.alpha};
            const pgmm = PrimitiveGaussian(T){.center = other.center, .angular = .{other.angular[0], other.angular[1] - 2, other.angular[2]}, .alpha = other.alpha};
            const pgmn = PrimitiveGaussian(T){.center = other.center, .angular = .{other.angular[0], other.angular[1], other.angular[2] - 2}, .alpha = other.alpha};

            const T0 = other.alpha * (2 * sum(T, &other.angular) + 3) * self.overlap(other);

            const T1 = -2 * std.math.pow(T, other.alpha, 2.0) * (self.overlap(pgpl) + self.overlap(pgpm) + self.overlap(pgpn));

            const T2_1 = if (other.angular[0] > 1) other.angular[0] * (other.angular[0] - 1) * self.overlap(pgml) else 0;
            const T2_2 = if (other.angular[1] > 1) other.angular[1] * (other.angular[1] - 1) * self.overlap(pgmm) else 0;
            const T2_3 = if (other.angular[2] > 1) other.angular[2] * (other.angular[2] - 1) * self.overlap(pgmn) else 0;

            return T0 + T1 - 0.5 * (T2_1 + T2_2 + T2_3);
        }

        /// Calculate the norm of the primitive gaussian.
        pub fn norm(self: PrimitiveGaussian(T)) T {
            const Nij = dfact(2 * self.angular[0] - 1) * std.math.sqrt(0.5 * std.math.pi / self.alpha) / std.math.pow(T, 4 * self.alpha, self.angular[0]);
            const Nkl = dfact(2 * self.angular[1] - 1) * std.math.sqrt(0.5 * std.math.pi / self.alpha) / std.math.pow(T, 4 * self.alpha, self.angular[1]);
            const Nmn = dfact(2 * self.angular[2] - 1) * std.math.sqrt(0.5 * std.math.pi / self.alpha) / std.math.pow(T, 4 * self.alpha, self.angular[2]);

            return std.math.sqrt(Nij * Nkl * Nmn);
        }

        /// Compute the nuclear integral between two primitive Gaussians.
        pub fn nuclear(self: PrimitiveGaussian(T), other: PrimitiveGaussian(T), system: ClassicalParticle(T)) T {
            var n: T = 0;

            for (0..system.atoms.?.len) |i| {

                const RPC: [3]T = .{
                    (self.alpha * self.center[0] + other.alpha * other.center[0]) / (self.alpha + other.alpha) - system.position.at(3 * i + 0),
                    (self.alpha * self.center[1] + other.alpha * other.center[1]) / (self.alpha + other.alpha) - system.position.at(3 * i + 1),
                    (self.alpha * self.center[2] + other.alpha * other.center[2]) / (self.alpha + other.alpha) - system.position.at(3 * i + 2)
                };

                for (0..@as(usize, @intFromFloat(self.angular[0] + other.angular[0])) + 1) |t| {

                    const Eij = hermite_coef(.{self.angular[0], other.angular[0]}, .{self.alpha, other.alpha}, self.center[0] - other.center[0], @floatFromInt(t));

                    for (0..@as(usize, @intFromFloat(self.angular[1] + other.angular[1])) + 1) |u| {

                        const Ekl = hermite_coef(.{self.angular[1], other.angular[1]}, .{self.alpha, other.alpha}, self.center[1] - other.center[1], @floatFromInt(u));

                        for (0..@as(usize, @intFromFloat(self.angular[2] + other.angular[2])) + 1) |v| {

                            const Emn = hermite_coef(.{self.angular[2], other.angular[2]}, .{self.alpha, other.alpha}, self.center[2] - other.center[2], @floatFromInt(v));

                            n -= @as(T, @floatFromInt(system.atoms.?[i])) * Eij * Ekl * Emn * hermite_integral([3]T{@floatFromInt(t), @floatFromInt(u), @floatFromInt(v)}, RPC, self.alpha + other.alpha, 0);
                        }
                    }
                }
            }

            return 2 * std.math.pi / (self.alpha + other.alpha) * n;
        }

        /// Compute the overlap integral between two primitive Gaussians.
        pub fn overlap(self: PrimitiveGaussian(T), other: PrimitiveGaussian(T)) T {
            const Sij = hermite_coef(.{self.angular[0], other.angular[0]}, .{self.alpha, other.alpha}, self.center[0] - other.center[0], 0);
            const Skl = hermite_coef(.{self.angular[1], other.angular[1]}, .{self.alpha, other.alpha}, self.center[1] - other.center[1], 0);
            const Smn = hermite_coef(.{self.angular[2], other.angular[2]}, .{self.alpha, other.alpha}, self.center[2] - other.center[2], 0);

            return Sij * Skl * Smn * std.math.pow(T, std.math.pi / (self.alpha + other.alpha), 1.5);
        }
    };
}
