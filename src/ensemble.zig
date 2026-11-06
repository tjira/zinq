//! Generates a representation of a classical trajectory ensemble with positions, momenta, and active states.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;

const Matrix = @import("tensor.zig").Matrix;
const Vector = @import("tensor.zig").Vector;

const syevd = @import("linear_algebra.zig").syevd;

const AU2K = @import("constants.zig").AU2K;

/// Generates a representation of a classical trajectory ensemble with positions, momenta, and active states.
pub fn Ensemble(comptime T: type) type {
    return struct {
        r: Matrix(T),
        p: Matrix(T),
        a: Matrix(T),
        m: []const T,

        s: Vector(usize),

        ndof: usize,

        /// Initializes the ensemble trajectories with allocated memory for positions, momenta, and forces.
        pub fn init(ndim: usize, ntraj: usize, mass: []const T, gpa: Allocator) !@This() {
            var r = try Matrix(T).init(ntraj, ndim, gpa);
            errdefer r.deinit(gpa);

            var p = try Matrix(T).init(ntraj, ndim, gpa);
            errdefer p.deinit(gpa);

            var a = try Matrix(T).init(ntraj, ndim, gpa);
            errdefer a.deinit(gpa);

            var s = try Vector(usize).init(ntraj, gpa);
            errdefer s.deinit(gpa);

            const m = try gpa.alloc(T, mass.len);
            errdefer gpa.free(m);

            @memcpy(m, mass);

            return .{ .r = r, .p = p, .a = a, .s = s, .m = m, .ndof = ndim };
        }

        /// Deallocates the positions, momenta, and force vectors of the trajectory ensemble.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.r.deinit(gpa);
            self.p.deinit(gpa);
            self.a.deinit(gpa);
            self.s.deinit(gpa);

            gpa.free(self.m);
        }

        /// Calculates the average classical kinetic energy of the ensemble trajectories.
        pub fn ekin(self: @This()) T {
            var sum: T = 0;

            for (0..self.p.nrow()) |i| for (0..self.p.ncol()) |j| {
                sum += self.p.at(i, j) * self.p.at(i, j) / (2 * self.m[j]);
            };

            return sum / @as(T, @floatFromInt(self.p.nrow()));
        }

        /// Computes the mean momentum vector averaged over all trajectories in the ensemble.
        pub fn mom(self: @This(), gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(self.p.ncol(), gpa);

            for (0..self.p.nrow()) |i| for (0..self.p.ncol()) |j| {
                value.ptr(j).* += self.p.at(i, j);
            };

            value.divs(@floatFromInt(self.p.nrow()));

            return value;
        }

        /// Computes the population fraction of trajectories occupying each electronic state.
        pub fn pop(self: @This(), nstate: usize, coefs: ?Matrix(Complex(T)), gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(nstate, gpa);

            if (coefs) |coef| {
                for (0..coef.nrow()) |i| for (0..nstate) |k| {
                    value.ptr(k).* += coef.at(i, k).squaredMagnitude();
                };

                value.divs(@floatFromInt(self.r.nrow()));
            }

            if (coefs == null) {
                for (0..self.s.length()) |i| {
                    value.ptr(self.s.at(i)).* += 1;
                }

                value.divs(@floatFromInt(self.s.length()));
            }

            return value;
        }

        /// Computes the mean position vector averaged over all trajectories in the ensemble.
        pub fn pos(self: @This(), gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(self.r.ncol(), gpa);

            for (0..self.r.nrow()) |i| for (0..self.r.ncol()) |j| {
                value.ptr(j).* += self.r.at(i, j);
            };

            value.divs(@floatFromInt(self.r.nrow()));

            return value;
        }

        /// Projects out overall angular momentum and rotational velocities from the ensemble.
        pub fn removeRotation(self: *@This()) usize {
            const natom = self.m.len / 3;

            if (natom < 2) return 0;

            var m_tot: T = 0;

            for (0..natom) |k| {
                m_tot += self.m[k * 3];
            }

            if (m_tot <= 0) return 0;

            var axes_removed: usize = 0;

            for (0..self.r.nrow()) |i| {
                const r_row = self.r.rowSlice(i);
                const p_row = self.p.rowSlice(i);

                const r_cm = computeCenterOfMass(T, r_row, self.m, m_tot);

                const l_and_i = computeInertiaTensor(T, r_row, p_row, self.m, r_cm);

                const ang = computeAngularVelocity(T, l_and_i.L, l_and_i.I);

                if (i == 0 or axes_removed == 0) axes_removed = ang.axes_removed;

                removeRotationalVelocity(T, p_row, r_row, self.m, r_cm, ang.omega);
            }

            self.ndof -|= axes_removed;

            return axes_removed;
        }

        /// Removes center-of-mass momentum drift from all trajectories in the ensemble.
        pub fn removeTranslation(self: *@This()) void {
            const natom, var m_tot: T = .{ self.m.len / 3, 0 };

            for (0..natom) |k| m_tot += self.m[k * 3];

            if (m_tot <= 0) return;

            for (0..self.r.nrow()) |i| {
                var p_cm = [3]T{ 0, 0, 0 };

                for (0..natom) |k| for (0..3) |c| {
                    p_cm[c] += self.p.at(i, k * 3 + c);
                };

                for (0..3) |c| {
                    p_cm[c] /= m_tot;
                }

                for (0..natom) |k| for (0..3) |c| {
                    self.p.ptr(i, k * 3 + c).* -= self.m[k * 3 + c] * p_cm[c];
                };
            }

            self.ndof -|= 3;
        }

        /// Rescales ensemble momenta to match the target kinetic temperature across active degrees of freedom.
        pub fn rescaleTemperature(self: *@This(), temp_k: T) void {
            if (self.ndof == 0) return;

            const target_ekin = 0.5 * @as(T, @floatFromInt(self.ndof)) * (temp_k / AU2K);

            for (0..self.p.nrow()) |i| {
                var cur_ekin: T = 0;

                for (0..self.p.ncol()) |j| {
                    cur_ekin += (self.p.at(i, j) * self.p.at(i, j)) / (2 * self.m[j]);
                }

                if (cur_ekin > 1e-12) {
                    const factor = @sqrt(target_ekin / cur_ekin);

                    for (0..self.p.ncol()) |j| {
                        self.p.ptr(i, j).* *= factor;
                    }
                }
            }
        }

        /// Samples initial positions and momenta from a Wigner-like Gaussian distribution or reference geometry.
        pub fn setGaussian(self: *@This(), ic: anytype, default_pos: ?[]const T) void {
            var split_mix = std.Random.SplitMix64.init(ic.seed);
            var rng = std.Random.DefaultPrng.init(split_mix.next());

            for (0..self.s.length()) |i| {
                self.s.ptr(i).* = ic.state;
            }

            const gamma = if (ic.distribution) |dist| switch (dist) {
                .wigner => |wig| wig.gamma,
                inline else => null,
            } else null;

            const temperature = if (ic.distribution) |dist| switch (dist) {
                .maxwell_boltzmann => |mb| mb.temperature,
                inline else => null,
            } else null;

            const random = rng.random();

            sampleGaussianCoordinates(T, &self.r, ic.position, default_pos, gamma, random);

            sampleMomenta(T, &self.p, ic.momentum, self.m, temperature, gamma, random);
        }

        /// Computes the ensemble-averaged electronic state coefficient magnitudes for nonadiabatic dynamics.
        pub fn stateCoef(self: @This(), nstate: usize, coefs: Matrix(Complex(T)), gpa: Allocator) !Vector(T) {
            var value = try Vector(T).initZero(nstate, gpa);

            for (0..coefs.nrow()) |i| for (0..nstate) |k| {
                value.ptr(k).* += coefs.at(i, k).squaredMagnitude();
            };

            value.divs(@floatFromInt(self.r.nrow()));

            return value;
        }

        /// Calculates the instantaneous kinetic temperature of the ensemble in Kelvin.
        pub fn temp(self: @This()) T {
            if (self.ndof == 0) return 0;

            return (2 * self.ekin() / @as(T, @floatFromInt(self.ndof))) * AU2K;
        }

        /// Diagonalizes the inertia tensor and projects angular momentum to obtain angular velocity.
        fn computeAngularVelocity(comptime U: type, L: [3]U, I_tensor: [9]U) struct { omega: [3]U, axes_removed: usize } {
            var ev: [3]U, var evec: [9]U = .{ undefined, undefined };

            var I_copy = I_tensor;

            var ev_vec = Vector(U).fromSlice(&ev);

            var evec_mat = Matrix(U).fromSlice(3, 3, &evec);

            syevd(U, &ev_vec, &evec_mat, Matrix(U).fromSlice(3, 3, &I_copy)) catch {
                return .{ .omega = .{ 0, 0, 0 }, .axes_removed = 0 };
            };

            var omega, var cur_axes: usize = .{ [3]U{ 0, 0, 0 }, 0 };

            for (0..3) |axis| if (ev[axis] > 1e-8) {
                cur_axes += 1;

                const l_proj = evec[axis] * L[0] + evec[3 + axis] * L[1] + evec[6 + axis] * L[2];

                for (0..3) |c| {
                    omega[c] += (l_proj / ev[axis]) * evec[c * 3 + axis];
                }
            };

            return .{ .omega = omega, .axes_removed = cur_axes };
        }

        /// Computes the mass-weighted center of mass coordinate vector.
        fn computeCenterOfMass(comptime U: type, r: []const U, m: []const U, m_tot: U) [3]U {
            const natom, var r_cm = .{ m.len / 3, [3]U{ 0, 0, 0 } };

            for (0..natom) |k| for (0..3) |c| {
                r_cm[c] += m[k * 3 + c] * r[k * 3 + c];
            };

            for (0..3) |c| r_cm[c] /= m_tot;

            return r_cm;
        }

        /// Calculates the angular momentum vector and inertia tensor relative to the center of mass.
        fn computeInertiaTensor(comptime U: type, r: []const U, p: []const U, m: []const U, r_cm: [3]U) struct { L: [3]U, I: [9]U } {
            const natom, var L = .{ m.len / 3, [3]U{ 0, 0, 0 } };

            var I_mat = [9]U{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            for (0..natom) |k| {
                const mk = m[k * 3];

                const rx = r[k * 3 + 0] - r_cm[0];
                const ry = r[k * 3 + 1] - r_cm[1];
                const rz = r[k * 3 + 2] - r_cm[2];

                const px = p[k * 3 + 0];
                const py = p[k * 3 + 1];
                const pz = p[k * 3 + 2];

                L[0] += ry * pz - rz * py;
                L[1] += rz * px - rx * pz;
                L[2] += rx * py - ry * px;

                I_mat[0] += mk * (ry * ry + rz * rz);
                I_mat[4] += mk * (rx * rx + rz * rz);
                I_mat[8] += mk * (rx * rx + ry * ry);

                I_mat[1] -= mk * rx * ry;
                I_mat[2] -= mk * rx * rz;
                I_mat[5] -= mk * ry * rz;
            }

            I_mat[3] = I_mat[1];
            I_mat[6] = I_mat[2];
            I_mat[7] = I_mat[5];

            return .{ .L = L, .I = I_mat };
        }

        /// Subtracts rotational linear momentum from trajectory coordinates.
        fn removeRotationalVelocity(comptime U: type, p: []U, r: []const U, m: []const U, r_cm: [3]U, omega: [3]U) void {
            const natom = m.len / 3;

            for (0..natom) |k| {
                const mk = m[k * 3];

                const rx = r[k * 3 + 0] - r_cm[0];
                const ry = r[k * 3 + 1] - r_cm[1];
                const rz = r[k * 3 + 2] - r_cm[2];

                p[k * 3 + 0] -= mk * (omega[1] * rz - omega[2] * ry);
                p[k * 3 + 1] -= mk * (omega[2] * rx - omega[0] * rz);
                p[k * 3 + 2] -= mk * (omega[0] * ry - omega[1] * rx);
            }
        }

        /// Samples ensemble Cartesian coordinates from reference geometry and optional harmonic widths.
        fn sampleGaussianCoordinates(comptime U: type, r: *Matrix(U), pos_opt: ?[]const f64, def_pos: ?[]const U, gamma: ?[]const f64, random: anytype) void {
            for (0..r.nrow()) |i| for (0..r.ncol()) |j| {
                const r0 = if (pos_opt) |r_pos| @as(U, @floatCast(r_pos[j])) else def_pos.?[j];

                if (gamma) |gam| {
                    const stdev = 1.0 / std.math.sqrt(2.0 * @as(U, @floatCast(gam[j])));

                    r.ptr(i, j).* = r0 + stdev * random.floatNorm(U);
                }

                if (gamma == null) {
                    r.ptr(i, j).* = r0;
                }
            };
        }

        /// Samples ensemble Cartesian momenta from thermal or harmonic phase space distributions.
        fn sampleMomenta(comptime U: type, p: *Matrix(U), mom_opt: ?[]const f64, m: []const U, temp_opt: ?f64, gamma: ?[]const f64, random: anytype) void {
            for (0..p.nrow()) |i| for (0..p.ncol()) |j| {
                const p0 = if (mom_opt) |p_mom| @as(U, @floatCast(p_mom[j])) else 0;

                if (temp_opt) |t_val| {
                    const stdev = std.math.sqrt(m[j] * @as(U, @floatCast(t_val)) / AU2K);

                    p.ptr(i, j).* = p0 + stdev * random.floatNorm(U);

                    continue;
                }

                if (gamma) |gam| {
                    const stdev = std.math.sqrt(@as(U, @floatCast(gam[j])) / 2.0);

                    p.ptr(i, j).* = p0 + stdev * random.floatNorm(U);

                    continue;
                }

                p.ptr(i, j).* = p0;
            };
        }
    };
}
