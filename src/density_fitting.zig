//! Implements auxiliary basis metric inversion and 3-center integral transformations for density fitting.

const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

const dot = @import("linear_algebra.zig").dot;
const exportAuxIfBuiltin = @import("molecular_integrals.zig").exportAuxIfBuiltin;
const gemm = @import("linear_algebra.zig").gemm;
const gemv = @import("linear_algebra.zig").gemv;
const geru = @import("linear_algebra.zig").geru;
const syevd = @import("linear_algebra.zig").syevd;

/// Configuration flags selecting memory allocation for density fitting workspace buffers.
pub const RiFlags = struct {
    coulomb: bool = false,
    exchange: bool = false,
    generalized: bool = false,
    response: bool = false,
};

/// Preallocated scratch buffers and auxiliary metric matrices for density fitting evaluation.
pub fn RiWorkspace(comptime T: type) type {
    return struct {
        dC_occ: ?Matrix(T) = null,
        W_resp: ?Matrix(T) = null,

        J_x_vec: ?Vector(T) = null,

        C_occ: ?Matrix(T) = null,
        P_tot: ?Vector(T) = null,
        q_aux: ?Vector(T) = null,
        q_der: ?Vector(T) = null,
        v_vec: ?Vector(T) = null,
        W_der: ?Matrix(T) = null,
        W_exc: ?Matrix(T) = null,

        /// Allocates intermediate tensors and auxiliary vectors for density fitting calculations.
        pub fn init(nbf: usize, naux: usize, nocc: usize, flags: RiFlags, gpa: Allocator) !@This() {
            var self: @This() = .{};

            if (flags.coulomb) {
                self.q_aux = try Vector(T).initZero(naux, gpa);
                errdefer if (self.q_aux) |*v| v.deinit(gpa);

                if (flags.generalized) {
                    const nbf_sq = nbf * nbf;

                    self.P_tot = try Vector(T).initZero(nbf_sq, gpa);
                    errdefer if (self.P_tot) |*v| v.deinit(gpa);

                    self.v_vec = try Vector(T).initZero(nbf_sq, gpa);
                    errdefer if (self.v_vec) |*v| v.deinit(gpa);
                }
            }

            if (flags.exchange) {
                const nspin = if (flags.generalized) 2 * nbf else nbf;

                self.W_exc = try Matrix(T).initZero(nspin, nocc * naux, gpa);
                errdefer if (self.W_exc) |*m| m.deinit(gpa);

                self.C_occ = try Matrix(T).initZero(nspin, nocc, gpa);
                errdefer if (self.C_occ) |*m| m.deinit(gpa);
            }

            if (flags.response) {
                const nspin = if (flags.generalized) 2 * nbf else nbf;

                self.q_der = try Vector(T).initZero(naux, gpa);
                errdefer if (self.q_der) |*v| v.deinit(gpa);

                self.W_der = try Matrix(T).initZero(nspin, nocc * naux, gpa);
                errdefer if (self.W_der) |*m| m.deinit(gpa);

                self.W_resp = try Matrix(T).initZero(nspin, nocc * naux, gpa);
                errdefer if (self.W_resp) |*m| m.deinit(gpa);

                self.dC_occ = try Matrix(T).initZero(nspin, nocc, gpa);
                errdefer if (self.dC_occ) |*m| m.deinit(gpa);

                if (flags.generalized) {
                    const nbf_sq = nbf * nbf;

                    self.J_x_vec = try Vector(T).initZero(nbf_sq, gpa);
                    errdefer if (self.J_x_vec) |*v| v.deinit(gpa);
                }
            }

            return self;
        }

        /// Deallocates all scratch tensors and auxiliary matrices stored in the density fitting workspace.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            if (self.dC_occ) |*m| m.deinit(gpa);
            if (self.W_resp) |*m| m.deinit(gpa);

            if (self.J_x_vec) |*v| v.deinit(gpa);

            if (self.C_occ) |*m| m.deinit(gpa);
            if (self.P_tot) |*v| v.deinit(gpa);
            if (self.q_aux) |*v| v.deinit(gpa);
            if (self.q_der) |*v| v.deinit(gpa);
            if (self.v_vec) |*v| v.deinit(gpa);
            if (self.W_der) |*m| m.deinit(gpa);
            if (self.W_exc) |*m| m.deinit(gpa);
        }

        /// Contracts the 3-center density fitting tensor with the density matrix to form the Coulomb Fock contribution.
        pub fn contractCoulomb(self: *@This(), F: *Matrix(T), B: Tensor(T, 3), P: Matrix(T)) void {
            const nbf = B.shape[0];

            const B_mat = B.asMatrix();

            if (P.shape[0] != nbf) {
                for (0..nbf) |mu| {
                    for (0..nbf) |nu| {
                        self.P_tot.?.ptr(mu * nbf + nu).* = P.at(mu, nu) + P.at(mu + nbf, nu + nbf);
                    }
                }

                gemv(T, &self.q_aux.?, B_mat, self.P_tot.?, 1, 0, true);

                gemv(T, &self.v_vec.?, B_mat, self.q_aux.?, 1, 0, false);

                for (0..nbf) |mu| {
                    for (0..nbf) |nu| {
                        const j_val = self.v_vec.?.at(mu * nbf + nu);

                        F.ptr(mu, nu).* += j_val;

                        F.ptr(mu + nbf, nu + nbf).* += j_val;
                    }
                }
            }

            if (P.shape[0] == nbf) {
                gemv(T, &self.q_aux.?, B_mat, P.asVector(), 1, 0, true);

                var F_vec = F.asVector();

                gemv(T, &F_vec, B_mat, self.q_aux.?, 1, 1, false);
            }
        }

        /// Contracts 3-center fitting tensors with molecular orbital coefficients to form the exchange Fock matrix.
        pub fn contractExchange(self: *@This(), F: *Matrix(T), B: Tensor(T, 3), C: Matrix(T), exch_factor: T) void {
            self.updateExchange(B, C);

            const scale: T = if (C.shape[0] != B.shape[0]) -exch_factor else -2.0 * exch_factor;

            gemm(T, F, self.W_exc.?, self.W_exc.?, scale, 1, false, true);
        }

        /// Contracts fitting tensors with perturbed MO coefficients to evaluate exchange response.
        pub fn contractExchangeResponse(self: *@This(), V_x: *Matrix(T), B: Tensor(T, 3), dC_x: Matrix(T)) void {
            const nbf, const naux = .{ B.shape[0], B.shape[2] };

            const nspin, const nocc = .{ self.dC_occ.?.shape[0], self.dC_occ.?.shape[1] };

            for (0..nspin) |nu| for (0..nocc) |i| {
                self.dC_occ.?.ptr(nu, i).* = dC_x.at(nu, i);
            };

            if (nspin == nbf) {
                for (0..nbf) |mu| {
                    const B_mu = B.slice(mu).asMatrix();

                    var W_resp_mu = self.W_resp.?.rowMatrix(mu, nocc, naux);

                    gemm(T, &W_resp_mu, self.dC_occ.?, B_mu, 1, 0, true, false);
                }
            }

            if (nspin != nbf) {
                const dC_alpha = self.dC_occ.?.sliceRows(0, nbf);

                const dC_beta = self.dC_occ.?.sliceRows(nbf, nbf);

                for (0..nbf) |mu| {
                    const B_mu = B.slice(mu).asMatrix();

                    var W_resp_a = self.W_resp.?.rowMatrix(mu, nocc, naux);

                    gemm(T, &W_resp_a, dC_alpha, B_mu, 1, 0, true, false);

                    var W_resp_b = self.W_resp.?.rowMatrix(nbf + mu, nocc, naux);

                    gemm(T, &W_resp_b, dC_beta, B_mu, 1, 0, true, false);
                }
            }

            gemm(T, V_x, self.W_resp.?, self.W_exc.?, -1.0, 1.0, false, true);
            gemm(T, V_x, self.W_exc.?, self.W_resp.?, -1.0, 1.0, false, true);
        }

        /// Contracts fitting tensor derivatives with density and orbitals to form the 2e Fock response.
        pub fn contractFockDerivative(self: *@This(), F_x: *Matrix(T), B: Tensor(T, 3), dB: Tensor(T, 3), P: Matrix(T), C: Matrix(T)) void {
            const nbf, const naux = .{ B.shape[0], B.shape[2] };

            const nocc = self.C_occ.?.shape[1];

            const B_mat, const dB_mat = .{ B.asMatrix(), dB.asMatrix() };

            if (P.shape[0] != nbf) {
                for (0..nbf) |mu| for (0..nbf) |nu| {
                    self.P_tot.?.ptr(mu * nbf + nu).* = P.at(mu, nu) + P.at(mu + nbf, nu + nbf);
                };

                gemv(T, &self.q_aux.?, B_mat, self.P_tot.?, 1, 0, true);

                gemv(T, &self.q_der.?, dB_mat, self.P_tot.?, 1, 0, true);

                gemv(T, &self.J_x_vec.?, dB_mat, self.q_aux.?, 1, 0, false);

                gemv(T, &self.J_x_vec.?, B_mat, self.q_der.?, 1, 1, false);

                for (0..nbf) |mu| for (0..nbf) |nu| {
                    const j_val = self.J_x_vec.?.at(mu * nbf + nu);

                    F_x.ptr(mu, nu).* += j_val;

                    F_x.ptr(mu + nbf, nu + nbf).* += j_val;
                };

                for (0..2 * nbf) |nu| for (0..nocc) |i| {
                    self.C_occ.?.ptr(nu, i).* = C.at(nu, i);
                };

                const C_occ_alpha = self.C_occ.?.sliceRows(0, nbf);

                const C_occ_beta = self.C_occ.?.sliceRows(nbf, nbf);

                for (0..nbf) |mu| {
                    const dB_mu = dB.slice(mu).asMatrix();

                    var W_der_a = self.W_der.?.rowMatrix(mu, nocc, naux);

                    gemm(T, &W_der_a, C_occ_alpha, dB_mu, 1, 0, true, false);

                    var W_der_b = self.W_der.?.rowMatrix(nbf + mu, nocc, naux);

                    gemm(T, &W_der_b, C_occ_beta, dB_mu, 1, 0, true, false);
                }

                gemm(T, F_x, self.W_der.?, self.W_exc.?, -1.0, 1.0, false, true);
                gemm(T, F_x, self.W_exc.?, self.W_der.?, -1.0, 1.0, false, true);
            }

            if (P.shape[0] == nbf) {
                gemv(T, &self.q_aux.?, B_mat, P.asVector(), 1, 0, true);

                gemv(T, &self.q_der.?, dB_mat, P.asVector(), 1, 0, true);

                var F_vec = F_x.asVector();

                gemv(T, &F_vec, dB_mat, self.q_aux.?, 1, 1, false);

                gemv(T, &F_vec, B_mat, self.q_der.?, 1, 1, false);

                for (0..nbf) |nu| for (0..nocc) |i| {
                    self.C_occ.?.ptr(nu, i).* = C.at(nu, i);
                };

                for (0..nbf) |mu| {
                    const dB_mu = dB.slice(mu).asMatrix();

                    var W_der_mu = self.W_der.?.rowMatrix(mu, nocc, naux);

                    gemm(T, &W_der_mu, self.C_occ.?, dB_mu, 1, 0, true, false);
                }

                gemm(T, F_x, self.W_der.?, self.W_exc.?, -1.0, 1.0, false, true);
                gemm(T, F_x, self.W_exc.?, self.W_der.?, -1.0, 1.0, false, true);
            }
        }

        /// Contracts 3-center and 2-center derivative integrals with density and fitting coefficients for gradients.
        pub fn contractGradient(self: *@This(), G: *Matrix(T), ints: anytype, C: Matrix(T), P: Matrix(T), exch_factor: T, gpa: Allocator) !void {
            const nbf, const naux, const nat = .{ ints.B.?.shape[0], ints.J.?.nrow(), G.nrow() };

            const nbf_sq = nbf * nbf;

            var J_inv_half = try invertMetric(T, ints.J.?, gpa);
            defer J_inv_half.deinit(gpa);

            const B_mat = ints.B.?.asMatrix();

            var w_aux = try Vector(T).initZero(naux, gpa);
            defer w_aux.deinit(gpa);

            const P_vec = if (P.shape[0] == nbf) P.asVector() else blk: {
                for (0..nbf) |mu| for (0..nbf) |nu| {
                    self.P_tot.?.ptr(mu * nbf + nu).* = P.at(mu, nu) + P.at(mu + nbf, nu + nbf);
                };

                break :blk self.P_tot.?;
            };

            gemv(T, &self.q_aux.?, B_mat, P_vec, 1, 0, true);

            gemv(T, &w_aux, J_inv_half, self.q_aux.?, 1, 0, false);

            var Y = try Matrix(T).initZero(nbf_sq, naux, gpa);
            defer Y.deinit(gpa);

            geru(T, &Y, P_vec, w_aux, 1.0);

            var Z = try Matrix(T).initZero(naux, naux, gpa);
            defer Z.deinit(gpa);

            geru(T, &Z, w_aux, w_aux, -0.5);

            if (exch_factor > 0) {
                self.updateExchange(ints.B.?, C);

                const nocc = self.C_occ.?.shape[1];

                var B_mo = try Matrix(T).initZero(nocc, nocc * naux, gpa);
                defer B_mo.deinit(gpa);

                gemm(T, &B_mo, self.C_occ.?, self.W_exc.?, 1, 0, true, false);

                const B_mo_mat = B_mo.reshape(nocc * nocc, naux);

                var w_mo = try Matrix(T).initZero(nocc * nocc, naux, gpa);
                defer w_mo.deinit(gpa);

                gemm(T, &w_mo, B_mo_mat, J_inv_half, 1, 0, false, false);

                const z_factor: T = if (P.shape[0] != nbf) 0.5 * exch_factor else 2.0 * exch_factor;

                gemm(T, &Z, w_mo, w_mo, z_factor, 1.0, true, false);

                var D_exc = try Matrix(T).initZero(nbf_sq, naux, gpa);
                defer D_exc.deinit(gpa);

                const y_factor: T = if (P.shape[0] != nbf) exch_factor else 4.0 * exch_factor;

                var T_i = try Matrix(T).initZero(nbf, naux, gpa);
                defer T_i.deinit(gpa);

                if (P.shape[0] == nbf) {
                    for (0..nocc) |i| {
                        const B_i = B_mo_mat.sliceRows(i * nocc, nocc);

                        gemm(T, &T_i, self.C_occ.?, B_i, 1, 0, false, false);

                        for (0..nbf) |mu| {
                            const c_val = y_factor * self.C_occ.?.at(mu, i);

                            for (0..nbf) |nu| for (0..naux) |Q| {
                                D_exc.ptr(mu * nbf + nu, Q).* += c_val * T_i.at(nu, Q);
                            };
                        }
                    }
                }

                if (P.shape[0] != nbf) {
                    var T_b_i = try Matrix(T).initZero(nbf, naux, gpa);
                    defer T_b_i.deinit(gpa);

                    const C_occ_alpha = self.C_occ.?.sliceRows(0, nbf);

                    const C_occ_beta = self.C_occ.?.sliceRows(nbf, nbf);

                    for (0..nocc) |i| {
                        const B_i = B_mo_mat.sliceRows(i * nocc, nocc);

                        gemm(T, &T_i, C_occ_alpha, B_i, 1, 0, false, false);

                        gemm(T, &T_b_i, C_occ_beta, B_i, 1, 0, false, false);

                        for (0..nbf) |mu| {
                            const ca_val = y_factor * C_occ_alpha.at(mu, i);

                            const cb_val = y_factor * C_occ_beta.at(mu, i);

                            for (0..nbf) |nu| for (0..naux) |Q| {
                                D_exc.ptr(mu * nbf + nu, Q).* += ca_val * T_i.at(nu, Q) + cb_val * T_b_i.at(nu, Q);
                            };
                        }
                    }
                }

                gemm(T, &Y, D_exc, J_inv_half, -1.0, 1.0, false, false);
            }

            const Y_vec = Y.asVector();
            const Z_vec = Z.asVector();

            for (0..nat) |A_atom| for (0..3) |alpha| {
                const coord = 3 * A_atom + alpha;

                const dI_c = ints.dI.?.slice(coord).asVector();
                const dJ_c = ints.dJ.?.slice(coord).asVector();

                const term1 = dot(T, Y_vec, dI_c);
                const term2 = dot(T, Z_vec, dJ_c);

                G.ptr(A_atom, alpha).* += term1 + term2;
            };
        }

        /// Transforms 3-center fitting tensors to half-MO basis scratch buffers.
        pub fn updateExchange(self: *@This(), B: Tensor(T, 3), C: Matrix(T)) void {
            const nbf, const naux = .{ B.shape[0], B.shape[2] };

            const nspin, const nocc = .{ self.C_occ.?.shape[0], self.C_occ.?.shape[1] };

            if (nspin == nbf) {
                for (0..nbf) |nu| {
                    for (0..nocc) |i| {
                        self.C_occ.?.ptr(nu, i).* = C.at(nu, i);
                    }
                }

                for (0..nbf) |mu| {
                    const B_mu = B.slice(mu).asMatrix();

                    var W_mu = self.W_exc.?.rowMatrix(mu, nocc, naux);

                    gemm(T, &W_mu, self.C_occ.?, B_mu, 1, 0, true, false);
                }
            }

            if (nspin != nbf) {
                for (0..2 * nbf) |nu| {
                    for (0..nocc) |i| {
                        self.C_occ.?.ptr(nu, i).* = C.at(nu, i);
                    }
                }

                const C_occ_alpha = self.C_occ.?.sliceRows(0, nbf);

                const C_occ_beta = self.C_occ.?.sliceRows(nbf, nbf);

                for (0..nbf) |mu| {
                    const B_mu = B.slice(mu).asMatrix();

                    var W_a_mu = self.W_exc.?.rowMatrix(mu, nocc, naux);

                    gemm(T, &W_a_mu, C_occ_alpha, B_mu, 1, 0, true, false);

                    var W_b_mu = self.W_exc.?.rowMatrix(nbf + mu, nocc, naux);

                    gemm(T, &W_b_mu, C_occ_beta, B_mu, 1, 0, true, false);
                }
            }
        }
    };
}

/// Evaluates and assigns auxiliary Coulomb integrals and derivatives for density fitting.
pub fn buildFittingIntegrals(comptime T: type, io: std.Io, ints: anytype, aux_basis: []const u8, needs_derivs: bool, nthreads: usize, gpa: Allocator) !void {
    const aux_path = try exportAuxIfBuiltin(io, aux_basis, gpa);

    defer if (std.mem.startsWith(u8, aux_basis, "builtin:")) {
        std.Io.Dir.cwd().deleteFile(io, aux_path) catch {};
    };

    if (ints.J == null) {
        ints.J = try ints.sys.coulomb2c(aux_path, nthreads, gpa);
    }

    if (ints.I == null and (ints.B == null or needs_derivs)) {
        ints.I = try ints.sys.coulomb3c(aux_path, nthreads, gpa);
    }

    if (ints.B == null) {
        var J_inv_half = try invertMetric(T, ints.J.?, gpa);
        defer J_inv_half.deinit(gpa);

        ints.B = try Tensor(T, 3).initZero(ints.I.?.shape, gpa);
        contractMetric(T, &ints.B.?, ints.I.?, J_inv_half);
    }

    if (needs_derivs) {
        if (ints.dJ == null) {
            ints.dJ = try ints.sys.coulomb2cD1(aux_path, nthreads, gpa);
        }

        if (ints.dI == null) {
            ints.dI = try ints.sys.coulomb3cD1(aux_path, nthreads, gpa);
        }
    }
}

/// Evaluates and contracts auxiliary Coulomb integrals to form the atomic orbital fitting tensor.
pub fn buildFittingTensor(comptime T: type, io: std.Io, sys: anytype, aux_basis: []const u8, nthreads: usize, gpa: Allocator) !Tensor(T, 3) {
    const aux_path = try exportAuxIfBuiltin(io, aux_basis, gpa);

    defer if (std.mem.startsWith(u8, aux_basis, "builtin:")) {
        std.Io.Dir.cwd().deleteFile(io, aux_path) catch {};
    };

    var J_aux = try sys.coulomb2c(aux_path, nthreads, gpa);
    defer J_aux.deinit(gpa);

    var J_inv_half = try invertMetric(T, J_aux, gpa);
    defer J_inv_half.deinit(gpa);

    var I_3c = try sys.coulomb3c(aux_path, nthreads, gpa);
    defer I_3c.deinit(gpa);

    var B_ao = try Tensor(T, 3).initZero(I_3c.shape, gpa);
    errdefer B_ao.deinit(gpa);

    contractMetric(T, &B_ao, I_3c, J_inv_half);

    return B_ao;
}

/// Evaluates Daleckii-Krein kernel matrix for the matrix inverse square root derivative.
pub fn buildMetricDerivKernel(comptime T: type, W: Vector(T), gpa: Allocator) !Matrix(T) {
    const naux = W.shape[0];

    var M = try Matrix(T).init(naux, naux, gpa);
    errdefer M.deinit(gpa);

    for (0..naux) |p| {
        const wp = W.at(p);

        const sqrt_p = if (wp > 1e-14) @sqrt(wp) else 0;

        for (0..naux) |q| {
            const wq = W.at(q);

            const sqrt_q = if (wq > 1e-14) @sqrt(wq) else 0;

            if (sqrt_p > 1e-14 and sqrt_q > 1e-14) {
                M.ptr(p, q).* = -1.0 / (sqrt_p * sqrt_q * (sqrt_p + sqrt_q));
            }

            if (sqrt_p <= 1e-14 or sqrt_q <= 1e-14) {
                M.ptr(p, q).* = 0;
            }
        }
    }

    return M;
}

/// Contracts 3-center Coulomb integrals with the inverse square root of the auxiliary metric using BLAS gemm.
pub fn contractMetric(comptime T: type, B_ao: *Tensor(T, 3), I_3c: Tensor(T, 3), J_inv_half: Matrix(T)) void {
    const I_mat = I_3c.asMatrix();

    var B_mat = B_ao.asMatrix();

    gemm(T, &B_mat, I_mat, J_inv_half, 1, 0, false, false);
}

/// Contracts 3-center integral and metric derivatives to evaluate fitting tensor gradient.
pub fn evaluateFittingTensorDeriv(comptime T: type, dB_ao: *Tensor(T, 3), dI_3c: Tensor(T, 3), I_3c: Tensor(T, 3), J_inv_half: Matrix(T), dJ_inv_half: Matrix(T)) void {
    const I_mat, const dI_mat = .{ I_3c.asMatrix(), dI_3c.asMatrix() };

    var dB_mat = dB_ao.asMatrix();

    gemm(T, &dB_mat, dI_mat, J_inv_half, 1, 0, false, false);
    gemm(T, &dB_mat, I_mat, dJ_inv_half, 1, 1, false, false);
}

/// Evaluates the nuclear derivative of the inverse square root metric via Daleckii-Krein.
pub fn evaluateMetricDeriv(comptime T: type, dJ_inv_half: *Matrix(T), dJ: Matrix(T), U: Matrix(T), M: Matrix(T), gpa: Allocator) !void {
    const naux = dJ.nrow();

    var tmp1 = try Matrix(T).init(naux, naux, gpa);
    defer tmp1.deinit(gpa);

    var tmp2 = try Matrix(T).init(naux, naux, gpa);
    defer tmp2.deinit(gpa);

    gemm(T, &tmp1, U, dJ, 1, 0, true, false);

    gemm(T, &tmp2, tmp1, U, 1, 0, false, false);

    for (0..naux) |p| for (0..naux) |q| {
        tmp1.ptr(p, q).* = M.at(p, q) * tmp2.at(p, q);
    };

    gemm(T, &tmp2, U, tmp1, 1, 0, false, false);

    gemm(T, dJ_inv_half, tmp2, U, 1, 0, false, true);
}

/// Inverts the auxiliary Coulomb metric matrix and evaluates its symmetric inverse square root via diagonalization.
pub fn invertMetric(comptime T: type, J: Matrix(T), gpa: Allocator) !Matrix(T) {
    const naux = J.nrow();

    var W = try Vector(T).init(naux, gpa);
    defer W.deinit(gpa);

    var U = try Matrix(T).init(naux, naux, gpa);
    defer U.deinit(gpa);

    var work = try Matrix(T).init(naux, naux, gpa);
    defer work.deinit(gpa);

    var J_inv_half = try Matrix(T).initZero(naux, naux, gpa);
    errdefer J_inv_half.deinit(gpa);

    try syevd(T, &W, &U, J);

    for (0..naux) |p| {
        for (0..naux) |k| {
            const w_val = W.at(k);

            const scale = if (w_val > 1e-14) 1.0 / @sqrt(w_val) else 0.0;

            work.ptr(p, k).* = U.at(p, k) * scale;
        }
    }

    gemm(T, &J_inv_half, work, U, 1, 0, false, true);

    return J_inv_half;
}
