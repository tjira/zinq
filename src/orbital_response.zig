//! Solves the Coupled-Perturbed Hartree-Fock equations to obtain the first-order orbital response to perturbations.

const std = @import("std");

const Allocator = std.mem.Allocator;

const Diis = @import("diis.zig").Diis;
const HartreeFockResult = @import("hartree_fock.zig").Result;
const Matrix = @import("tensor.zig").Matrix;
const FittingDerivativeWorkspace = @import("density_fitting.zig").FittingDerivativeWorkspace;
const RiFlags = @import("density_fitting.zig").RiFlags;
const RiWorkspace = @import("density_fitting.zig").RiWorkspace;
const Tensor = @import("tensor.zig").Tensor;
const Vector = @import("tensor.zig").Vector;

const ao2mo_1e = @import("integral_transform.zig").ao2mo_1e;
const contractFockDerivative = @import("two_electron_contraction.zig").contractFockDerivative;
const contractResponseHessian = @import("two_electron_contraction.zig").contractResponseHessian;
const gemm = @import("linear_algebra.zig").gemm;
const mo2ao_1e = @import("integral_transform.zig").mo2ao_1e;
const printf = @import("read_write.zig").printf;

/// Scratch memory buffers for evaluating non-relaxed perturbed molecular orbital coefficients.
pub fn NonRelaxedOrbitalWorkspace(comptime T: type) type {
    return struct {
        temp: Matrix(T),
        S_mo: Matrix(T),

        /// Allocates intermediate scratch matrices for non-relaxed orbital coefficient evaluation.
        pub fn init(nbf: usize, nmo: usize, gpa: Allocator) !@This() {
            var temp = try Matrix(T).init(nbf, nmo, gpa);
            errdefer temp.deinit(gpa);

            var S_mo = try Matrix(T).init(nmo, nmo, gpa);
            errdefer S_mo.deinit(gpa);

            return .{ .temp = temp, .S_mo = S_mo };
        }

        /// Deallocates intermediate scratch matrices.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.S_mo.deinit(gpa);
            self.temp.deinit(gpa);
        }

        /// Evaluates non-relaxed perturbed coefficients dC = -0.5 * C * (C^T * dS * C).
        pub fn evaluate(self: *@This(), dC: *Matrix(T), C: Matrix(T), dS_ao: Matrix(T)) void {
            gemm(T, &self.temp, dS_ao, C, 1, 0, false, false);

            gemm(T, &self.S_mo, C, self.temp, 1, 0, true, false);

            gemm(T, dC, C, self.S_mo, -0.5, 0, false, false);
        }
    };
}

/// Scratch memory buffers for storing perturbed density, potential, and orbital rotation matrices during CPHF iterations.
fn CphfWorkspace(comptime T: type) type {
    return struct {
        F_x: Matrix(T),
        P_x: Matrix(T),
        U_x: Matrix(T),
        U_p: Matrix(T),
        V_x: Matrix(T),

        F_x_mo: Matrix(T),
        S_x_mo: Matrix(T),
        V_x_mo: Matrix(T),

        dC_x: Matrix(T),

        /// Allocates intermediate scratch matrices required for CPHF orbital response calculations.
        pub fn init(nbf: usize, gpa: Allocator) !@This() {
            var F_x = try Matrix(T).init(nbf, nbf, gpa);
            errdefer F_x.deinit(gpa);

            var P_x = try Matrix(T).init(nbf, nbf, gpa);
            errdefer P_x.deinit(gpa);

            var U_x = try Matrix(T).init(nbf, nbf, gpa);
            errdefer U_x.deinit(gpa);

            var U_p = try Matrix(T).init(nbf, nbf, gpa);
            errdefer U_p.deinit(gpa);

            var V_x = try Matrix(T).init(nbf, nbf, gpa);
            errdefer V_x.deinit(gpa);

            var F_x_mo = try Matrix(T).init(nbf, nbf, gpa);
            errdefer F_x_mo.deinit(gpa);

            var S_x_mo = try Matrix(T).init(nbf, nbf, gpa);
            errdefer S_x_mo.deinit(gpa);

            var V_x_mo = try Matrix(T).init(nbf, nbf, gpa);
            errdefer V_x_mo.deinit(gpa);

            var dC_x = try Matrix(T).init(nbf, nbf, gpa);
            errdefer dC_x.deinit(gpa);

            return .{
                .F_x = F_x,
                .P_x = P_x,
                .U_x = U_x,
                .U_p = U_p,
                .V_x = V_x,

                .F_x_mo = F_x_mo,
                .S_x_mo = S_x_mo,
                .V_x_mo = V_x_mo,

                .dC_x = dC_x,
            };
        }

        /// Deallocates all intermediate scratch matrices allocated for CPHF iterations.
        pub fn deinit(self: *@This(), gpa: Allocator) void {
            self.F_x.deinit(gpa);
            self.P_x.deinit(gpa);
            self.U_x.deinit(gpa);
            self.U_p.deinit(gpa);
            self.V_x.deinit(gpa);

            self.F_x_mo.deinit(gpa);
            self.S_x_mo.deinit(gpa);
            self.V_x_mo.deinit(gpa);

            self.dC_x.deinit(gpa);
        }
    };
}

/// Computes the first-order derivative of molecular orbital coefficients and energies using CPHF.
pub fn orbitalResponse(comptime T: type, io: std.Io, hfres: HartreeFockResult(T), opt: anytype, log: bool, gpa: Allocator) !struct { Tensor(T, 3), Matrix(T) } {
    const nbf = hfres.C_mo.shape[0];

    var dC = try Tensor(T, 3).initZero(.{ 3 * hfres.ints.sys.atoms.len, nbf, nbf }, gpa);
    errdefer dC.deinit(gpa);

    var de = try Matrix(T).initZero(3 * hfres.ints.sys.atoms.len, nbf, gpa);
    errdefer de.deinit(gpa);

    try solveCPHF(T, io, &dC, &de, hfres, opt, log, gpa);

    return .{ dC, de };
}

/// Solves the Coupled-Perturbed Hartree-Fock equations to obtain the first-order orbital response to perturbations.
pub fn solveCPHF(comptime T: type, io: std.Io, dC: *Tensor(T, 3), de: *Matrix(T), hfres: HartreeFockResult(T), opt: anytype, log: bool, gpa: Allocator) !void {
    const generalized = hfres.ints.sys.nbf != hfres.C_mo.shape[0];

    const dS = hfres.ints.dS orelse return error.MissingIntegrals;
    const dH = hfres.ints.dH orelse return error.MissingIntegrals;

    if (hfres.ints.B == null and (hfres.ints.g == null or hfres.ints.dg == null)) {
        return error.MissingIntegrals;
    }

    const C, const P, const e = .{ hfres.C_mo, hfres.P_ao, hfres.e_mo };

    const nocc, const nbf = .{ if (generalized) hfres.ints.sys.nel else hfres.ints.sys.nel / 2, C.shape[0] };

    const exch_factor: T = if (generalized) 1 else 0.5;
    const dens_factor: T = if (generalized) 1 else 2.0;

    var ws = try CphfWorkspace(T).init(nbf, gpa);
    defer ws.deinit(gpa);

    var ri_ws: ?RiWorkspace(T) = null;
    defer if (ri_ws) |*w| w.deinit(gpa);

    var deriv_ws: ?FittingDerivativeWorkspace(T) = null;
    defer if (deriv_ws) |*w| w.deinit(gpa);

    if (hfres.ints.B) |B| {
        const J = hfres.ints.J orelse return error.MissingIntegrals;

        _ = hfres.ints.dI orelse return error.MissingIntegrals;
        _ = hfres.ints.dJ orelse return error.MissingIntegrals;

        _ = hfres.ints.I orelse return error.MissingIntegrals;

        deriv_ws = try FittingDerivativeWorkspace(T).init(J, B.shape, gpa);

        const flags: RiFlags = .{
            .coulomb = true,
            .exchange = true,
            .generalized = generalized,
            .response = true,
        };

        ri_ws = try RiWorkspace(T).init(hfres.ints.sys.nbf, J.nrow(), nocc, flags, gpa);

        ri_ws.?.updateExchange(B, C);
    }

    for (0..3 * hfres.ints.sys.atoms.len) |c| {
        if (log) {
            try printf(io, "\nSOLVING CPHF FOR COORDINATE {d}\n", .{c + 1});

            try printf(io, "{s:4} {s:9} {s:9} {s:9}\n", .{ "ITER", "RMS(DU)", "RMS(DP)", "TIME" });
        }

        var diis_acc: ?Diis(T) = if (opt.diis) |d| if (d > 0) Diis(T).init(d) else null else null;
        defer if (diis_acc) |*d| d.deinit(gpa);

        ws.F_x.zero();

        for (0..nbf) |mu| for (0..nbf) |nu| {
            ws.F_x.ptr(mu, nu).* = dH.at(.{ c, mu, nu });
        };

        if (hfres.ints.B) |B| {
            _ = try deriv_ws.?.evaluateCoord(hfres.ints, c, gpa);

            ri_ws.?.contractFockDerivative(&ws.F_x, B, deriv_ws.?.dB_ao, P, C);
        }

        if (hfres.ints.dg) |dg| {
            contractFockDerivative(T, &ws.F_x, dg, c, P, exch_factor);
        }

        ws.F_x.symmetrize();

        try ao2mo_1e(T, &ws.F_x_mo, ws.F_x, C, gpa);

        try ao2mo_1e(T, &ws.S_x_mo, dS.slice(c).asMatrix(), C, gpa);

        initOrbitalRotations(T, &ws.U_x, ws.F_x_mo, ws.S_x_mo, e, nocc);

        ws.P_x.zero();

        for (0..opt.iterations) |i| {
            var timer = std.Io.Timestamp.now(io, .real);

            if (opt.diis) |d| if (d > 0) for (0..nbf) |p| for (0..nbf) |q| {
                ws.U_p.ptr(p, q).* = ws.U_x.at(p, q);
            };

            try mo2ao_1e(T, &ws.dC_x, ws.U_x, C);

            const rms_dp = computeResponseDensity(T, &ws.P_x, ws.dC_x, C, nocc, dens_factor);

            ws.V_x.zero();

            if (hfres.ints.B) |B| {
                ri_ws.?.contractResponseHessian(&ws.V_x, B, ws.P_x, ws.dC_x);
            }

            if (hfres.ints.g) |g| {
                contractResponseHessian(T, &ws.V_x, g, ws.P_x, exch_factor);
            }

            ws.V_x.symmetrize();

            try ao2mo_1e(T, &ws.V_x_mo, ws.V_x, C, gpa);

            const rms_du = updateOrbitalRotations(T, &ws.U_x, ws.F_x_mo, ws.V_x_mo, ws.S_x_mo, e, nocc);

            if (diis_acc) |*d| {
                var e_diis = try Matrix(T).init(nbf, nbf, gpa);
                errdefer e_diis.deinit(gpa);

                for (0..nbf) |p| for (0..nbf) |q| {
                    e_diis.ptr(p, q).* = ws.U_x.at(p, q) - ws.U_p.at(p, q);
                };

                try d.update(&ws.U_x, e_diis, false, gpa);

                for (0..nocc) |k| for (nocc..nbf) |a| {
                    ws.U_x.ptr(k, a).* = -ws.S_x_mo.at(k, a) - ws.U_x.at(a, k);
                };
            }

            const elapsed = timer.untilNow(io, .real);

            if (log) {
                const fmt = "{d:4} {e:9.3} {e:9.3} {s:9} {s}\n";

                const time_str = try std.fmt.allocPrint(gpa, "{f}", .{elapsed});
                defer gpa.free(time_str);

                const diis_tag = if (diis_acc) |d| (if (d.err_hist.items.len >= 2) "DIIS" else "") else "";

                try printf(io, fmt, .{ i + 1, rms_du, rms_dp, time_str, diis_tag });
            }

            if (rms_du < opt.threshold and rms_dp < opt.threshold) {
                break;
            }
        } else return error.CphfDidNotConverge;

        try finalizeCoordinateResponse(T, dC, de, &ws, c, C, e);
    }
}

/// Initializes the orbital rotation matrix from MO overlap and bare Fock derivatives.
fn initOrbitalRotations(comptime T: type, U_x: *Matrix(T), F_x_mo: Matrix(T), S_x_mo: Matrix(T), e: Vector(T), nocc: usize) void {
    const nbf = U_x.nrow();

    for (0..nbf) |p| for (0..nbf) |q| {
        U_x.ptr(p, q).* = -0.5 * S_x_mo.at(p, q);
    };

    for (nocc..nbf) |a| for (0..nocc) |i| {
        const diff = e.at(a) - e.at(i);

        if (@abs(diff) > 1e-12) {
            const u_ai = (e.at(i) * S_x_mo.at(a, i) - F_x_mo.at(a, i)) / diff;

            U_x.ptr(a, i).*, U_x.ptr(i, a).* = .{ u_ai, -S_x_mo.at(i, a) - u_ai };
        }
    };
}

/// Forms the first-order response density matrix from perturbed orbital coefficients.
fn computeResponseDensity(comptime T: type, P_x: *Matrix(T), dC_x: Matrix(T), C: Matrix(T), nocc: usize, dens_factor: T) T {
    const nbf = P_x.nrow();

    var sum_sq_dp: T = 0;

    for (0..nbf) |mu| for (0..nbf) |nu| {
        var sum: T = 0;

        for (0..nocc) |k| {
            sum += dens_factor * (dC_x.at(mu, k) * C.at(nu, k) + C.at(mu, k) * dC_x.at(nu, k));
        }

        const diff = sum - P_x.at(mu, nu);

        sum_sq_dp += diff * diff;

        P_x.ptr(mu, nu).* = sum;
    };

    return @sqrt(sum_sq_dp / @as(T, @floatFromInt(nbf * nbf)));
}

/// Updates occupied-virtual orbital rotations using the response potential and overlap derivative.
fn updateOrbitalRotations(comptime T: type, U_x: *Matrix(T), F_x_mo: Matrix(T), V_x_mo: Matrix(T), S_x_mo: Matrix(T), e: Vector(T), nocc: usize) T {
    const nbf = U_x.nrow();

    var sum_sq_du: T = 0;

    for (nocc..nbf) |a| for (0..nocc) |k| {
        const prev, const delta_e = .{ U_x.at(a, k), e.at(a) - e.at(k) };

        const f = F_x_mo.at(a, k);
        const v = V_x_mo.at(a, k);
        const s = S_x_mo.at(a, k);

        U_x.ptr(a, k).* = if (@abs(delta_e) > 1e-12) (e.at(k) * s - f - v) / delta_e else -0.5 * s;

        const diff = U_x.at(a, k) - prev;

        sum_sq_du += diff * diff;
    };

    const rms_du = @sqrt(sum_sq_du / @as(T, @floatFromInt((nbf - nocc) * nocc)));

    for (0..nocc) |k| for (nocc..nbf) |a| {
        U_x.ptr(k, a).* = -S_x_mo.at(k, a) - U_x.at(a, k);
    };

    return rms_du;
}

/// Computes complete orbital rotation matrices, back-transforms to AO basis, and evaluates energy derivatives.
fn finalizeCoordinateResponse(comptime T: type, dC: *Tensor(T, 3), de: *Matrix(T), ws: *CphfWorkspace(T), c: usize, C: Matrix(T), e: Vector(T)) !void {
    const nbf = ws.U_x.nrow();

    for (0..nbf) |p| for (0..nbf) |q| {
        const delta_e = e.at(q) - e.at(p);

        const f = ws.F_x_mo.at(p, q);
        const v = ws.V_x_mo.at(p, q);
        const s = ws.S_x_mo.at(p, q);

        ws.U_x.ptr(p, q).* = if (@abs(delta_e) > 1e-12) (f + v - e.at(q) * s) / delta_e else -0.5 * s;
    };

    try mo2ao_1e(T, &ws.dC_x, ws.U_x, C);

    for (0..nbf) |mu| for (0..nbf) |p| {
        dC.ptr(.{ c, mu, p }).* = ws.dC_x.at(mu, p);
    };

    for (0..nbf) |p| {
        de.ptr(c, p).* = ws.F_x_mo.at(p, p) + ws.V_x_mo.at(p, p) - ws.S_x_mo.at(p, p) * e.at(p);
    }
}
