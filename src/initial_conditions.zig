//! Samples initial phase space positions and momenta from classical thermal and quantum Wigner distributions.

const std = @import("std");

const Allocator = std.mem.Allocator;

const electronic_structure = @import("electronic_structure.zig");

const Ensemble = @import("ensemble.zig").Ensemble;
const Matrix = @import("tensor.zig").Matrix;
const MethodOptions = @import("potential.zig").MethodOptions;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Vector = @import("tensor.zig").Vector;

const calculateNormalModes = @import("frequency_analysis.zig").calculateNormalModes;
const stripSecondaryOptions = @import("nuclear_derivative.zig").stripSecondaryOptions;

const AMU2AU = @import("constants.zig").AMU2AU;
const AU2K = @import("constants.zig").AU2K;

/// Initial phase space parameters and distribution configuration for trajectory sampling.
pub const InitialConditions = struct {
    distribution: ?DistributionOptions = null,
    momentum: ?[]const f64 = null,
    position: ?[]const f64 = null,
    seed: u32 = 1,
    state: u32 = 0,
};

/// Tagged union selecting thermal Maxwell-Boltzmann or quantum harmonic Wigner distributions.
pub const DistributionOptions = union(enum) {
    maxwell_boltzmann: MaxwellBoltzmannOptions,
    wigner: WignerOptions,
};

/// Configuration options for thermal Maxwell-Boltzmann momentum sampling.
pub const MaxwellBoltzmannOptions = struct {
    temperature: f64 = 298.15,
};

/// Parameters for quantum harmonic oscillator Wigner phase space distribution sampling.
pub const WignerOptions = struct {
    gamma: ?[]const f64 = null,
    method: ?MethodOptions = null,
    temperature: ?f64 = null,
};

/// Extracts or configures electronic structure method options for ab initio Wigner Hessian calculations.
pub fn getWignerMethod(opt: anytype) MethodOptions {
    if (opt.initial_conditions.distribution.?.wigner.method) |m| {
        return stripSecondaryOptions(m, .{ .hessian = true });
    }

    var method = stripSecondaryOptions(opt.potential.ab_initio, .{});

    switch (method) {
        inline else => |*m| m.hessian = .{ .numeric = .{} },
    }

    return method;
}

/// Executes ab initio electronic structure calculation to obtain Hessian and samples Wigner distribution.
pub fn sampleWignerAbInitio(comptime T: type, io: std.Io, ensemble: *Ensemble(T), opt: anytype, gpa: Allocator) !void {
    switch (getWignerMethod(opt)) {
        inline else => |m_opt| {
            var res = try electronic_structure.run(T, io, m_opt, true, gpa);
            defer res.deinit(gpa);

            const sys = electronic_structure.getSystem(T, res);

            try sampleWignerFromHessian(T, ensemble, res.hess[0], sys, opt, gpa);
        },
    }
}

/// Samples ensemble positions and momenta from normal modes computed from a molecular Hessian.
pub fn sampleWignerFromHessian(comptime T: type, ensemble: *Ensemble(T), hess: Matrix(T), sys: MolecularSystem(T), opt: anytype, gpa: Allocator) !void {
    var nm = try calculateNormalModes(T, hess, sys, gpa);
    defer nm.deinit(gpa);

    const ndim, const ntraj = .{ sys.atoms.len * 3, ensemble.r.nrow() };

    const temp = opt.initial_conditions.distribution.?.wigner.temperature orelse 0;

    var split_mix = std.Random.SplitMix64.init(opt.initial_conditions.seed);
    var rng = std.Random.DefaultPrng.init(split_mix.next());

    const sigma_Q = try gpa.alloc(T, ndim);
    defer gpa.free(sigma_Q);

    const sigma_P = try gpa.alloc(T, ndim);
    defer gpa.free(sigma_P);

    computeModeWidths(T, sigma_Q, sigma_P, nm.freqs, temp);

    const Q = try gpa.alloc(T, ndim);
    defer gpa.free(Q);

    const P_Q = try gpa.alloc(T, ndim);
    defer gpa.free(P_Q);

    const random = rng.random();

    for (0..ntraj) |i| {
        ensemble.s.ptr(i).* = opt.initial_conditions.state;

        sampleNormalDisplacements(T, Q, P_Q, sigma_Q, sigma_P, random);

        transformModesToCartesian(T, ensemble, sys, nm.modes, Q, P_Q, i);
    }
}

/// Computes quantum harmonic oscillator normal coordinate and momentum standard deviations.
fn computeModeWidths(comptime T: type, sigma_Q: []T, sigma_P: []T, freqs: Vector(T), temp: f64) void {
    for (0..freqs.length()) |k| {
        const omega_au = freqs.at(k) / @sqrt(AMU2AU);

        if (omega_au <= 1e-5) {
            sigma_Q[k] = 0;
            sigma_P[k] = 0;

            continue;
        }

        if (temp > 1e-6) {
            const coth = 1.0 / std.math.tanh(0.5 * (AU2K / @as(T, @floatCast(temp))) * omega_au);

            sigma_Q[k] = std.math.sqrt(0.5 / omega_au * coth);
            sigma_P[k] = std.math.sqrt(0.5 * omega_au * coth);
        }

        if (temp <= 1e-6) {
            sigma_Q[k] = std.math.sqrt(0.5 / omega_au);
            sigma_P[k] = std.math.sqrt(0.5 * omega_au);
        }
    }
}

/// Draws uncorrelated Gaussian random displacements in normal mode coordinate space.
fn sampleNormalDisplacements(comptime T: type, Q: []T, P_Q: []T, sigma_Q: []const T, sigma_P: []const T, random: anytype) void {
    for (0..Q.len) |k| {
        Q[k] = if (sigma_Q[k] > 0) sigma_Q[k] * random.floatNorm(T) else 0;

        P_Q[k] = if (sigma_P[k] > 0) sigma_P[k] * random.floatNorm(T) else 0;
    }
}

/// Transforms normal mode displacements and momenta into mass-weighted Cartesian coordinates.
fn transformModesToCartesian(comptime T: type, ensemble: *Ensemble(T), sys: MolecularSystem(T), modes: Matrix(T), Q: []const T, P_Q: []const T, traj: usize) void {
    const ndim = sys.atoms.len * 3;

    for (0..ndim) |a| {
        var delta_r: T, var p_a: T = .{ 0, 0 };

        for (0..ndim) |k| {
            const u_ak = modes.at(a, k);

            delta_r += u_ak * Q[k];

            p_a += u_ak * P_Q[k];
        }

        const sqrt_m = std.math.sqrt(ensemble.m[a]);

        ensemble.r.ptr(traj, a).* = sys.coors[a] + delta_r / sqrt_m;

        ensemble.p.ptr(traj, a).* = p_a * sqrt_m;
    }
}
