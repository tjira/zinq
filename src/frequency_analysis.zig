//! Harmonic vibrational frequency analysis from molecular Hessians.

const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;
const Vector = @import("tensor.zig").Vector;

const eigh = @import("linear_algebra.zig").eigh;
const eighSlice = @import("linear_algebra.zig").eighSlice;
const getMass = @import("constant.zig").getMass;
const mm = @import("linear_algebra.zig").mm;
const printf = @import("read_write.zig").printf;

const a0 = @import("constant.zig").a0;
const h = @import("constant.zig").h;
const kB = @import("constant.zig").kB;
const mu = @import("constant.zig").mu;
const R = @import("constant.zig").R;

const CAL2J = @import("constant.zig").CAL2J;
const AU2CM = @import("constant.zig").AU2CM;
const AU2K = @import("constant.zig").AU2K;
const CM2AU = @import("constant.zig").CM2AU;

/// Options for thermochemical analysis specifying the thermodynamic temperature and pressure.
pub const Options = struct {
    pressure: f64 = 101325.0,
    temperature: f64 = 298.15,
};

/// Thermodynamic state functions derived from molecular partition functions under the RRHO approximation.
pub fn Thermochemistry(comptime T: type) type {
    return struct {
        E_thrm: T,
        E_zpve: T,
        G_corr: T,
        H_corr: T,
        M_mass: T,
        S_elec: T,
        S_rota: T,
        S_tran: T,
        S_vibr: T,
    };
}

/// Calculates mass-weighted harmonic vibrational frequencies and projects out translations and rotations.
pub fn calculateHarmonicFrequencies(comptime T: type, hessian: Matrix(T), sys: MolecularSystem(T), gpa: Allocator) !Vector(T) {
    std.debug.assert(hessian.nrow() == sys.atoms.len * 3);

    const atoms = sys.atoms;
    const coors = sys.coors;

    var HM = try Matrix(T).init(atoms.len * 3, atoms.len * 3, gpa);
    defer HM.deinit(gpa);

    for (0..atoms.len) |i| for (0..3) |xyz_i| for (0..atoms.len) |j| for (0..3) |xyz_j| {
        const row = i * 3 + xyz_i;
        const col = j * 3 + xyz_j;

        const mass_i = try getMass(T, atoms[i]);
        const mass_j = try getMass(T, atoms[j]);

        HM.ptr(row, col).* = hessian.at(row, col) / std.math.sqrt(mass_i * mass_j);
    };

    if (atoms.len > 0) {
        var masses = try gpa.alloc(T, atoms.len);
        defer gpa.free(masses);

        for (0..atoms.len) |i| {
            masses[i] = try getMass(T, atoms[i]);
        }

        var com = [3]T{ 0, 0, 0 };

        for (0..atoms.len) |i| for (0..3) |c_idx| {
            com[c_idx] += masses[i] * coors[i * 3 + c_idx];
        };

        var total_mass: T = 0;

        for (0..atoms.len) |i| {
            total_mass += masses[i];
        }

        for (0..3) |c_idx| {
            com[c_idx] /= total_mass;
        }

        var I_tensor = [9]T{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        for (0..atoms.len) |i| {
            const x = coors[i * 3 + 0] - com[0];
            const y = coors[i * 3 + 1] - com[1];
            const z = coors[i * 3 + 2] - com[2];

            I_tensor[0 * 3 + 0] += masses[i] * (y * y + z * z);
            I_tensor[1 * 3 + 1] += masses[i] * (x * x + z * z);
            I_tensor[2 * 3 + 2] += masses[i] * (x * x + y * y);

            I_tensor[0 * 3 + 1] -= masses[i] * x * y;
            I_tensor[1 * 3 + 0] -= masses[i] * x * y;

            I_tensor[0 * 3 + 2] -= masses[i] * x * z;
            I_tensor[2 * 3 + 0] -= masses[i] * x * z;

            I_tensor[1 * 3 + 2] -= masses[i] * y * z;
            I_tensor[2 * 3 + 1] -= masses[i] * y * z;
        }

        var inertia_ev: [3]T, var inertia_vecs: [9]T = .{ undefined, undefined };

        try eighSlice(T, &inertia_ev, &inertia_vecs, &I_tensor);

        var G = try gpa.alloc(T, atoms.len * 3);
        defer gpa.free(G);

        for (0..atoms.len) |i| {
            const rx = coors[i * 3 + 0] - com[0];
            const ry = coors[i * 3 + 1] - com[1];
            const rz = coors[i * 3 + 2] - com[2];

            for (0..3) |k| {
                const iv_0 = inertia_vecs[0 * 3 + k];
                const iv_1 = inertia_vecs[1 * 3 + k];
                const iv_2 = inertia_vecs[2 * 3 + k];

                G[i * 3 + k] = rx * iv_0 + ry * iv_1 + rz * iv_2;
            }
        }

        var TR = try gpa.alloc(T, atoms.len * 18);
        defer gpa.free(TR);

        @memset(TR, 0);

        for (0..atoms.len) |i| {
            const mass_12 = @sqrt(masses[i]);

            TR[(i * 3 + 0) * 6 + 0] = mass_12;
            TR[(i * 3 + 1) * 6 + 1] = mass_12;
            TR[(i * 3 + 2) * 6 + 2] = mass_12;

            for (0..3) |j| {
                const g1 = G[i * 3 + 1];
                const g2 = G[i * 3 + 2];
                const g0 = G[i * 3 + 0];

                const v2 = inertia_vecs[j * 3 + 2];
                const v1 = inertia_vecs[j * 3 + 1];
                const v0 = inertia_vecs[j * 3 + 0];

                TR[(i * 3 + j) * 6 + 3] = mass_12 * (g1 * v2 - g2 * v1);
                TR[(i * 3 + j) * 6 + 4] = mass_12 * (g2 * v0 - g0 * v2);
                TR[(i * 3 + j) * 6 + 5] = mass_12 * (g0 * v1 - g1 * v0);
            }
        }

        var U_tr = try gpa.alloc(T, atoms.len * 3 * 6);
        defer gpa.free(U_tr);

        @memset(U_tr, 0);

        var temp = try gpa.alloc(T, atoms.len * 3);
        defer gpa.free(temp);

        var num_tr: usize = 0;

        for (0..6) |col| {
            for (0..atoms.len * 3) |row| {
                temp[row] = TR[row * 6 + col];
            }

            for (0..num_tr) |j| {
                var dot: T = 0;

                for (0..atoms.len * 3) |row| {
                    dot += temp[row] * U_tr[row * 6 + j];
                }

                for (0..atoms.len * 3) |row| {
                    temp[row] -= dot * U_tr[row * 6 + j];
                }
            }

            var norm: T = 0;

            for (0..atoms.len * 3) |row| {
                norm += temp[row] * temp[row];
            }

            norm = @sqrt(norm);

            if (norm > 1e-6) {
                for (0..atoms.len * 3) |row| {
                    U_tr[row * 6 + num_tr] = temp[row] / norm;
                }

                num_tr += 1;
            }
        }

        var P_mat = try Matrix(T).init(atoms.len * 3, atoms.len * 3, gpa);
        defer P_mat.deinit(gpa);

        for (0..atoms.len * 3) |i| for (0..atoms.len * 3) |j| {
            var sum: T = 0;

            for (0..num_tr) |k| {
                sum += U_tr[i * 6 + k] * U_tr[j * 6 + k];
            }

            P_mat.ptr(i, j).* = (if (i == j) @as(T, 1) else @as(T, 0)) - sum;
        };

        var H_temp = try Matrix(T).init(atoms.len * 3, atoms.len * 3, gpa);
        defer H_temp.deinit(gpa);

        mm(T, &H_temp, P_mat, HM, 1, 0, false, false);
        mm(T, &HM, H_temp, P_mat, 1, 0, false, false);
    }

    var w = try Vector(T).init(atoms.len * 3, gpa);
    errdefer w.deinit(gpa);

    var u = try Matrix(T).init(hessian.nrow(), hessian.nrow(), gpa);
    defer u.deinit(gpa);

    try eigh(T, &w, &u, HM);

    for (0..hessian.nrow()) |i| {
        w.ptr(i).* = std.math.sign(w.at(i)) * std.math.sqrt(@abs(w.at(i)));
    }

    return w;
}

/// Evaluates thermodynamic state functions and partition functions under the ideal gas and RRHO approximations.
pub fn calculateThermochemistry(comptime T: type, opt: Options, sys: MolecularSystem(T), freqs: Vector(T), multiplicity: u32, gpa: Allocator) !Thermochemistry(T) {
    const temp, const pres, var total_mass: T = .{ opt.temperature, opt.pressure, 0 };

    var masses = try gpa.alloc(T, sys.atoms.len);
    defer gpa.free(masses);

    for (0..sys.atoms.len) |i| {
        masses[i] = try getMass(T, sys.atoms[i]);

        total_mass += masses[i];
    }

    var com = [3]T{ 0, 0, 0 };

    for (0..sys.atoms.len) |i| for (0..3) |c_idx| {
        com[c_idx] += masses[i] * sys.coors[i * 3 + c_idx];
    };

    for (0..3) |c_idx| {
        com[c_idx] /= total_mass;
    }

    var I_tensor = [9]T{ 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    for (0..sys.atoms.len) |i| {
        const x = sys.coors[i * 3 + 0] - com[0];
        const y = sys.coors[i * 3 + 1] - com[1];
        const z = sys.coors[i * 3 + 2] - com[2];

        I_tensor[0 * 3 + 0] += masses[i] * (y * y + z * z);
        I_tensor[1 * 3 + 1] += masses[i] * (x * x + z * z);
        I_tensor[2 * 3 + 2] += masses[i] * (x * x + y * y);

        I_tensor[0 * 3 + 1] -= masses[i] * x * y;
        I_tensor[1 * 3 + 0] -= masses[i] * x * y;

        I_tensor[0 * 3 + 2] -= masses[i] * x * z;
        I_tensor[2 * 3 + 0] -= masses[i] * x * z;

        I_tensor[1 * 3 + 2] -= masses[i] * y * z;
        I_tensor[2 * 3 + 1] -= masses[i] * y * z;
    }

    var inertia_ev: [3]T, var inertia_vecs: [9]T = .{ undefined, undefined };

    try eighSlice(T, &inertia_ev, &inertia_vecs, &I_tensor);

    const conv_I = mu * a0 * a0;

    const I0 = inertia_ev[0] * conv_I;
    const I1 = inertia_ev[1] * conv_I;
    const I2 = inertia_ev[2] * conv_I;

    var E_rot: T, var S_rota_over_kB: T = .{ 0, 0 };

    if (sys.atoms.len > 1) {
        const is_linear = inertia_ev[0] < 1e-4 * inertia_ev[2] or inertia_ev[0] < 1e-6;

        if (is_linear) {
            E_rot, const I_lin = .{ temp / AU2K, (I1 + I2) / 2 };

            const q_rot = (8 * std.math.pi * std.math.pi * I_lin * kB * temp) / (h * h);

            if (q_rot > 0) S_rota_over_kB = @log(q_rot) + 1;
        }

        if (!is_linear) {
            E_rot, const prefactor = .{ 1.5 * temp / AU2K, 8 * std.math.pi * std.math.pi * kB * temp / (h * h) };

            const q_rot = std.math.sqrt(std.math.pi) * std.math.sqrt(prefactor * prefactor * prefactor * I0 * I1 * I2);

            if (q_rot > 0) S_rota_over_kB = @log(q_rot) + 1.5;
        }
    }

    const lambda_term = (2 * std.math.pi * total_mass * mu * kB * temp) / (h * h);

    var E_zpve: T, var E_vib: T, var S_vibr_over_kB: T = .{ 0, 0, 0 };

    const S_tran_over_kB = @log(std.math.pow(T, lambda_term, 1.5) * (kB * temp / pres)) + 2.5;

    for (0..freqs.length()) |i| {
        const wn = AU2CM * freqs.at(i);

        if (wn > 1) {
            const eps = wn * CM2AU;

            E_zpve += 0.5 * eps;

            const x = (eps * AU2K) / temp;

            if (x < 100) {
                E_vib += eps / (@exp(x) - 1);

                S_vibr_over_kB += x / (@exp(x) - 1) - @log(1 - @exp(-x));
            }
        }
    }

    const S_elec_over_kB = @log(@as(T, @floatFromInt(@max(1, multiplicity))));

    const S_totl_over_kB = S_tran_over_kB + S_rota_over_kB + S_vibr_over_kB + S_elec_over_kB;

    const E_thrm = E_zpve + 1.5 * temp / AU2K + E_rot + E_vib;

    const H_corr, const G_corr = .{ E_thrm + temp / AU2K, E_thrm + temp / AU2K - (temp / AU2K) * S_totl_over_kB };

    return .{
        .M_mass = total_mass,

        .E_thrm = E_thrm,
        .E_zpve = E_zpve,
        .G_corr = G_corr,
        .H_corr = H_corr,

        .S_elec = S_elec_over_kB * R / CAL2J,
        .S_rota = S_rota_over_kB * R / CAL2J,
        .S_tran = S_tran_over_kB * R / CAL2J,
        .S_vibr = S_vibr_over_kB * R / CAL2J,
    };
}

/// Formats and prints harmonic frequencies in wavenumbers (cm^-1) to the output stream.
pub fn printHarmonicFrequencies(comptime T: type, io: std.Io, freqs: Vector(T), method_str: []const u8) !void {
    try printf(io, "\n{s} HARMONIC VIBRATIONAL FREQUENCIES (cm^-1)\n", .{method_str});

    for (0..freqs.length()) |i| {
        try printf(io, "MODE {d:02}: {d:13.4}\n", .{ i + 1, AU2CM * freqs.at(i) });
    }
}

/// Formats and prints thermochemical properties and thermodynamic functions to the output stream.
pub fn printThermochemistry(comptime T: type, io: std.Io, opt: Options, th: Thermochemistry(T), energy: T, method_str: []const u8) !void {
    try printf(io, "\n{s} SYSTEM CONDITIONS\n", .{method_str});

    try printf(io, "T_TEMP: {d:20.14} K\nP_PRES: {d:20.14} kPa\nM_MASS: {d:20.14} amu\n", .{ opt.temperature, opt.pressure / 1000, th.M_mass });

    try printf(io, "\n{s} THERMAL CORRECTIONS\n", .{method_str});

    try printf(io, "E_ZPVE: {d:20.14} Eh\n", .{th.E_zpve});
    try printf(io, "E_THRM: {d:20.14} Eh\n", .{th.E_thrm});
    try printf(io, "H_CORR: {d:20.14} Eh\n", .{th.H_corr});
    try printf(io, "G_CORR: {d:20.14} Eh\n", .{th.G_corr});

    try printf(io, "\n{s} ENTROPIES (cal/(mol K))\n", .{method_str});

    const S_totl = th.S_tran + th.S_rota + th.S_vibr + th.S_elec;

    try printf(io, "S_TRAN: {d:20.14}\n", .{th.S_tran});
    try printf(io, "S_ROTA: {d:20.14}\n", .{th.S_rota});
    try printf(io, "S_VIBR: {d:20.14}\n", .{th.S_vibr});
    try printf(io, "S_ELEC: {d:20.14}\n", .{th.S_elec});
    try printf(io, "S_TOTL: {d:20.14}\n", .{S_totl});

    try printf(io, "\n{s} SUM OF ELECTRONIC AND THERMAL ENERGIES\n", .{method_str});

    try printf(io, "E_ELEC: {d:20.14} Eh\n", .{energy});
    try printf(io, "E_ZERO: {d:20.14} Eh\n", .{energy + th.E_zpve});
    try printf(io, "E_TOTL: {d:20.14} Eh\n", .{energy + th.E_thrm});
    try printf(io, "H_TOTL: {d:20.14} Eh\n", .{energy + th.H_corr});
    try printf(io, "G_TOTL: {d:20.14} Eh\n", .{energy + th.G_corr});
}
