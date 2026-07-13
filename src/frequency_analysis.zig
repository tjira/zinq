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

const AN2M = @import("constant.zig").AN2M;
const AN2SM = @import("constant.zig").AN2SM;
const AU2CM = @import("constant.zig").AU2CM;

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

        for (0..atoms.len) |i| for (0..3) |c| {
            com[c] += masses[i] * coors[i * 3 + c];
        };

        var total_mass: T = 0;

        for (0..atoms.len) |i| {
            total_mass += masses[i];
        }

        for (0..3) |c| {
            com[c] /= total_mass;
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
                G[i * 3 + k] = rx * inertia_vecs[0 * 3 + k] + ry * inertia_vecs[1 * 3 + k] + rz * inertia_vecs[2 * 3 + k];
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

pub fn printHarmonicFrequencies(comptime T: type, io: std.Io, freqs: Vector(T), method_str: []const u8) !void {
    try printf(io, "\n{s} HARMONIC VIBRATIONAL FREQUENCIES (cm^-1)\n", .{method_str});

    for (0..freqs.length()) |i| {
        try printf(io, "MODE {d:3}: {d:12.4}\n", .{ i + 1, AU2CM * freqs.at(i) });
    }
}
