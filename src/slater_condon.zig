//! Evaluates Hamiltonian matrix elements between Slater determinants using the Slater-Condon rules.

const std = @import("std");

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;
const Value = @import("value.zig").Value;

const calculateDeterminantSpin = @import("determinant.zig").calculateDeterminantSpin;
const findDeterminantDifferences = @import("determinant.zig").findDeterminantDifferences;
const primType = @import("value.zig").primType;
const signOfExcitations = @import("determinant.zig").signOfExcitations;

/// Adds an S^2 penalty operator to the CI Hamiltonian to shift non-target spin states to higher energies.
pub fn applySpinPenalty(comptime T: type, H: *Matrix(T), dets: []const []const usize, lambda: f64, s2_target: f64) void {
    const val_lambda = Value(T).fromFloat(lambda);

    for (0..dets.len) |i| {
        const s2_diag = calculateDeterminantSpin(dets[i]);

        const shift_diag = val_lambda.mul(Value(T).fromFloat(s2_diag - s2_target));

        H.ptr(i, i).* = Value(T).init(H.at(i, i)).add(shift_diag).val;

        for (i + 1..dets.len) |j| {
            const diff = findDeterminantDifferences(dets[i], dets[j]) orelse continue;

            if (diff.len != 2) continue;

            const p0, const p1 = .{ diff.a[0] / 2, diff.a[1] / 2 };
            const q0, const q1 = .{ diff.b[0] / 2, diff.b[1] / 2 };

            if (p0 == q0 and p1 == q1 and (diff.a[0] % 2 != diff.a[1] % 2) and (diff.b[0] % 2 != diff.a[0] % 2)) {
                const s_A = signOfExcitations(dets[i], diff.a[0..2]);
                const s_B = signOfExcitations(dets[j], diff.b[0..2]);

                const phase = @as(f64, @floatFromInt(s_A * s_B));

                const shift_off = val_lambda.mul(Value(T).fromFloat(phase));

                H.ptr(i, j).* = Value(T).init(H.at(i, j)).add(shift_off).val;
                H.ptr(j, i).* = Value(T).init(H.at(j, i)).add(shift_off).val;
            }
        }
    }
}

/// Assembles the electronic Hamiltonian matrix across Slater determinants via the Slater-Condon rules.
pub fn buildHamiltonianMatrix(comptime T: type, dets: []const []const usize, h_ms: Matrix(T), g_ms: Tensor(T, 4), gpa: Allocator) !Matrix(T) {
    var H = try Matrix(T).init(dets.len, dets.len, gpa);
    errdefer H.deinit(gpa);

    for (0..dets.len) |i| for (i..dets.len) |j| {
        const val = slater(T, dets[i], dets[j], h_ms, g_ms);

        H.ptr(i, j).* = val;
        H.ptr(j, i).* = val;
    };

    return H;
}

/// Computes the Hamiltonian matrix element between two Slater determinants using the Slater-Condon rules.
pub fn slater(comptime T: type, A: []const usize, B: []const usize, h_ms: Matrix(T), g_ms: Tensor(T, 4)) T {
    const diff = findDeterminantDifferences(A, B) orelse return Value(T).fromFloat(0).val;

    const s_A = signOfExcitations(A, diff.a[0..diff.len]);
    const s_B = signOfExcitations(B, diff.b[0..diff.len]);

    const phase = @as(primType(T), @floatFromInt(s_A * s_B));

    return evaluateSlaterCondon(T, A, diff.a, diff.b, diff.len, phase, h_ms, g_ms);
}

/// Evaluates Hamiltonian matrix elements between Slater determinants using the Slater-Condon excitation rules.
fn evaluateSlaterCondon(comptime T: type, A: []const usize, diff_a: [2]usize, diff_b: [2]usize, len: usize, phase: primType(T), h_ms: Matrix(T), g_ms: Tensor(T, 4)) T {
    return switch (len) {
        0 => slaterDiagonal(T, h_ms, g_ms, A),
        1 => slaterSingle(T, h_ms, g_ms, A, diff_a[0], diff_b[0], phase),
        2 => slaterDouble(T, g_ms, diff_a[0], diff_a[1], diff_b[0], diff_b[1], phase),
        else => Value(T).fromFloat(0).val,
    };
}

/// Evaluates diagonal Hamiltonian matrix elements between identical Slater determinants via Slater-Condon rules.
fn slaterDiagonal(comptime T: type, h_ms: Matrix(T), g_ms: Tensor(T, 4), det: []const usize) T {
    var sum = Value(T).fromFloat(0);

    for (0..det.len) |i| {
        sum = sum.add(Value(T).init(h_ms.at(det[i], det[i])));
    }

    for (0..det.len) |i| for (i + 1..det.len) |j| {
        const p, const q = .{ det[i], det[j] };

        const coulomb = Value(T).init(g_ms.at(.{ p, q, p, q }));

        const exchange = Value(T).init(g_ms.at(.{ p, q, q, p }));

        sum = sum.add(coulomb.sub(exchange));
    };

    return sum.val;
}

/// Evaluates Hamiltonian matrix elements between Slater determinants differing by two spin-orbitals.
fn slaterDouble(comptime T: type, g_ms: Tensor(T, 4), p: usize, q: usize, r: usize, s: usize, phase: primType(T)) T {
    const direct = Value(T).init(g_ms.at(.{ p, q, r, s }));

    const exchange = Value(T).init(g_ms.at(.{ p, q, s, r }));

    return direct.sub(exchange).muls(phase).val;
}

/// Evaluates Hamiltonian matrix elements between Slater determinants differing by a single spin-orbital.
fn slaterSingle(comptime T: type, h_ms: Matrix(T), g_ms: Tensor(T, 4), det: []const usize, p: usize, q: usize, phase: primType(T)) T {
    var elem = Value(T).init(h_ms.at(p, q));

    for (0..det.len) |idx| if (det[idx] != p) {
        const k = det[idx];

        const coulomb = Value(T).init(g_ms.at(.{ p, k, q, k }));

        const exchange = Value(T).init(g_ms.at(.{ p, k, k, q }));

        elem = elem.add(coulomb.sub(exchange));
    };

    return elem.muls(phase).val;
}
