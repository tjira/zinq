//! Evaluates Hamiltonian matrix elements between Slater determinants using the Slater-Condon rules.

const Matrix = @import("tensor.zig").Matrix;
const Tensor = @import("tensor.zig").Tensor;
const Value = @import("value.zig").Value;

const findDeterminantDifferences = @import("determinant.zig").findDeterminantDifferences;
const primType = @import("value.zig").primType;
const signOfExcitations = @import("determinant.zig").signOfExcitations;

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
