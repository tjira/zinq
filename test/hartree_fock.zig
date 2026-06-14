const std = @import("std");
const zinq = @import("zinq");

const TEST_TOLERANCE = 1e-12;

test "Restricted Hartree-Fock on Water (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "example/basis/sto-3g.g94",
        .gradient = true,
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0000002231697467, -0.0000014871958312,  0.0000005408963109,
        -0.0000007203022483,  0.0000006261685166, -0.0000005408894763,
         0.0000004971324895,  0.0000008610273134, -0.0000000000068363,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.96590121729727, hf.energy, 1e-12);

    try std.testing.expectApproxEqAbs(expected_grad[0], hf.G.?.at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], hf.G.?.at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], hf.G.?.at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], hf.G.?.at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], hf.G.?.at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], hf.G.?.at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], hf.G.?.at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], hf.G.?.at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], hf.G.?.at(2, 2), TEST_TOLERANCE);
}

test "Generalized Hartree-Fock on Water (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "example/basis/sto-3g.g94",
        .generalized = true,
        .gradient = true,
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0000002299250568, -0.0000014848555010,  0.0000005435981549,
        -0.0000007236326030,  0.0000006247526345, -0.0000005421458191,
         0.0000004937075367,  0.0000008601028569, -0.0000000014523428,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.96590121729730, hf.energy, 1e-12);

    try std.testing.expectApproxEqAbs(expected_grad[0], hf.G.?.at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], hf.G.?.at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], hf.G.?.at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], hf.G.?.at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], hf.G.?.at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], hf.G.?.at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], hf.G.?.at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], hf.G.?.at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], hf.G.?.at(2, 2), TEST_TOLERANCE);
}
