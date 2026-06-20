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
         0.0000002178410248, -0.0000014890494290,  0.0000005387672022,
        -0.0000007176761931,  0.0000006273176250, -0.0000005399082537,
         0.0000004998351610,  0.0000008617318090,  0.0000000011410497,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.9659012172971900, hf.energy[0], TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(expected_grad[0], hf.gradient[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], hf.gradient[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], hf.gradient[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], hf.gradient[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], hf.gradient[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], hf.gradient[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], hf.gradient[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], hf.gradient[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], hf.gradient[0].at(2, 2), TEST_TOLERANCE);
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
         0.0000002178410003, -0.0000014890494329,  0.0000005387672015,
        -0.0000007176761749,  0.0000006273176159, -0.0000005399082350,
         0.0000004998351675,  0.0000008617318237,  0.0000000011410507,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.9659012172972300, hf.energy[0], TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(expected_grad[0], hf.gradient[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], hf.gradient[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], hf.gradient[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], hf.gradient[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], hf.gradient[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], hf.gradient[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], hf.gradient[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], hf.gradient[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], hf.gradient[0].at(2, 2), TEST_TOLERANCE);
}

test "Restricted Kohn-Sham DFT with LDA Functional on Water (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "example/basis/sto-3g.g94",
        .dft = .{
            .exchange = "lda_x",
            .correlation = "lda_c_vwn",
        },
        .gradient = false,
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-74.7399385581956500, hf.energy[0], TEST_TOLERANCE);
}

test "Generalized Kohn-Sham DFT with LDA Functional on Water (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "example/basis/sto-3g.g94",
        .dft = .{
            .exchange = "lda_x",
            .correlation = "lda_c_vwn",
        },
        .generalized = true,
        .gradient = false,
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-74.7399385581956500, hf.energy[0], TEST_TOLERANCE);
}
