const std = @import("std");
const zinq = @import("zinq");

const TEST_TOLERANCE = 1e-8;

test "Restricted Hartree-Fock on Water (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .gradient = .{ .analytic = .{} },
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
        .basis = "builtin:sto-3g",
        .generalized = true,
        .gradient = .{ .analytic = .{} },
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

test "Restricted Kohn-Sham DFT with LDA (VWN5) Functional on Water (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "lda_x",
            .correlation = "lda_c_vwn",
        },
        .gradient = null,
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-74.7399385581956500, hf.energy[0], TEST_TOLERANCE);
}

test "Generalized Kohn-Sham DFT with LDA (VWN5) Functional on Water (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "lda_x",
            .correlation = "lda_c_vwn",
        },
        .generalized = true,
        .gradient = null,
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-74.7399385581956500, hf.energy[0], TEST_TOLERANCE);
}

test "Restricted Kohn-Sham DFT with GGA (PBE) Functional on Water (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "gga_x_pbe",
            .correlation = "gga_c_pbe",
        },
        .gradient = null,
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.23434062179214, hf.energy[0], TEST_TOLERANCE);
}

test "Generalized Kohn-Sham DFT with GGA (PBE) Functional on Water (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "gga_x_pbe",
            .correlation = "gga_c_pbe",
        },
        .generalized = true,
        .gradient = null,
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.23434062179214, hf.energy[0], TEST_TOLERANCE);
}

test "Restricted Kohn-Sham DFT with meta-GGA (TPSS) Functional on Water (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "mgga_x_tpss",
            .correlation = "mgga_c_tpss",
        },
        .gradient = null,
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.33570008867521, hf.energy[0], TEST_TOLERANCE);
}

test "Generalized Kohn-Sham DFT with meta-GGA (TPSS) Functional on Water (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "mgga_x_tpss",
            .correlation = "mgga_c_tpss",
        },
        .generalized = true,
        .gradient = null,
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.33570008867527, hf.energy[0], TEST_TOLERANCE);
}

test "Restricted Kohn-Sham DFT with meta-GGA (SCAN) Functional on Water (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "mgga_x_scan",
            .correlation = "mgga_c_scan",
        },
        .gradient = null,
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.30199624497159, hf.energy[0], TEST_TOLERANCE);
}

test "Generalized Kohn-Sham DFT with meta-GGA (SCAN) Functional on Water (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "mgga_x_scan",
            .correlation = "mgga_c_scan",
        },
        .generalized = true,
        .gradient = null,
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.3019962449716, hf.energy[0], TEST_TOLERANCE);
}

test "Restricted Kohn-Sham DFT with Hybrid GGA (B3LYP) Functional on Water (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange_correlation = "hyb_gga_xc_b3lyp",
        },
        .gradient = null,
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.32009990219437, hf.energy[0], TEST_TOLERANCE);
}

test "Generalized Kohn-Sham DFT with Hybrid GGA (B3LYP) Functional on Water (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange_correlation = "hyb_gga_xc_b3lyp",
        },
        .generalized = true,
        .gradient = null,
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.32009990219433, hf.energy[0], TEST_TOLERANCE);
}
