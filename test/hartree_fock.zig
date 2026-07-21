const std = @import("std");
const zinq = @import("zinq");

const TEST_TOLERANCE = 1e-8;

test "Restricted Hartree-Fock on Water with Analytical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0000002178410248, -0.0000014890494290,  0.0000005387672022,
        -0.0000007176761931,  0.0000006273176250, -0.0000005399082537,
         0.0000004998351610,  0.0000008617318090,  0.0000000011410497,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.9659012172971900, res.energy[0], TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(expected_grad[0], res.grad[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], res.grad[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], res.grad[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], res.grad[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], res.grad[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], res.grad[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], res.grad[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], res.grad[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], res.grad[0].at(2, 2), TEST_TOLERANCE);
}

test "Generalized Hartree-Fock on Water with Analytical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0000002178410003, -0.0000014890494329,  0.0000005387672015,
        -0.0000007176761749,  0.0000006273176159, -0.0000005399082350,
         0.0000004998351675,  0.0000008617318237,  0.0000000011410507,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.9659012172972300, res.energy[0], TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(expected_grad[0], res.grad[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], res.grad[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], res.grad[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], res.grad[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], res.grad[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], res.grad[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], res.grad[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], res.grad[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], res.grad[0].at(2, 2), TEST_TOLERANCE);
}

test "Restricted Kohn-Sham DFT with LDA (VWN5) Functional on Water (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "lda_x",
            .correlation = "lda_c_vwn",
        },
        .gradient = null,
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-74.7399385581956500, res.energy[0], TEST_TOLERANCE);
}

test "Generalized Kohn-Sham DFT with LDA (VWN5) Functional on Water (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "lda_x",
            .correlation = "lda_c_vwn",
        },
        .generalized = true,
        .gradient = null,
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-74.7399385581956500, res.energy[0], TEST_TOLERANCE);
}

test "Restricted Kohn-Sham DFT with GGA (PBE) Functional on Water (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "gga_x_pbe",
            .correlation = "gga_c_pbe",
        },
        .gradient = null,
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.23434062179214, res.energy[0], TEST_TOLERANCE);
}

test "Generalized Kohn-Sham DFT with GGA (PBE) Functional on Water (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "gga_x_pbe",
            .correlation = "gga_c_pbe",
        },
        .generalized = true,
        .gradient = null,
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.23434062179214, res.energy[0], TEST_TOLERANCE);
}

test "Restricted Kohn-Sham DFT with meta-GGA (TPSS) Functional on Water (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "mgga_x_tpss",
            .correlation = "mgga_c_tpss",
        },
        .gradient = null,
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.33570008867521, res.energy[0], TEST_TOLERANCE);
}

test "Generalized Kohn-Sham DFT with meta-GGA (TPSS) Functional on Water (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "mgga_x_tpss",
            .correlation = "mgga_c_tpss",
        },
        .generalized = true,
        .gradient = null,
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.33570008867527, res.energy[0], TEST_TOLERANCE);
}

test "Restricted Kohn-Sham DFT with meta-GGA (SCAN) Functional on Water (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "mgga_x_scan",
            .correlation = "mgga_c_scan",
        },
        .gradient = null,
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.30199624497159, res.energy[0], TEST_TOLERANCE);
}

test "Generalized Kohn-Sham DFT with meta-GGA (SCAN) Functional on Water (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "mgga_x_scan",
            .correlation = "mgga_c_scan",
        },
        .generalized = true,
        .gradient = null,
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.3019962449716, res.energy[0], TEST_TOLERANCE);
}

test "Restricted Kohn-Sham DFT with Hybrid GGA (B3LYP) Functional on Water (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange_correlation = "hyb_gga_xc_b3lyp",
        },
        .gradient = null,
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.32009990219437, res.energy[0], TEST_TOLERANCE);
}

test "Generalized Kohn-Sham DFT with Hybrid GGA (B3LYP) Functional on Water (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange_correlation = "hyb_gga_xc_b3lyp",
        },
        .generalized = true,
        .gradient = null,
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.32009990219433, res.energy[0], TEST_TOLERANCE);
}

test "Restricted Hartree-Fock on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0000002138733635, -0.0000014864554032,  0.0000005435651929,
        -0.0000007190692486,  0.0000006274092357, -0.0000005378808510,
         0.0000005037747997,  0.0000008640199667,  0.0000000021316282,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.9659012172972400, res.energy[0], TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(expected_grad[0], res.grad[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], res.grad[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], res.grad[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], res.grad[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], res.grad[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], res.grad[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], res.grad[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], res.grad[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], res.grad[0].at(2, 2), TEST_TOLERANCE);
}

test "Generalized Hartree-Fock on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0000002195577053, -0.0000014878764887,  0.0000005393019364,
        -0.0000007233325050,  0.0000006238565220, -0.0000005350386800,
         0.0000005051958851,  0.0000008633094239, -0.0000000035527137,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.9659012172971900, res.energy[0], TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(expected_grad[0], res.grad[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], res.grad[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], res.grad[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], res.grad[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], res.grad[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], res.grad[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], res.grad[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], res.grad[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], res.grad[0].at(2, 2), TEST_TOLERANCE);
}

test "Restricted Kohn-Sham DFT with LDA (VWN5) Functional on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "lda_x",
            .correlation = "lda_c_vwn",
        },
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0504874833495705,  0.0174896662485935,  0.0201944729383285,
        -0.0298922536501323,  0.0154032370858204, -0.0193921955826681,
        -0.0205952410681221, -0.0328929097292985, -0.0008022837505450,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.7399385581955600, res.energy[0], TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(expected_grad[0], res.grad[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], res.grad[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], res.grad[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], res.grad[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], res.grad[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], res.grad[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], res.grad[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], res.grad[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], res.grad[0].at(2, 2), TEST_TOLERANCE);
}

test "Generalized Kohn-Sham DFT with LDA (VWN5) Functional on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .dft = .{
            .exchange = "lda_x",
            .correlation = "lda_c_vwn",
        },
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0504874854811987,  0.0174896634064226,  0.0201944665434439,
        -0.0298922564923032,  0.0154032377963631, -0.0193921977142963,
        -0.0205952353837802, -0.0328929097292985, -0.0008022780662031,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.7399385581955000, res.energy[0], TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(expected_grad[0], res.grad[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], res.grad[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], res.grad[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], res.grad[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], res.grad[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], res.grad[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], res.grad[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], res.grad[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], res.grad[0].at(2, 2), TEST_TOLERANCE);
}

test "Restricted Kohn-Sham DFT with GGA (PBE) Functional on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "gga_x_pbe",
            .correlation = "gga_c_pbe",
        },
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0584284151727843,  0.0202408919847130,  0.0233700383489577,
        -0.0344129986729058,  0.0168859003224497, -0.0220802768069461,
        -0.0240154150787930, -0.0371267958598764, -0.0012897622525543,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.2343406217922000, res.energy[0], TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(expected_grad[0], res.grad[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], res.grad[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], res.grad[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], res.grad[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], res.grad[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], res.grad[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], res.grad[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], res.grad[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], res.grad[0].at(2, 2), TEST_TOLERANCE);
}

test "Generalized Kohn-Sham DFT with GGA (PBE) Functional on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .dft = .{
            .exchange = "gga_x_pbe",
            .correlation = "gga_c_pbe",
        },
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0584284173044125,  0.0202408848792857,  0.0233700390595004,
        -0.0344130015150768,  0.0168859003224497, -0.0220802725436897,
        -0.0240154122366221, -0.0371267951493337, -0.0012897665158107,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.2343406217922400, res.energy[0], TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(expected_grad[0], res.grad[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], res.grad[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], res.grad[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], res.grad[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], res.grad[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], res.grad[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], res.grad[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], res.grad[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], res.grad[0].at(2, 2), TEST_TOLERANCE);
}

test "Restricted Kohn-Sham DFT with meta-GGA (TPSS) Functional on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "mgga_x_tpss",
            .correlation = "mgga_c_tpss",
        },
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0600934107808371,  0.0208183351446678,  0.0240339282697732,
        -0.0353932449570493,  0.0173579849160888, -0.0227054535173465,
        -0.0247001850084416, -0.0381763221923848, -0.0013284804367686,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.3357000886753000, res.energy[0], TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(expected_grad[0], res.grad[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], res.grad[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], res.grad[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], res.grad[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], res.grad[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], res.grad[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], res.grad[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], res.grad[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], res.grad[0].at(2, 2), TEST_TOLERANCE);
}

test "Generalized Kohn-Sham DFT with meta-GGA (TPSS) Functional on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .dft = .{
            .exchange = "mgga_x_tpss",
            .correlation = "mgga_c_tpss",
        },
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0600934093597516,  0.0208183379868387,  0.0240339325330297,
        -0.0353932428254211,  0.0173579778106614, -0.0227054592016884,
        -0.0247001857189844, -0.0381763257450984, -0.0013284854105677,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.3357000886754100, res.energy[0], TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(expected_grad[0], res.grad[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], res.grad[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], res.grad[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], res.grad[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], res.grad[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], res.grad[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], res.grad[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], res.grad[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], res.grad[0].at(2, 2), TEST_TOLERANCE);
}

test "Restricted Kohn-Sham DFT with meta-GGA (SCAN) Functional on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "mgga_x_scan",
            .correlation = "mgga_c_scan",
        },
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0446175661750203,  0.0155010482671969,  0.0178289710106583,
        -0.0260442448052345,  0.0114168855702701, -0.0162821514493317,
        -0.0185733298962987, -0.0269179274425824, -0.0015468117453565,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.3019962449717000, res.energy[0], TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(expected_grad[0], res.grad[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], res.grad[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], res.grad[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], res.grad[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], res.grad[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], res.grad[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], res.grad[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], res.grad[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], res.grad[0].at(2, 2), TEST_TOLERANCE);
}

test "Generalized Kohn-Sham DFT with meta-GGA (SCAN) Functional on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .dft = .{
            .exchange = "mgga_x_scan",
            .correlation = "mgga_c_scan",
        },
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0446175697277340,  0.0155010482671969,  0.0178289703001155,
        -0.0260442440946917,  0.0114168791753855, -0.0162821478966180,
        -0.0185733320279269, -0.0269179295742106, -0.0015468060610146,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.3019962449717500, res.energy[0], TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(expected_grad[0], res.grad[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], res.grad[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], res.grad[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], res.grad[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], res.grad[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], res.grad[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], res.grad[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], res.grad[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], res.grad[0].at(2, 2), TEST_TOLERANCE);
}

test "Restricted Kohn-Sham DFT with Hybrid GGA (B3LYP) Functional on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange_correlation = "hyb_gga_xc_b3lyp",
        },
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0468003861442412,  0.0162124187852442,  0.0187193258227580,
        -0.0277282929062039,  0.0143766300197967, -0.0180137874394859,
        -0.0190720868431526, -0.0305890488050409, -0.0007055433570713,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.3200999021943400, res.energy[0], TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(expected_grad[0], res.grad[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], res.grad[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], res.grad[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], res.grad[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], res.grad[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], res.grad[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], res.grad[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], res.grad[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], res.grad[0].at(2, 2), TEST_TOLERANCE);
}

test "Generalized Kohn-Sham DFT with Hybrid GGA (B3LYP) Functional on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .dft = .{
            .exchange_correlation = "hyb_gga_xc_b3lyp",
        },
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0468003911180403,  0.0162124173641587,  0.0187193265333008,
        -0.0277282957483749,  0.0143766378357668, -0.0180137831762295,
        -0.0190720982118364, -0.0305890516472118, -0.0007055369621867,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.3200999021943000, res.energy[0], TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(expected_grad[0], res.grad[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], res.grad[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], res.grad[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], res.grad[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], res.grad[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], res.grad[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], res.grad[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], res.grad[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], res.grad[0].at(2, 2), TEST_TOLERANCE);
}
