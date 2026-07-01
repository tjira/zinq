const std = @import("std");
const zinq = @import("zinq");

const TEST_TOLERANCE = 1e-8;

test "Restricted Hartree-Fock on Water with Analytical Gradient (STO-3G)" {
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

test "Generalized Hartree-Fock on Water with Analytical Gradient (STO-3G)" {
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

test "Restricted Hartree-Fock on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .gradient = .{ .numeric = .{} },
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0000002188471626, -0.0000014885870314,  0.0000005393019364,
        -0.0000007204903341,  0.0000006288303211, -0.0000005400124792,
         0.0000005009326287,  0.0000008590461675, -0.0000000028421709,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.9659012172972400, hf.energy[0], TEST_TOLERANCE);

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

test "Generalized Hartree-Fock on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .gradient = .{ .numeric = .{} },
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0000002145839062, -0.0000014928502878,  0.0000005421441074,
        -0.0000007126743640,  0.0000006238565220, -0.0000005364597655,
         0.0000005023537142,  0.0000008519407402,  0.0000000007105427,
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

test "Restricted Kohn-Sham DFT with LDA (VWN5) Functional on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "lda_x",
            .correlation = "lda_c_vwn",
        },
        .gradient = .{ .numeric = .{} },
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0504705710113740,  0.0174838149291645,  0.0201877078609414,
        -0.0298822371291862,  0.0153980849404434, -0.0193856969588069,
        -0.0205883374349014, -0.0328818948958087, -0.0008020194286473,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.7399385581955600, hf.energy[0], TEST_TOLERANCE);

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

test "Generalized Kohn-Sham DFT with LDA (VWN5) Functional on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .dft = .{
            .exchange = "lda_x",
            .correlation = "lda_c_vwn",
        },
        .gradient = .{ .numeric = .{} },
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0504705774062586,  0.0174838099553654,  0.0201877057293132,
        -0.0298822413924427,  0.0153980813877297, -0.0193856948271787,
        -0.0205883267767604, -0.0328818956063515, -0.0008020222708183,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.7399385581955000, hf.energy[0], TEST_TOLERANCE);

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

test "Restricted Kohn-Sham DFT with GGA (PBE) Functional on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "gga_x_pbe",
            .correlation = "gga_c_pbe",
        },
        .gradient = .{ .numeric = .{} },
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0584088503785551,  0.0202341183808130,  0.0233622174050652,
        -0.0344014694064754,  0.0168802529287859, -0.0220728750832677,
        -0.0240073667612251, -0.0371143535460305, -0.0012893451639684,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.2343406217922000, hf.energy[0], TEST_TOLERANCE);

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

test "Generalized Kohn-Sham DFT with GGA (PBE) Functional on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .dft = .{
            .exchange = "gga_x_pbe",
            .correlation = "gga_c_pbe",
        },
        .gradient = .{ .numeric = .{} },
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0584088517996406,  0.0202341162491848,  0.0233622166945224,
        -0.0344014708275608,  0.0168802543498714, -0.0220728800570669,
        -0.0240073674717678, -0.0371143528354878, -0.0012893437428829,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.2343406217922400, hf.energy[0], TEST_TOLERANCE);

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

test "Restricted Kohn-Sham DFT with meta-GGA (TPSS) Functional on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "mgga_x_tpss",
            .correlation = "mgga_c_tpss",
        },
        .gradient = .{ .numeric = .{} },
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0600732938949022,  0.0208113682731437,  0.0240258835049190,
        -0.0353813831566185,  0.0173521669921684, -0.0226978414730183,
        -0.0246919100277410, -0.0381635310020556, -0.0013280356370160,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.3357000886753000, hf.energy[0], TEST_TOLERANCE);

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

test "Generalized Kohn-Sham DFT with meta-GGA (TPSS) Functional on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .dft = .{
            .exchange = "mgga_x_tpss",
            .correlation = "mgga_c_tpss",
        },
        .gradient = .{ .numeric = .{} },
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0600732860789321,  0.0208113711153146,  0.0240258884787181,
        -0.0353813852882468,  0.0173521748081384, -0.0226978428941038,
        -0.0246919093171982, -0.0381635324231411, -0.0013280349264733,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.3357000886754100, hf.energy[0], TEST_TOLERANCE);

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

test "Restricted Kohn-Sham DFT with meta-GGA (SCAN) Functional on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "mgga_x_scan",
            .correlation = "mgga_c_scan",
        },
        .gradient = .{ .numeric = .{} },
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0446026206191163,  0.0154958492259993,  0.0178230003200497,
        -0.0260355200509821,  0.0114130585870953, -0.0162766909284073,
        -0.0185671005681343, -0.0269089149185220, -0.0015463072600141,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.3019962449717000, hf.energy[0], TEST_TOLERANCE);

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

test "Generalized Kohn-Sham DFT with meta-GGA (SCAN) Functional on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .dft = .{
            .exchange = "mgga_x_scan",
            .correlation = "mgga_c_scan",
        },
        .gradient = .{ .numeric = .{} },
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0446026220402018,  0.0154958442522002,  0.0178229981884215,
        -0.0260355236036958,  0.0114130628503517, -0.0162766887967791,
        -0.0185671076735616, -0.0269089127868938, -0.0015463086810996,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.3019962449717500, hf.energy[0], TEST_TOLERANCE);

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

test "Restricted Kohn-Sham DFT with Hybrid GGA (B3LYP) Functional on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange_correlation = "hyb_gga_xc_b3lyp",
        },
        .gradient = .{ .numeric = .{} },
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0467847087293194,  0.0162069952125421,  0.0187130503093158,
        -0.0277190054021048,  0.0143718075662491, -0.0180077435629755,
        -0.0190657075904710, -0.0305788034893339, -0.0007053017725411,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.3200999021943400, hf.energy[0], TEST_TOLERANCE);

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

test "Generalized Kohn-Sham DFT with Hybrid GGA (B3LYP) Functional on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .dft = .{
            .exchange_correlation = "hyb_gga_xc_b3lyp",
        },
        .gradient = .{ .numeric = .{} },
    };

    var hf = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer hf.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0467847030449775,  0.0162069937914566,  0.0187130545725722,
        -0.0277190004283057,  0.0143718168033047, -0.0180077456946037,
        -0.0190657132748129, -0.0305788056209622, -0.0007053003514557,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.3200999021943000, hf.energy[0], TEST_TOLERANCE);

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
