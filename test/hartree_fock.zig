const std = @import("std");
const zinq = @import("zinq");

const TEST_TOLERANCE = 1e-8;

test "Restricted Hartree-Fock on Water with Analytical Gradient (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
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
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
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
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "lda_x",
            .correlation = "lda_c_vwn",
        },
        .gradient = null,
    };

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-74.7399385581956500, res.energy[0], TEST_TOLERANCE);
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

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-74.7399385581956500, res.energy[0], TEST_TOLERANCE);
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

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.23434062179214, res.energy[0], TEST_TOLERANCE);
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

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.23434062179214, res.energy[0], TEST_TOLERANCE);
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

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.33570008867521, res.energy[0], TEST_TOLERANCE);
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

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.33570008867527, res.energy[0], TEST_TOLERANCE);
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

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.30199624497159, res.energy[0], TEST_TOLERANCE);
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

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.3019962449716, res.energy[0], TEST_TOLERANCE);
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

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.32009990219437, res.energy[0], TEST_TOLERANCE);
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

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.32009990219433, res.energy[0], TEST_TOLERANCE);
}

test "Restricted Hartree-Fock on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0000002188471626, -0.0000014907186596,  0.0000005321965091,
        -0.0000007219114195,  0.0000006266986929, -0.0000005421441074,
         0.0000005044853424,  0.0000008640199667,  0.0000000063948846,
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
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0000002138733635, -0.0000014857448605,  0.0000005336175946,
        -0.0000007283063042,  0.0000006345146630, -0.0000005371703082,
         0.0000005044853424,  0.0000008675726804, -0.0000000028421709,
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
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "lda_x",
            .correlation = "lda_c_vwn",
        },
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0504874797968569,  0.0174896769067345,  0.0201944786226704,
        -0.0298922572028459,  0.0154032449017905, -0.0193921820823562,
        -0.0205952417786648, -0.0328929154136404, -0.0008022851716305,
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

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0504874876128269,  0.0174896740645636,  0.0201944743594140,
        -0.0298922586239314,  0.0154032399279913, -0.0193921856350698,
        -0.0205952545684340, -0.0328929104398412, -0.0008022894348869,
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
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "gga_x_pbe",
            .correlation = "gga_c_pbe",
        },
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0584284023830151,  0.0202408827476575,  0.0233700404805859,
        -0.0344129944096494,  0.0168858953486506, -0.0220802739647752,
        -0.0240154136577075, -0.0371267873333636, -0.0012897714896098,
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

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0584283995408441,  0.0202408834582002,  0.0233700426122141,
        -0.0344130015150768,  0.0168858917959369, -0.0220802796491171,
        -0.0240154115260793, -0.0371267901755346, -0.0012897700685244,
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
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "mgga_x_tpss",
            .correlation = "mgga_c_tpss",
        },
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0600934065175807,  0.0208183436711806,  0.0240339431911707,
        -0.0353932321672801,  0.0173579806528323, -0.0227054556489747,
        -0.0247001715081296, -0.0381763236134702, -0.0013284839894823,
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

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0600934065175807,  0.0208183465133516,  0.0240339424806280,
        -0.0353932328778228,  0.0173579806528323, -0.0227054506751756,
        -0.0247001644027023, -0.0381763200607566, -0.0013284839894823,
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
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "mgga_x_scan",
            .correlation = "mgga_c_scan",
        },
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0446175619117639,  0.0155010454250259,  0.0178289624841455,
        -0.0260442483579482,  0.0114168713594154, -0.0162821422122761,
        -0.0185733227908713, -0.0269179274425824, -0.0015468216929548,
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

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0446175683066485,  0.0155010482671969,  0.0178289624841455,
        -0.0260442448052345,  0.0114168692277872, -0.0162821422122761,
        -0.0185733249224995, -0.0269179288636678, -0.0015468280878395,
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
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange_correlation = "hyb_gga_xc_b3lyp",
        },
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0468003840126130,  0.0162124173641587,  0.0187193172962452,
        -0.0277282936167467,  0.0143766264670830, -0.0180137817551440,
        -0.0190720918169518, -0.0305890502261263, -0.0007055447781568,
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
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .dft = .{
            .exchange_correlation = "hyb_gga_xc_b3lyp",
        },
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0468003811704420,  0.0162124194957869,  0.0187193229805871,
        -0.0277282921956612,  0.0143766243354548, -0.0180137753602594,
        -0.0190720911064091, -0.0305890509366691, -0.0007055433570713,
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
