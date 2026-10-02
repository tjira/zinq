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

test "Restricted Hartree-Fock on Water with Integral-Direct Analytical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .integral_direct = true,
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

test "Generalized Hartree-Fock on Water with Integral-Direct Analytical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .integral_direct = true,
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

test "Restricted Hartree-Fock on Water with Integral-Direct SCF and Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .integral_direct = true,
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

test "Generalized Hartree-Fock on Water with Integral-Direct SCF and Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .integral_direct = true,
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
         0.0504874961393398,  0.0174896761961918,  0.0201944679645294,
        -0.0298922493868758,  0.0154032370858204, -0.0193921962932109,
        -0.0205952318310665, -0.0328929132820122, -0.0008022809083741,
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

test "Restricted Kohn-Sham DFT with LDA (VWN5) Functional on Water with Integral-Direct SCF and Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "lda_x",
            .correlation = "lda_c_vwn",
        },
        .integral_direct = true,
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0504874961393398,  0.0174896761961918,  0.0201944679645294,
        -0.0298922493868758,  0.0154032370858204, -0.0193921962932109,
        -0.0205952318310665, -0.0328929132820122, -0.0008022809083741,
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
         0.0504874918760834,  0.0174896747751063,  0.0201944658329012,
        -0.0298922500974186,  0.0154032434807050, -0.0193921955826681,
        -0.0205952318310665, -0.0328929090187557, -0.0008022858821732,
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

test "Generalized Kohn-Sham DFT with LDA (VWN5) Functional on Water with Integral-Direct SCF and Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .dft = .{
            .exchange = "lda_x",
            .correlation = "lda_c_vwn",
        },
        .integral_direct = true,
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0504874932971688,  0.0174896712223926,  0.0201944651223585,
        -0.0298922515185041,  0.0154032406385340, -0.0193921962932109,
        -0.0205952332521520, -0.0328929111503840, -0.0008022865927160,
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

test "Restricted Kohn-Sham DFT with GGA (PBE) Functional on Water with Integral-Direct SCF and Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "gga_x_pbe",
            .correlation = "gga_c_pbe",
        },
        .integral_direct = true,
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

test "Generalized Kohn-Sham DFT with GGA (PBE) Functional on Water with Integral-Direct SCF and Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .dft = .{
            .exchange = "gga_x_pbe",
            .correlation = "gga_c_pbe",
        },
        .integral_direct = true,
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
         0.0600934242811491,  0.0208183443817234,  0.0240339254276023,
        -0.0353932414043356,  0.0173579842055460, -0.0227054570700602,
        -0.0247001693765014, -0.0381763285872694, -0.0013284818578541,
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

test "Restricted Kohn-Sham DFT with meta-GGA (TPSS) Functional on Water with Integral-Direct SCF and Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "mgga_x_tpss",
            .correlation = "mgga_c_tpss",
        },
        .integral_direct = true,
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0600934242811491,  0.0208183465133516,  0.0240339254276023,
        -0.0353932442465066,  0.0173579827844605, -0.0227054563595175,
        -0.0247001693765014, -0.0381763285872694, -0.0013284818578541,
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
         0.0600934214389781,  0.0208183500660653,  0.0240339218748886,
        -0.0353932414043356,  0.0173579778106614, -0.0227054520962611,
        -0.0247001700870442, -0.0381763214818420, -0.0013284711997130,
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

test "Generalized Kohn-Sham DFT with meta-GGA (TPSS) Functional on Water with Integral-Direct SCF and Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .dft = .{
            .exchange = "mgga_x_tpss",
            .correlation = "mgga_c_tpss",
        },
        .integral_direct = true,
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0600934164651790,  0.0208183543293217,  0.0240339261381450,
        -0.0353932399832502,  0.0173579827844605, -0.0227054556489747,
        -0.0247001757713861, -0.0381763207712993, -0.0013284676469993,
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
         0.0446175789647896,  0.0155010376090559,  0.0178289567998036,
        -0.0260442519106618,  0.0114168820175564, -0.0162821372384769,
        -0.0185733377122688, -0.0269179231793260, -0.0015468096137283,
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

test "Restricted Kohn-Sham DFT with meta-GGA (SCAN) Functional on Water with Integral-Direct SCF and Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "mgga_x_scan",
            .correlation = "mgga_c_scan",
        },
        .integral_direct = true,
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0446175796753323,  0.0155010404512268,  0.0178289567998036,
        -0.0260442526212046,  0.0114168756226718, -0.0162821393701051,
        -0.0185733348700978, -0.0269179203371550, -0.0015468074821001,
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
         0.0446175832280460,  0.0155010468461114,  0.0178289596419745,
        -0.0260442490684909,  0.0114168841491846, -0.0162821372384769,
        -0.0185733391333542, -0.0269179231793260, -0.0015468145875275,
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

test "Generalized Kohn-Sham DFT with meta-GGA (SCAN) Functional on Water with Integral-Direct SCF and Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .dft = .{
            .exchange = "mgga_x_scan",
            .correlation = "mgga_c_scan",
        },
        .integral_direct = true,
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0446175803858750,  0.0155010482671969,  0.0178289589314318,
        -0.0260442462263200,  0.0114168848597274, -0.0162821358173915,
        -0.0185733412649824, -0.0269179196266123, -0.0015468117453565,
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
         0.0468003797493566,  0.0162124280222997,  0.0187193229805871,
        -0.0277282971694603,  0.0143766307303395, -0.0180137831762295,
        -0.0190720932380373, -0.0305890601737246, -0.0007055369621867,
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

test "Restricted Kohn-Sham DFT with Hybrid GGA (B3LYP) Functional on Water with Integral-Direct SCF and Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange_correlation = "hyb_gga_xc_b3lyp",
        },
        .integral_direct = true,
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0468003847231557,  0.0162124258906715,  0.0187193229805871,
        -0.0277282936167467,  0.0143766271776258, -0.0180137838867722,
        -0.0190720896853236, -0.0305890587526392, -0.0007055426465286,
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
         0.0468003847231557,  0.0162124266012142,  0.0187193165857025,
        -0.0277282985905458,  0.0143766435201087, -0.0180137838867722,
        -0.0190720868431526, -0.0305890530682973, -0.0007055398043576,
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

test "Generalized Kohn-Sham DFT with Hybrid GGA (B3LYP) Functional on Water with Integral-Direct SCF and Numerical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .dft = .{
            .exchange_correlation = "hyb_gga_xc_b3lyp",
        },
        .integral_direct = true,
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.hartree_fock.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0468003854336985,  0.0162124244695860,  0.0187193244016726,
        -0.0277282936167467,  0.0143766385463095, -0.0180137838867722,
        -0.0190720882642381, -0.0305890466734127, -0.0007055426465286,
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
