const std = @import("std");
const zinq = @import("zinq");

const TEST_TOLERANCE = 1e-8;

test "Restricted Hartree-Fock on Water with Analytical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.electronic_structure.run(f64, std.testing.io, opt, false, std.testing.allocator);
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

    var mul = try zinq.population_analysis.mulliken(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer mul.deinit(std.testing.allocator);

    var low = try zinq.population_analysis.lowdin(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer low.deinit(std.testing.allocator);

    var may = try zinq.population_analysis.mayer(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer may.deinit(std.testing.allocator);

    var wib = try zinq.population_analysis.wiberg(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer wib.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(-0.3305298848381675, mul.at(0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.1652652608338326, mul.at(1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.1652646240043356, mul.at(2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(-0.2315525783308345, low.at(0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.1157764533114921, low.at(1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.1157761250193441, low.at(2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs( 0.0000000000000000, may.at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9641856021324611, may.at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9641858126236432, may.at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9641856021324611, may.at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, may.at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0085017914290628, may.at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9641858126236432, may.at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0085017914290628, may.at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, may.at(2, 2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs( 0.0000000000000000, wib.at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9863699855869255, wib.at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9863700616038187, wib.at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9863699855869255, wib.at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, wib.at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0002258272716871, wib.at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9863700616038190, wib.at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0002258272716871, wib.at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, wib.at(2, 2), TEST_TOLERANCE);
    // zig fmt: on
}

test "Restricted Hartree-Fock on Water with Integral-Direct Analytical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .integral_direct = true,
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.electronic_structure.run(f64, std.testing.io, opt, false, std.testing.allocator);
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

    var mul = try zinq.population_analysis.mulliken(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer mul.deinit(std.testing.allocator);

    var low = try zinq.population_analysis.lowdin(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer low.deinit(std.testing.allocator);

    var may = try zinq.population_analysis.mayer(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer may.deinit(std.testing.allocator);

    var wib = try zinq.population_analysis.wiberg(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer wib.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(-0.3305298848381870, mul.at(0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.1652652608338461, mul.at(1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.1652646240043425, mul.at(2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(-0.2315525783308488, low.at(0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.1157764533115030, low.at(1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.1157761250193499, low.at(2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs( 0.0000000000000000, may.at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9641856021324600, may.at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9641858126236413, may.at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9641856021324600, may.at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, may.at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0085017914290621, may.at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9641858126236413, may.at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0085017914290621, may.at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, may.at(2, 2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs( 0.0000000000000000, wib.at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9863699855869258, wib.at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9863700616038175, wib.at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9863699855869253, wib.at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, wib.at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0002258272716870, wib.at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9863700616038178, wib.at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0002258272716870, wib.at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, wib.at(2, 2), TEST_TOLERANCE);
    // zig fmt: on
}

test "Generalized Hartree-Fock on Water with Analytical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.electronic_structure.run(f64, std.testing.io, opt, false, std.testing.allocator);
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

    var mul = try zinq.population_analysis.mulliken(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer mul.deinit(std.testing.allocator);

    var low = try zinq.population_analysis.lowdin(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer low.deinit(std.testing.allocator);

    var may = try zinq.population_analysis.mayer(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer may.deinit(std.testing.allocator);

    var wib = try zinq.population_analysis.wiberg(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer wib.deinit(std.testing.allocator);

    const s2 = try zinq.spin_analysis.calculateTotalSpin(f64, res.P_ao, res.ints.S.?, std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(-0.3305298846562366, mul.at(0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.1652652607291367, mul.at(1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.1652646239270967, mul.at(2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(-0.2315525781757657, low.at(0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.1157764532222183, low.at(1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.1157761249535468, low.at(2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs( 0.0000000000000000, may.at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9641856021661122, may.at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9641858126482160, may.at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9641856021661122, may.at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, may.at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0085017914300194, may.at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9641858126482160, may.at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0085017914300194, may.at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, may.at(2, 2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs( 0.0000000000000000, wib.at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9863699856074000, wib.at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9863700616188548, wib.at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9863699856073997, wib.at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, wib.at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0002258272718870, wib.at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9863700616188548, wib.at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0002258272718870, wib.at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, wib.at(2, 2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(-0.0000000000000009, s2, TEST_TOLERANCE);
    // zig fmt: on
}

test "Generalized Hartree-Fock on Water with Integral-Direct Analytical Gradient (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .generalized = true,
        .integral_direct = true,
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.electronic_structure.run(f64, std.testing.io, opt, false, std.testing.allocator);
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

    var mul = try zinq.population_analysis.mulliken(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer mul.deinit(std.testing.allocator);

    var low = try zinq.population_analysis.lowdin(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer low.deinit(std.testing.allocator);

    var may = try zinq.population_analysis.mayer(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer may.deinit(std.testing.allocator);

    var wib = try zinq.population_analysis.wiberg(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer wib.deinit(std.testing.allocator);

    const s2 = try zinq.spin_analysis.calculateTotalSpin(f64, res.P_ao, res.ints.S.?, std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(-0.3305298846562312, mul.at(0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.1652652607291363, mul.at(1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.1652646239270970, mul.at(2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(-0.2315525781757621, low.at(0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.1157764532222181, low.at(1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.1157761249535471, low.at(2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs( 0.0000000000000000, may.at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9641856021661107, may.at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9641858126482153, may.at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9641856021661107, may.at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, may.at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0085017914300196, may.at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9641858126482153, may.at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0085017914300196, may.at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, may.at(2, 2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs( 0.0000000000000000, wib.at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9863699856073977, wib.at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9863700616188535, wib.at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9863699856073980, wib.at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, wib.at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0002258272718870, wib.at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.9863700616188538, wib.at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0002258272718870, wib.at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, wib.at(2, 2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs( 0.0000000000000000, s2, TEST_TOLERANCE);
    // zig fmt: on
}

test "Restricted Hartree-Fock on Water with RI Analytical Gradient (cc-pVDZ)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:cc-pvdz",
        .aux_basis = "builtin:cc-pvdz-rifit",
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.electronic_structure.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
        -0.0540505549919206, -0.0187267254417092, -0.0216176929728003,
         0.0320177257340895, -0.0165742310010988,  0.0207924883279960,
         0.0220328292578451,  0.0353009564428026,  0.0008252046448018,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-76.0242031339119300, res.energy[0], TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(expected_grad[0], res.grad[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], res.grad[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], res.grad[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], res.grad[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], res.grad[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], res.grad[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], res.grad[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], res.grad[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], res.grad[0].at(2, 2), TEST_TOLERANCE);

    var mul = try zinq.population_analysis.mulliken(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer mul.deinit(std.testing.allocator);

    var low = try zinq.population_analysis.lowdin(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer low.deinit(std.testing.allocator);

    var may = try zinq.population_analysis.mayer(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer may.deinit(std.testing.allocator);

    var wib = try zinq.population_analysis.wiberg(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer wib.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(-0.3326962850644861, mul.at(0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.1663478145515493, mul.at(1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.1663484705129412, mul.at(2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(-0.4972446997003832, low.at(0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.2486220445757987, low.at(1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.2486226551245864, low.at(2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs( 0.0000000000000000, may.at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 1.0098464941689850, may.at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 1.0098459189558386, may.at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 1.0098464941689846, may.at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, may.at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0060585928588789, may.at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 1.0098459189558386, may.at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0060585928588789, may.at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, may.at(2, 2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs( 0.0000000000000000, wib.at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 1.0435251369776446, wib.at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 1.0435238302424540, wib.at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 1.0435251369776450, wib.at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, wib.at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0416054940822643, wib.at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 1.0435238302424543, wib.at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0416054940822644, wib.at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, wib.at(2, 2), TEST_TOLERANCE);
    // zig fmt: on
}

test "Generalized Hartree-Fock on Water with RI Analytical Gradient (cc-pVDZ)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:cc-pvdz",
        .aux_basis = "builtin:cc-pvdz-rifit",
        .generalized = true,
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.electronic_structure.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
        -0.0540505574435590, -0.0187267262920912, -0.0216176939530583,
         0.0320177266955608, -0.0165742292061106,  0.0207924882905328,
         0.0220328307480071,  0.0353009554981973,  0.0008252056625232,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-76.0242031339116700, res.energy[0], TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(expected_grad[0], res.grad[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], res.grad[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], res.grad[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], res.grad[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], res.grad[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], res.grad[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], res.grad[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], res.grad[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], res.grad[0].at(2, 2), TEST_TOLERANCE);

    var mul = try zinq.population_analysis.mulliken(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer mul.deinit(std.testing.allocator);

    var low = try zinq.population_analysis.lowdin(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer low.deinit(std.testing.allocator);

    var may = try zinq.population_analysis.mayer(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer may.deinit(std.testing.allocator);

    var wib = try zinq.population_analysis.wiberg(f64, res.ints.sys, res.P_ao, res.ints.S.?, std.testing.allocator);
    defer wib.deinit(std.testing.allocator);

    const s2 = try zinq.spin_analysis.calculateTotalSpin(f64, res.P_ao, res.ints.S.?, std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(-0.3326962762135395, mul.at(0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.1663478101252444, mul.at(1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.1663484660882926, mul.at(2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(-0.4972446971793669, low.at(0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.2486220433147943, low.at(1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.2486226538645747, low.at(2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs( 0.0000000000000000, may.at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 1.0098464975058894, may.at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 1.0098459222922247, may.at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 1.0098464975058885, may.at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, may.at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0060585915046183, may.at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 1.0098459222922247, may.at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0060585915046183, may.at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, may.at(2, 2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs( 0.0000000000000000, wib.at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 1.0435251382670320, wib.at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 1.0435238315308506, wib.at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 1.0435251382670316, wib.at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, wib.at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0416054935122587, wib.at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 1.0435238315308510, wib.at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0416054935122588, wib.at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs( 0.0000000000000000, wib.at(2, 2), TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs( 0.0000000000000018, s2, TEST_TOLERANCE);
    // zig fmt: on
}
