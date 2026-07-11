const std = @import("std");
const zinq = @import("zinq");

const TEST_TOLERANCE = 1e-8;

test "Restricted CIS on Water with Analytical Gradient (STO-3G)" {
    const opt = zinq.ConfigurationInteractionOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
        },
        .excitations = &.{1},
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.configuration_interaction_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0000002177925702, -0.0000014890587308,  0.0000005387456539,
        -0.0000007176513073,  0.0000006272955517, -0.0000005398894336,
         0.0000004998587402,  0.0000008617631853,  0.0000000011437817,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.9659012172972900, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5822026772564000, res.energy[1], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5822026772563200, res.energy[2], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5822026772562900, res.energy[3], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5070515329961900, res.energy[4], TEST_TOLERANCE);

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

test "Generalized CIS on Water with Analytical Gradient (STO-3G)" {
    const opt = zinq.ConfigurationInteractionOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
            .generalized = true,
        },
        .excitations = &.{1},
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.configuration_interaction_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0000002177925587, -0.0000014890587346,  0.0000005387456521,
        -0.0000007176513022,  0.0000006272955480, -0.0000005398894318,
         0.0000004998587437,  0.0000008617631866,  0.0000000011437803,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.9659012172973100, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5822026772563600, res.energy[1], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5822026772563500, res.energy[2], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5822026772562800, res.energy[3], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5070515329962000, res.energy[4], TEST_TOLERANCE);

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

test "Restricted CISD on Water with Analytical Gradient (STO-3G)" {
    const opt = zinq.ConfigurationInteractionOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
            .diis = null,
        },
        .excitations = &.{ 1, 2 },
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.configuration_interaction_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0470446555188239,  0.0162970018327955,  0.0188163571895635,
        -0.0276008508277128,  0.0130359747527521, -0.0175627729362029,
        -0.0194438046909178, -0.0293329765854609, -0.0012535842532794,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.0195290856878400, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6224197038007400, res.energy[1], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6224197038007200, res.energy[2], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6224197038007200, res.energy[3], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5642383868465500, res.energy[4], TEST_TOLERANCE);

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

test "Generalized CISD on Water with Analytical Gradient (STO-3G)" {
    const opt = zinq.ConfigurationInteractionOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
            .generalized = true,
            .diis = null,
        },
        .excitations = &.{ 1, 2 },
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.configuration_interaction_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0470446554994146,  0.0162970018260687,  0.0188163571818012,
        -0.0276008508198990,  0.0130359747659291, -0.0175627729360961,
        -0.0194438046796355, -0.0293329765918243, -0.0012535842456223,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.0195290856886600, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6224197033478200, res.energy[1], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6224197033477600, res.energy[2], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6224197033477000, res.energy[3], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5642383862388600, res.energy[4], TEST_TOLERANCE);

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

test "Restricted CIS on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.ConfigurationInteractionOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
        },
        .excitations = &.{1},
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.configuration_interaction_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0000002238209618, -0.0000014892975742,  0.0000005393019364,
        -0.0000007197797913,  0.0000006217248938, -0.0000005378808510,
         0.0000005037747997,  0.0000008654410522, -0.0000000007105427,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.9659012172973000, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5822026772563200, res.energy[1], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5822026772562800, res.energy[2], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5822026772562500, res.energy[3], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5070515329962400, res.energy[4], TEST_TOLERANCE);

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

test "Generalized CIS on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.ConfigurationInteractionOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
            .generalized = true,
        },
        .excitations = &.{1},
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.configuration_interaction_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0000002209787908, -0.0000014900081169,  0.0000005364597655,
        -0.0000007126743640,  0.0000006259881502, -0.0000005393019364,
         0.0000004909850304,  0.0000008760991932,  0.0000000000000000,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.9659012172972600, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5822026772563200, res.energy[1], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5822026772563000, res.energy[2], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5822026772562200, res.energy[3], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5070515329962000, res.energy[4], TEST_TOLERANCE);

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

test "Restricted CISD on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.ConfigurationInteractionOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
            .diis = null,
        },
        .excitations = &.{ 1, 2 },
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.configuration_interaction_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0470446593681117,  0.0162969996608808,  0.0188163504333261,
        -0.0276008449873188,  0.0130359715910799, -0.0175627711485049,
        -0.0194438023015664, -0.0293329740941317, -0.0012535892324195,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.0195290856878500, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6224197038008200, res.energy[1], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6224197038007400, res.energy[2], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6224197038007400, res.energy[3], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5642383868466000, res.energy[4], TEST_TOLERANCE);

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

test "Generalized CISD on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.ConfigurationInteractionOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
            .generalized = true,
            .diis = null,
        },
        .excitations = &.{ 1, 2 },
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.configuration_interaction_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0470446693157100,  0.0162970025030518,  0.0188163525649543,
        -0.0276008513822035,  0.0130359730121654, -0.0175627661747058,
        -0.0194438023015664, -0.0293329712519608, -0.0012535821269921,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.0195290856886500, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6224197033478300, res.energy[1], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6224197033477300, res.energy[2], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6224197033476600, res.energy[3], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5642383862388800, res.energy[4], TEST_TOLERANCE);

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
