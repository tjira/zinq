const std = @import("std");
const zinq = @import("zinq");

const TEST_TOLERANCE = 1e-8;

test "Restricted CIS on Water with Analytical Gradient (STO-3G)" {
    const opt = zinq.configuration_interaction.Options{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
        },
        .space = .{ .excitations = &.{1} },
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.electronic_structure.run(f64, std.testing.io, opt, false, std.testing.allocator);
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
    const opt = zinq.configuration_interaction.Options{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
            .generalized = true,
        },
        .space = .{ .excitations = &.{1} },
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.electronic_structure.run(f64, std.testing.io, opt, false, std.testing.allocator);
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
    const opt = zinq.configuration_interaction.Options{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
            .diis = null,
        },
        .space = .{ .excitations = &.{ 1, 2 } },
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.electronic_structure.run(f64, std.testing.io, opt, false, std.testing.allocator);
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
    const opt = zinq.configuration_interaction.Options{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
            .generalized = true,
            .diis = null,
        },
        .space = .{ .excitations = &.{ 1, 2 } },
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.electronic_structure.run(f64, std.testing.io, opt, false, std.testing.allocator);
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

test "Restricted CASCI(2, 2) on Water with Analytical Gradient (STO-3G)" {
    const opt = zinq.configuration_interaction.Options{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
        },
        .space = .{
            .cas = .{
                .active_electrons = 2,
                .active_orbitals = 2,
            },
        },
        .nstate = 2,
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.electronic_structure.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
        -0.0000502940916258, -0.0000190526019582, -0.0000196455021886,
         0.0000092064654927,  0.0000889232831214, -0.0000210670833004,
         0.0000410876244518, -0.0000698706818709,  0.0000407125851463,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.9671384470571100, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5822026773010000, res.energy[1], TEST_TOLERANCE);

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

test "Restricted CIS on Water with RI Analytical Gradient (STO-3G)" {
    const opt = zinq.configuration_interaction.Options{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
            .aux_basis = "builtin:cc-pvdz-rifit",
        },
        .aux_basis = "builtin:cc-pvdz-rifit",
        .space = .{ .excitations = &.{1} },
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.electronic_structure.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
        -0.0000247323239280, -0.0000101329469743, -0.0000094402779061,
         0.0000179260986877, -0.0000270966851443,  0.0000167846808210,
         0.0000068062254714,  0.0000372296315669, -0.0000073444029505,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.9663802872875000, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5825939796854700, res.energy[1], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5825939796853700, res.energy[2], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5825939796853700, res.energy[3], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5074825820719000, res.energy[4], TEST_TOLERANCE);

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

test "Generalized CIS on Water with RI Analytical Gradient (STO-3G)" {
    const opt = zinq.configuration_interaction.Options{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
            .aux_basis = "builtin:cc-pvdz-rifit",
            .generalized = true,
        },
        .aux_basis = "builtin:cc-pvdz-rifit",
        .space = .{ .excitations = &.{1} },
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.electronic_structure.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
        -0.0000247323239249, -0.0000101329469765, -0.0000094402779085,
         0.0000179260986897, -0.0000270966851472,  0.0000167846808241,
         0.0000068062254724,  0.0000372296315629, -0.0000073444029511,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.9663802872874600, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5825939796579400, res.energy[1], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5825939796579300, res.energy[2], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5825939796578700, res.energy[3], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5074825820441200, res.energy[4], TEST_TOLERANCE);

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

test "Restricted CISD on Water with RI Analytical Gradient (STO-3G)" {
    const opt = zinq.configuration_interaction.Options{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
            .aux_basis = "builtin:cc-pvdz-rifit",
            .diis = null,
        },
        .aux_basis = "builtin:cc-pvdz-rifit",
        .space = .{ .excitations = &.{ 1, 2 } },
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.electronic_structure.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0470166048852008,  0.0162872837653640,  0.0188051380905998,
        -0.0275808928809093,  0.0130100141665803, -0.0175453003571556,
        -0.0194357120041803, -0.0292972979323280, -0.0012598377333928,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.0200022365328700, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6228130080259000, res.energy[1], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6228130080258600, res.energy[2], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6228130080258400, res.energy[3], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5646741873121400, res.energy[4], TEST_TOLERANCE);

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

test "Generalized CISD on Water with RI Analytical Gradient (STO-3G)" {
    const opt = zinq.configuration_interaction.Options{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
            .aux_basis = "builtin:cc-pvdz-rifit",
            .generalized = true,
            .diis = null,
        },
        .aux_basis = "builtin:cc-pvdz-rifit",
        .space = .{ .excitations = &.{ 1, 2 } },
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.electronic_structure.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0470166048852008,  0.0162872837653640,  0.0188051380905998,
        -0.0275808928809093,  0.0130100141665803, -0.0175453003571556,
        -0.0194357120041803, -0.0292972979323280, -0.0012598377333928,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.0200022365338000, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6228130075728700, res.energy[1], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6228130075728600, res.energy[2], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6228130075728200, res.energy[3], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5646741867043400, res.energy[4], TEST_TOLERANCE);

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

test "Restricted CASCI(2, 2) on Water with RI Analytical Gradient (STO-3G)" {
    const opt = zinq.configuration_interaction.Options{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
            .aux_basis = "builtin:cc-pvdz-rifit",
        },
        .aux_basis = "builtin:cc-pvdz-rifit",
        .space = .{
            .cas = .{
                .active_electrons = 2,
                .active_orbitals = 2,
            },
        },
        .nstate = 2,
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.electronic_structure.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
        -0.0000742766950915, -0.0000273612692497, -0.0000292375686269,
         0.0000273299060955,  0.0000612216403478, -0.0000040091005398,
         0.0000469467875536, -0.0000338603723510,  0.0000332466687910,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.9676161849660600, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5825939796852300, res.energy[1], TEST_TOLERANCE);

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
