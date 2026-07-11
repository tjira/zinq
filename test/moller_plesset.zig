const std = @import("std");
const zinq = @import("zinq");

const TEST_TOLERANCE = 1e-8;

test "Restricted MP2 on Water with Analytical Gradient (STO-3G)" {
    const opt = zinq.MollerPlessetOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
        },
        .order = 2,
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.moller_plesset_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0338613302971151,  0.0117296326295283,  0.0135435816121511,
        -0.0195068891589080,  0.0075152397256791, -0.0119223300334759,
        -0.0143544411382030, -0.0192448723551982, -0.0016212515786747,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.0048550569061500, res.energy[0], TEST_TOLERANCE);

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

test "Generalized MP2 on Water with Analytical Gradient (STO-3G)" {
    const opt = zinq.MollerPlessetOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
            .generalized = true,
        },
        .order = 2,
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.moller_plesset_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0338613302577487,  0.0117296326231899,  0.0135435815942967,
        -0.0195068891388165,  0.0075152397035204, -0.0119223300170307,
        -0.0143544411189378, -0.0192448723267111, -0.0016212515772662,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.0048550569062000, res.energy[0], TEST_TOLERANCE);

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

test "Restricted MP3 on Water with Analytical Gradient (STO-3G)" {
    const opt = zinq.MollerPlessetOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
        },
        .order = 3,
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.moller_plesset_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0435344876800621,  0.0150808935630933,  0.0174124402917781,
        -0.0253267507773931,  0.0109477345767154, -0.0158229658344019,
        -0.0182077368946440, -0.0260286281367898, -0.0015894744447869,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.0154555119455900, res.energy[0], TEST_TOLERANCE);

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

test "Generalized MP3 on Water with Analytical Gradient (STO-3G)" {
    const opt = zinq.MollerPlessetOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
            .generalized = true,
        },
        .order = 3,
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.moller_plesset_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0435344876800520,  0.0150808935630996,  0.0174124402917784,
        -0.0253267507773966,  0.0109477345767133, -0.0158229658343919,
        -0.0182077368995235, -0.0260286281327251, -0.0015894744447874,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.0154555119455300, res.energy[0], TEST_TOLERANCE);

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

test "Restricted MP4 on Water with Analytical Gradient (STO-3G)" {
    const opt = zinq.MollerPlessetOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
        },
        .order = 4,
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.moller_plesset_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0467931751650953,  0.0162098685616673,  0.0187157778766019,
        -0.0273563596661054,  0.0124626194249257, -0.0172750206199374,
        -0.0194368154901675, -0.0286724879832730, -0.0014407572432108,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.0187039890576300, res.energy[0], TEST_TOLERANCE);

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

test "Generalized MP4 on Water with Analytical Gradient (STO-3G)" {
    const opt = zinq.MollerPlessetOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
            .generalized = true,
        },
        .order = 4,
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.moller_plesset_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0467931751650850,  0.0162098685616736,  0.0187157778766021,
        -0.0273563596661088,  0.0124626194249235, -0.0172750206199273,
        -0.0194368154958567, -0.0286724879788217, -0.0014407572432113,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.0187039890575700, res.energy[0], TEST_TOLERANCE);

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

test "Restricted MP2 on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.MollerPlessetOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
        },
        .order = 2,
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.moller_plesset_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0338644397857024,  0.0117304537639029,  0.0135461810657489,
        -0.0195097676680689,  0.0075144001243643, -0.0119243821927739,
        -0.0143567952193280, -0.0192460682058027, -0.0016223197008003,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.0048550569062100, res.energy[0], TEST_TOLERANCE);

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

test "Generalized MP2 on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.MollerPlessetOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
            .generalized = true,
        },
        .order = 2,
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.moller_plesset_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0338644454700443,  0.0117172881175520,  0.0135484945928965,
        -0.0195056593099707,  0.0075213982597688, -0.0119255922470529,
        -0.0143508238181766, -0.0192413814659176, -0.0016228284493991,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.0048550569061500, res.energy[0], TEST_TOLERANCE);

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

test "Restricted MP3 on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.MollerPlessetOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
        },
        .order = 3,
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.moller_plesset_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0435384414743112,  0.0150821115596500,  0.0174157847254719,
        -0.0253304520470010,  0.0109466050446372, -0.0158256163729220,
        -0.0182107520174668, -0.0260302449817118, -0.0015908234729523,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.0154555119455900, res.energy[0], TEST_TOLERANCE);

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

test "Generalized MP3 on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.MollerPlessetOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
            .generalized = true,
        },
        .order = 3,
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.moller_plesset_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0435384471586531,  0.0150648055807778,  0.0174178865108843,
        -0.0253245737269481,  0.0109559174177321, -0.0158266423966325,
        -0.0182026184347706, -0.0260239403360174, -0.0015908739214865,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.0154555119455300, res.energy[0], TEST_TOLERANCE);

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

test "Restricted MP4 on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.MollerPlessetOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
        },
        .order = 4,
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.moller_plesset_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0467974231810331,  0.0162112414159310,  0.0187193855083478,
        -0.0273603497191743,  0.0124613862340084, -0.0172778818807728,
        -0.0194400634256908, -0.0286742611876889, -0.0014422013805415,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.0187039890576300, res.energy[0], TEST_TOLERANCE);

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

test "Generalized MP4 on Water with Numerical Gradient (STO-3G)" {
    const opt = zinq.MollerPlessetOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
            .generalized = true,
        },
        .order = 4,
        .gradient = .{ .numeric = .{} },
    };

    var res = try zinq.moller_plesset_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0467974288653750,  0.0161924504027411,  0.0187213139213327,
        -0.0273537786199540,  0.0124715299421041, -0.0172787800067908,
        -0.0194311169821049, -0.0286673540017546, -0.0014420024285755,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-75.0187039890575700, res.energy[0], TEST_TOLERANCE);

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
