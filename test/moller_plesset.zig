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
         0.0338613268979771,  0.0117296536927824,  0.0135435655579386,
        -0.0195068821540190,  0.0075152307488224, -0.0119223237504684,
        -0.0143544433228726, -0.0192448844416049, -0.0016212453601838,
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
         0.0338613219241779,  0.0117296593771243,  0.0135435669790240,
        -0.0195068885489036,  0.0075152385647925, -0.0119223187766693,
        -0.0143544440334153, -0.0192448808888912, -0.0016212545972394,
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
         0.0435344823301875,  0.0150809249532813,  0.0174124203056181,
        -0.0253267387506639,  0.0109477220178178, -0.0158229561009193,
        -0.0182077421584381, -0.0260286469710991, -0.0015894677574124,
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
         0.0435344773563884,  0.0150809306376232,  0.0174124217267035,
        -0.0253267451455486,  0.0109477291232452, -0.0158229511271202,
        -0.0182077428689809, -0.0260286434183854, -0.0015894769944680,
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
         0.0467931691616741,  0.0162099041745023,  0.0187157567665963,
        -0.0273563458108583,  0.0124626055253430, -0.0172750091564922,
        -0.0194368219297303, -0.0286725096998452, -0.0014407504522751,
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
         0.0467931641878749,  0.0162099098588442,  0.0187157581876818,
        -0.0273563522057429,  0.0124626126307703, -0.0172750041826930,
        -0.0194368226402730, -0.0286725061471316, -0.0014407596893307,
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
