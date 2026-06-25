const std = @import("std");
const zinq = @import("zinq");

const TEST_TOLERANCE = 1e-8;

test "Restricted MP2 on Water (STO-3G)" {
    const opt = zinq.MollerPlessetOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "example/basis/sto-3g.g94",
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

    try std.testing.expectApproxEqAbs(expected_grad[0], res.gradient[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], res.gradient[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], res.gradient[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], res.gradient[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], res.gradient[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], res.gradient[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], res.gradient[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], res.gradient[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], res.gradient[0].at(2, 2), TEST_TOLERANCE);
}

test "Generalized MP2 on Water (STO-3G)" {
    const opt = zinq.MollerPlessetOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "example/basis/sto-3g.g94",
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

    try std.testing.expectApproxEqAbs(expected_grad[0], res.gradient[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[1], res.gradient[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[2], res.gradient[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[3], res.gradient[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[4], res.gradient[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[5], res.gradient[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[6], res.gradient[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[7], res.gradient[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_grad[8], res.gradient[0].at(2, 2), TEST_TOLERANCE);
}
