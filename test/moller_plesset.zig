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
        .gradient = true,
    };

    var res = try zinq.moller_plesset_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0338122794133439,  0.0117126467140770,  0.0135239610741029,
        -0.0195946586449228,  0.0081070984905233, -0.0121370677177705,
        -0.0142176207684284, -0.0198197452045953, -0.0013868933563342,
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
        .gradient = true,
    };

    var res = try zinq.moller_plesset_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
         0.0338122794133231,  0.0117126467140746,  0.0135239610741036,
        -0.0195946586449066,  0.0081070984905145, -0.0121370677177529,
        -0.0142176207684237, -0.0198197452045824, -0.0013868933563335,
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
