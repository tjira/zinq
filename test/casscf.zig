//! Validates Complete Active Space Self-Consistent Field energies, analytic gradients, and nonadiabatic couplings.

const std = @import("std");
const zinq = @import("zinq");

const TEST_TOLERANCE = 1e-8;

test "CASSCF(2, 2) on Water with Analytical Gradient (STO-3G)" {
    const opt = zinq.casscf.Options{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
        },
        .active_electrons = 2,
        .active_orbitals = 2,
        .nstate = 2,
        .gradient = .{ .analytic = .{} },
    };

    var res = try zinq.casscf.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_grad = [_]f64{
        -0.0034503490701551, -0.0011984580335957, -0.0013791045175785,
         0.0020502567616785, -0.0010960352028815,  0.0013414661508448,
         0.0014000923084667,  0.0022944932364693,  0.0000376383667331,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.9607593092007400, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5994698247005200, res.energy[1], TEST_TOLERANCE);

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

test "CASSCF(2, 2) on Water with Numerical Nonadiabatic Coupling Vectors (STO-3G)" {
    const opt = zinq.casscf.Options{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "builtin:sto-3g",
        },
        .active_electrons = 2,
        .active_orbitals = 2,
        .nstate = 2,
        .nacv = .{ .numeric = .{} },
    };

    var res = try zinq.casscf.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    // zig fmt: off
    const expected_nacv = [_]f64{
        -0.0000000020580018,  0.0000000048290234,  0.0000000016759061,
         0.0000000034045267,  0.0000000005827526,  0.0000000034607395,
         0.0000000008063784, -0.0000000005184780,  0.0000000010674211,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.9607593092007400, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5994698247005200, res.energy[1], TEST_TOLERANCE);

    try std.testing.expectApproxEqAbs(expected_nacv[0], res.nacv[0].at(0, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_nacv[1], res.nacv[0].at(0, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_nacv[2], res.nacv[0].at(0, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_nacv[3], res.nacv[0].at(1, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_nacv[4], res.nacv[0].at(1, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_nacv[5], res.nacv[0].at(1, 2), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_nacv[6], res.nacv[0].at(2, 0), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_nacv[7], res.nacv[0].at(2, 1), TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(expected_nacv[8], res.nacv[0].at(2, 2), TEST_TOLERANCE);
}
