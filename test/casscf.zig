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
         0.0003680201154452,  0.0001259449883660,  0.0001476416756960,
        -0.0002498870943641,  0.0002759871083754, -0.0002046069694993,
        -0.0001181330210927, -0.0004019320967471,  0.0000569652938030,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.9671422267358900, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5823356594350700, res.energy[1], TEST_TOLERANCE);

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
         0.0000000071154264, -0.0000000053290852, -0.0000000073865912,
         0.0000000022172157,  0.0000000019238401, -0.0000000030282720,
        -0.0000000033168700,  0.0000000079157723,  0.0000000006134596,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.9671422267358900, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5823356594350700, res.energy[1], TEST_TOLERANCE);

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
