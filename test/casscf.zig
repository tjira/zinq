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
        -0.0042017622262023, -0.0014578213925164, -0.0016799164817509,
         0.0024191397991716, -0.0009288866354294,  0.0014776447663356,
         0.0017826224270249,  0.0023867080279418,  0.0002022717154132,
    };
    // zig fmt: on

    try std.testing.expectApproxEqAbs(-74.9609874914254400, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6000541696338800, res.energy[1], TEST_TOLERANCE);

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
