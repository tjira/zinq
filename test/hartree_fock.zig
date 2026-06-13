const std = @import("std");
const zinq = @import("zinq");

const TEST_TOLERANCE = 1e-12;

test "Restricted Hartree-Fock on Water (STO-3G)" {
    const opt = zinq.HartreeFockOptions{
        .system = "example/molecule/water.xyz",
        .basis = "example/basis/sto-3g.g94",
    };

    var ints = try zinq.hartree_fock_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer ints.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-74.96590121729727, ints.energy, 1e-12);
}
