const std = @import("std");
const zinq = @import("zinq");

const TEST_TOLERANCE = 1e-12;

test "Restricted MP2 on Water (STO-3G)" {
    const opt = zinq.MollerPlessetOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "example/basis/sto-3g.g94",
        },
        .order = 2,
    };

    var res = try zinq.moller_plesset_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.0048550571061800, res.energy[0], TEST_TOLERANCE);
}

test "Generalized MP2 on Water (STO-3G)" {
    const opt = zinq.MollerPlessetOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "example/basis/sto-3g.g94",
            .generalized = true,
        },
        .order = 2,
    };

    var res = try zinq.moller_plesset_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.0048550573543900, res.energy[0], TEST_TOLERANCE);
}
