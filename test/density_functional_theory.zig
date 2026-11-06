//! Validates ground-state Kohn-Sham density functional theory electronic energies for local density approximations.

const std = @import("std");
const zinq = @import("zinq");

const TEST_TOLERANCE = 1e-8;

test "Restricted LDA on Water (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .dft = .{
            .exchange = "lda_x",
            .correlation = "lda_c_vwn",
        },
    };

    var res = try zinq.electronic_structure.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-74.7399385581954800, res.energy[0], TEST_TOLERANCE);
}

test "Restricted LDA on Water with Integral-Direct Evaluation (STO-3G)" {
    const opt = zinq.hartree_fock.Options{
        .system = "example/molecule/water.xyz",
        .basis = "builtin:sto-3g",
        .integral_direct = true,
        .dft = .{
            .exchange = "lda_x",
            .correlation = "lda_c_vwn",
        },
    };

    var res = try zinq.electronic_structure.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-74.7399385581954200, res.energy[0], TEST_TOLERANCE);
}
