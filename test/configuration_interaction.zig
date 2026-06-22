const std = @import("std");
const zinq = @import("zinq");

const TEST_TOLERANCE = 1e-8;

test "Restricted CIS on Water (STO-3G)" {
    const opt = zinq.ConfigurationInteractionOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "example/basis/sto-3g.g94",
        },
        .excitations = &.{1},
    };

    var res = try zinq.configuration_interaction_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-74.9659012172972900, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5822026772564000, res.energy[1], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5822026772563200, res.energy[2], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5822026772562900, res.energy[3], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5070515329961900, res.energy[4], TEST_TOLERANCE);
}

test "Generalized CIS on Water (STO-3G)" {
    const opt = zinq.ConfigurationInteractionOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "example/basis/sto-3g.g94",
            .generalized = true,
        },
        .excitations = &.{1},
    };

    var res = try zinq.configuration_interaction_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-74.9659012172973100, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5822026772563600, res.energy[1], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5822026772563500, res.energy[2], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5822026772562800, res.energy[3], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5070515329962000, res.energy[4], TEST_TOLERANCE);
}

test "Restricted CISD on Water (STO-3G)" {
    const opt = zinq.ConfigurationInteractionOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "example/basis/sto-3g.g94",
            .diis = null,
        },
        .excitations = &.{ 1, 2 },
    };

    var res = try zinq.configuration_interaction_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.0195290856878400, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6224197038007400, res.energy[1], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6224197038007200, res.energy[2], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6224197038007200, res.energy[3], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5642383868465500, res.energy[4], TEST_TOLERANCE);
}

test "Generalized CISD on Water (STO-3G)" {
    const opt = zinq.ConfigurationInteractionOptions{
        .hartree_fock = .{
            .system = "example/molecule/water.xyz",
            .basis = "example/basis/sto-3g.g94",
            .generalized = true,
            .diis = null,
        },
        .excitations = &.{ 1, 2 },
    };

    var res = try zinq.configuration_interaction_run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer res.deinit(std.testing.allocator);

    try std.testing.expectApproxEqAbs(-75.0195290856886600, res.energy[0], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6224197033478200, res.energy[1], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6224197033477600, res.energy[2], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.6224197033477000, res.energy[3], TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(-74.5642383862388600, res.energy[4], TEST_TOLERANCE);
}
