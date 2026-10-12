// zig fmt: off
const std  = @import("std" );
const zinq = @import("zinq");
// zig fmt: on

const TEST_TOLERANCE = 1e-8;

test "Adiabatic Landau--Zener on Tully's First Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_1 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .landau_zener = .{} } },
        .mass = &.{2000},
        .trajectories = 10,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 10.6656268035362718, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 13.1521976635855786, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.5000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.5000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      -0.0000000017158227, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0439197242839718, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Landau--Zener on Tully's First Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_1 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .landau_zener = .{} } },
        .mass = &.{2000},
        .trajectories = 10,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 10.4102050656993903, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 12.7282698987647933, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.6000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.4000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0019998953645219, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0418823550035157, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Landau--Zener on Tully's Second Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_2 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .landau_zener = .{} } },
        .mass = &.{2000},
        .trajectories = 10,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 12.7927660164116048, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 14.6354249793054265, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  1.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      -0.0000000196917084, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0539195441611981, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Landau--Zener on Tully's Second Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_2 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .landau_zener = .{} } },
        .mass = &.{2000},
        .trajectories = 10,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 12.5004176736399657, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 14.5927482922727982, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  1.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0536100151621981, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Landau--Zener on Tully's Third Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_3 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .landau_zener = .{} } },
        .mass = &.{2000},
        .trajectories = 10,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 26.9222894860173838, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 31.8271555187349904, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  1.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      -0.2000008999539862, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.2533197213564702, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Landau--Zener on Tully's Third Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_3 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .landau_zener = .{} } },
        .mass = &.{2000},
        .trajectories = 10,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 11.9222407208208416, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 14.6354276821528622, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  1.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0006000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0539195684542596, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Landau--Zener on Time-Linear Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{0}, .momentum = &.{0}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .time_linear = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .landau_zener = .{} } },
        .mass = &.{1},
        .trajectories = 10,
        .iterations = 2000,
        .time_step = 0.01,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), -10.933662718651488, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), -0.7260653668261500, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.2000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.8000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      60.0119988002399367, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       1.6087210520063333, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Landau--Zener on Time-Linear Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{0}, .momentum = &.{0}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .time_linear = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .landau_zener = .{} } },
        .mass = &.{1},
        .trajectories = 10,
        .iterations = 2000,
        .time_step = 0.01,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), -7.3223471593512031, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), -0.3645723178471355, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.7000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.3000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      40.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.8077216762263024, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Fewest Switches on Tully's First Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_1 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .fewest_switches = .{} } },
        .mass = &.{2000},
        .trajectories = 10,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 11.0050887212904378, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 13.7769409725990499, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.7000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.3000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      -0.0039999961708517, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0479196129299625, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Fewest Switches on Tully's First Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_1 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .fewest_switches = .{} } },
        .mass = &.{2000},
        .trajectories = 10,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 11.4462517262895229, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 14.0675748664262166, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.1000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.9000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      -0.0079999981562800, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0508921161863623, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Fewest Switches on Tully's Second Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_2 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .fewest_switches = .{} } },
        .mass = &.{2000},
        .trajectories = 10,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 12.6758586536157249, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 14.6354351149612150, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  1.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      -0.0000000427207902, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0539196098911001, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Fewest Switches on Tully's Second Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_2 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .fewest_switches = .{} } },
        .mass = &.{2000},
        .trajectories = 10,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0),  9.8316493874426403, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 10.7497391119596024, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.9000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.1000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0049999999999527, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0440895818063596, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Fewest Switches on Tully's Third Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_3 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .fewest_switches = .{} } },
        .mass = &.{2000},
        .trajectories = 10,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 15.1345197513123360, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 17.6684970542176103, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.9000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.1000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      -0.1400606365636651, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.1933795329088452, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Fewest Switches on Tully's Third Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_3 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .fewest_switches = .{} } },
        .mass = &.{2000},
        .trajectories = 10,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 11.9222407208208416, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 14.6354276821528622, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.5000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.5000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0539195684542596, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Fewest Switches on Time-Linear Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{0}, .momentum = &.{0}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .time_linear = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .fewest_switches = .{} } },
        .mass = &.{1},
        .trajectories = 10,
        .iterations = 2000,
        .time_step = 0.01,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), -3.4436523009974729, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0),  0.0548681083088737, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.3000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.7000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      40.0079992001599578, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       3.0243765838306489, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Fewest Switches on Time-Linear Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{0}, .momentum = &.{0}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .time_linear = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .fewest_switches = .{} } },
        .mass = &.{1},
        .trajectories = 10,
        .iterations = 2000,
        .time_step = 0.01,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0),  -7.3223471593512031, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0),  -0.3645723178471355, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),   1.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),   0.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      100.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,        0.8077216762263024, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Ehrenfest on Tully's First Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_1 = .{} },
        .nonadiabatic = .{ .ehrenfest = .{} },
        .mass = &.{2000},
        .trajectories = 10,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 11.1050744301319391, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 13.7392283350592006, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.6778064362111490, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.3221935637888466, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      -0.0035561249904751, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0475356050762153, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Ehrenfest on Tully's First Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_1 = .{} },
        .nonadiabatic = .{ .ehrenfest = .{} },
        .mass = &.{2000},
        .trajectories = 10,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 11.1050744301319355, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 13.7392283350592006, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.3221935637888462, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.6778064362111489, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      -0.0035561249904751, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0475356050762153, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Ehrenfest on Tully's Second Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_2 = .{} },
        .nonadiabatic = .{ .ehrenfest = .{} },
        .mass = &.{2000},
        .trajectories = 10,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 12.0237319718250149, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 13.5994442246659926, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.8505797952958785, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.1494202047033702, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0074709545291668, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0467114470472358, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Ehrenfest on Tully's Second Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_2 = .{} },
        .nonadiabatic = .{ .ehrenfest = .{} },
        .mass = &.{2000},
        .trajectories = 10,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 12.0235650532486886, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 13.5993029433022237, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.8503511626097490, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.1496488373895047, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0074723223268263, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0467100979913131, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Ehrenfest on Tully's Third Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_3 = .{} },
        .nonadiabatic = .{ .ehrenfest = .{} },
        .mass = &.{2000},
        .trajectories = 10,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 17.2706908790505942, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 21.0981204670464066, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.6450608001291388, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.3549391869530560, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      -0.0580245580115124, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.1113458759664723, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Ehrenfest on Tully's Third Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_3 = .{} },
        .nonadiabatic = .{ .ehrenfest = .{} },
        .mass = &.{2000},
        .trajectories = 10,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0),  2.5501452145589107, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), -1.4160861088221239, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.3801717018587997, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.6198282912718285, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0355182152434366, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0189996160465501, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Ehrenfest on Time-Linear Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{0}, .momentum = &.{0}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .time_linear = .{} },
        .nonadiabatic = .{ .ehrenfest = .{} },
        .mass = &.{1},
        .trajectories = 10,
        .iterations = 2000,
        .time_step = 0.01,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), -7.3223471593512031, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), -0.3645723178471355, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.2845956092560151, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.7153646808023986, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      43.0855216747033865, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.8077216762263024, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Ehrenfest on Time-Linear Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{0}, .momentum = &.{0}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .time_linear = .{} },
        .nonadiabatic = .{ .ehrenfest = .{} },
        .mass = &.{1},
        .trajectories = 10,
        .iterations = 2000,
        .time_step = 0.01,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), -7.3223471593512031, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), -0.3645723178471355, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.7245718462070273, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.2753884438513731, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      44.3968192572406792, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.8077216762263024, TEST_TOLERANCE);
    // zig fmt: on
}
