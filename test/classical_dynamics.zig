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
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 11.3312896120968050, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 13.5971917050606630, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.5300000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.4700000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      -0.0006001827820393, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0471450477819875, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Landau--Zener on Tully's First Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_1 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .landau_zener = .{} } },
        .mass = &.{2000},
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 11.6323986087966680, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 14.0180950224960980, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.3300000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.6700000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      -0.0033999979510066, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0498692969839902, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Landau--Zener on Tully's Second Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_2 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .landau_zener = .{} } },
        .mass = &.{2000},
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 12.8036285444219670, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 14.2573240059221970, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.9200000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.0800000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0040021822029992, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0525425178725186, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Landau--Zener on Tully's Second Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_2 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .landau_zener = .{} } },
        .mass = &.{2000},
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 12.5626221661956060, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 14.1973190568882540, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.9200000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.0800000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0039968873400233, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0522379841290380, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Landau--Zener on Tully's Third Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_3 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .landau_zener = .{} } },
        .mass = &.{2000},
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 27.6826368484236400, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 31.9935479192494720, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  1.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      -0.2000008999363736, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.2559452886436811, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Landau--Zener on Tully's Third Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_3 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .landau_zener = .{} } },
        .mass = &.{2000},
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 12.5220947185096630, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 15.0091217734431250, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  1.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0006000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0565447217553644, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Landau--Zener on Time-Linear Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{0}, .momentum = &.{0}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .time_linear = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .landau_zener = .{} } },
        .mass = &.{1},
        .trajectories = 100,
        .iterations = 2000,
        .time_step = 0.01,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), -0.8878446095195264, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), -0.0988554174204397, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.3100000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.6900000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      38.0075992401519700, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       1.6941659415406718, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Landau--Zener on Time-Linear Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{0}, .momentum = &.{0}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .time_linear = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .landau_zener = .{} } },
        .mass = &.{1},
        .trajectories = 100,
        .iterations = 2000,
        .time_step = 0.01,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0),  0.1908475272076625, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0),  0.0091217734431337, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.7200000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.2800000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      44.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.4526169090816229, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Fewest Switches on Tully's First Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_1 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .fewest_switches = .{} } },
        .mass = &.{2000},
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 11.5803442185294190, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 13.9584695001026870, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.6500000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.3500000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      -0.0029999904136419, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0495447657622636, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Fewest Switches on Tully's First Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_1 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .fewest_switches = .{} } },
        .mass = &.{2000},
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 11.2240753628929290, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 13.5200327736925910, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.3800000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.6200000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      -0.0023999550866957, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0464549357190951, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Fewest Switches on Tully's Second Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_2 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .fewest_switches = .{} } },
        .mass = &.{2000},
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 12.4182303937708480, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 13.7647161981900310, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.8800000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.1200000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0059957358450081, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0505489932626763, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Fewest Switches on Tully's Second Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_2 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .fewest_switches = .{} } },
        .mass = &.{2000},
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 11.5148707789799050, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 12.6126032704844920, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.9000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.1000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0032035221290807, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0472945844383648, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Fewest Switches on Tully's Third Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_3 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .fewest_switches = .{} } },
        .mass = &.{2000},
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 12.5289939577292980, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 13.3625346365611270, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.7400000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.2600000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      -0.1199315061175614, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.1758759550682386, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Fewest Switches on Tully's Third Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_3 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .fewest_switches = .{} } },
        .mass = &.{2000},
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 12.5220947185096630, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 15.0091217734431250, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.5500000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.4500000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0000600000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0565447217553644, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Fewest Switches on Time-Linear Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{0}, .momentum = &.{0}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .time_linear = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .fewest_switches = .{} } },
        .mass = &.{1},
        .trajectories = 100,
        .iterations = 2000,
        .time_step = 0.01,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0),  1.6763424299984810, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0),  0.1635017328919781, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.3600000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.6400000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      28.0055994401119800, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       2.4083820506614500, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Fewest Switches on Time-Linear Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{0}, .momentum = &.{0}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .time_linear = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .fewest_switches = .{} } },
        .mass = &.{1},
        .trajectories = 100,
        .iterations = 2000,
        .time_step = 0.01,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0),  0.1908475272076625, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0),  0.0091217734431337, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.7800000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.2200000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      56.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.4526169090816229, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Ehrenfest on Tully's First Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_1 = .{} },
        .nonadiabatic = .{ .ehrenfest = .{} },
        .mass = &.{2000},
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 11.6802042578144580, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 14.0982468202290420, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.6647537702384002, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.3352462297615999, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      -0.0032950679439719, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0499006402144838, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Ehrenfest on Tully's First Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_1 = .{} },
        .nonadiabatic = .{ .ehrenfest = .{} },
        .mass = &.{2000},
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 11.6802042578144580, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 14.0982468202290420, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.3352462297615999, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.6647537702384002, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      -0.0032950679439719, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0499006402144838, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Ehrenfest on Tully's Second Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_2 = .{} },
        .nonadiabatic = .{ .ehrenfest = .{} },
        .mass = &.{2000},
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 12.5712447447864760, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 13.9874375088259220, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.8481717205288644, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.1518282794703638, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0075913139300451, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0492208778278723, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Ehrenfest on Tully's Second Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_2 = .{} },
        .nonadiabatic = .{ .ehrenfest = .{} },
        .mass = &.{2000},
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 12.5711515914182850, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 13.9873060948859570, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.8482161667419055, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.1517838332573222, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0075923546826893, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0492198809756093, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Ehrenfest on Tully's Third Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_3 = .{} },
        .nonadiabatic = .{ .ehrenfest = .{} },
        .mass = &.{2000},
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 17.8174027527849040, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 21.2281186407351600, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.6418720815441132, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.3581279051462304, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      -0.0567490723865715, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.1126957002466461, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Ehrenfest on Tully's Third Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 0 },
        .potential = .{ .tully_3 = .{} },
        .nonadiabatic = .{ .ehrenfest = .{} },
        .mass = &.{2000},
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 4.0059001662911500, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 1.8447064645109301, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0), 0.5496561089271500, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1), 0.4503438824941452, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      0.0453131113438013, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,      0.0118295743021424, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Ehrenfest on Time-Linear Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{0}, .momentum = &.{0}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .time_linear = .{} },
        .nonadiabatic = .{ .ehrenfest = .{} },
        .mass = &.{1},
        .trajectories = 100,
        .iterations = 2000,
        .time_step = 0.01,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0),  0.1908475272076625, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0),  0.0091217734431337, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.2845956092560016, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.7153646808023937, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      43.0855216747040760, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.4526169090816228, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Ehrenfest on Time-Linear Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{0}, .momentum = &.{0}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .time_linear = .{} },
        .nonadiabatic = .{ .ehrenfest = .{} },
        .mass = &.{1},
        .trajectories = 100,
        .iterations = 2000,
        .time_step = 0.01,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0),  0.1908475272076625, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0),  0.0091217734431337, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.7245718462070273, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.2753884438513730, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      44.3968192572406650, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.4526169090816228, TEST_TOLERANCE);
    // zig fmt: on
}
