// zig fmt: off
const std  = @import("std" );
const zinq = @import("zinq");
// zig fmt: on

const TEST_TOLERANCE = 1e-8;

test "Adiabatic Landau--Zener on Tully's First Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .tully_1 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .landau_zener = .{} } },
        .mass = 2000,
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 13.5131323261878540, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 16.1427158175045270, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.4600000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.5400000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0007999958390219, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0657445231400523, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Landau--Zener on Tully's First Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .tully_1 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .landau_zener = .{} } },
        .mass = 2000,
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 13.6655776269540750, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 16.3260684484533800, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.4700000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.5300000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,      -0.0005999997545603, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0672098969720091, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Landau--Zener on Time-Linear Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{0}, .momentum = &.{0}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .time_linear = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .landau_zener = .{} } },
        .mass = 1,
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
        .mass = 1,
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
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .tully_1 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .fewest_switches = .{} } },
        .mass = 2000,
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 13.4849803222532180, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 16.0917854719325900, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.4400000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.5600000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0011999958406232, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0653446603416648, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Fewest Switches on Tully's First Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .tully_1 = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .fewest_switches = .{} } },
        .mass = 2000,
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 13.7651613977406570, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 16.4323403120651750, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.6700000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.3300000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0034000002314440, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0680141209505341, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Fewest Switches on Time-Linear Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{0}, .momentum = &.{0}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .time_linear = .{} },
        .nonadiabatic = .{ .surface_hopping = .{ .fewest_switches = .{} } },
        .mass = 1,
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
        .mass = 1,
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
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .tully_1 = .{} },
        .nonadiabatic = .{ .ehrenfest = .{} },
        .mass = 2000,
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 13.4470357090589620, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 16.0371507259824180, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.4024981905838427, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.5975018094161568, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0019500354483309, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0645290837637736, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic Ehrenfest on Tully's First Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .tully_1 = .{} },
        .nonadiabatic = .{ .ehrenfest = .{} },
        .mass = 2000,
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.classical_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.pos.?.at(0), 13.4470357090589620, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.mom.?.at(0), 16.0371507259824180, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(0),  0.5975018094161568, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.pop.?.at(1),  0.4024981905838427, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.epot.?,       0.0019500354483309, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.ekin.?,       0.0645290837637736, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic Ehrenfest on Time-Linear Potential" {
    const opt = zinq.classical_dynamics.Options{
        .initial_conditions = .{ .position = &.{0}, .momentum = &.{0}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .time_linear = .{} },
        .nonadiabatic = .{ .ehrenfest = .{} },
        .mass = 1,
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
        .mass = 1,
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
