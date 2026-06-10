// zig fmt: off
const std  = @import("std" );
const zinq = @import("zinq");
// zig fmt: on

const TEST_TOLERANCE = 1e-12;

test "Adiabatic Landau--Zener on Tully's First Potential" {
    const opt = zinq.ClassicalDynamicsOptions{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .tully_1 = .{} },
        .surface_hopping = .{ .landau_zener = .{} },
        .mass = 2000,
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);

    const output = try zinq.classical_dynamics_run(f64, std.testing.io, opt, false, std.testing.allocator, arena.allocator());

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.pos.?.at(0), 13.5131323261878540, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.mom.?.at(0), 16.1427158175045270, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.pop.?.at(0),  0.4600000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.pop.?.at(1),  0.5400000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.epot.?,       0.0007999958390219, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.ekin.?,       0.0657445231400523, TEST_TOLERANCE);
    // zig fmt: on

    arena.deinit();
}

test "Diabatic Landau--Zener on Tully's First Potential" {
    const opt = zinq.ClassicalDynamicsOptions{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .tully_1 = .{} },
        .surface_hopping = .{ .landau_zener = .{} },
        .mass = 2000,
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);

    const output = try zinq.classical_dynamics_run(f64, std.testing.io, opt, false, std.testing.allocator, arena.allocator());

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.pos.?.at(0), 13.6655776269540750, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.mom.?.at(0), 16.3260684484533800, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.pop.?.at(0),  0.4700000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.pop.?.at(1),  0.5300000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.epot.?,      -0.0005999997545603, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.ekin.?,       0.0672098969720091, TEST_TOLERANCE);
    // zig fmt: on

    arena.deinit();
}

test "Adiabatic Landau--Zener on Time-Linear Potential" {
    const opt = zinq.ClassicalDynamicsOptions{
        .initial_conditions = .{ .position = &.{0}, .momentum = &.{0}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .time_linear = .{} },
        .surface_hopping = .{ .landau_zener = .{} },
        .mass = 1,
        .trajectories = 100,
        .iterations = 2000,
        .time_step = 0.01,
        .adiabatic = true,
    };

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);

    const output = try zinq.classical_dynamics_run(f64, std.testing.io, opt, false, std.testing.allocator, arena.allocator());

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.pos.?.at(0), -0.8878446095195264, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.mom.?.at(0), -0.0988554174204397, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.pop.?.at(0),  0.3100000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.pop.?.at(1),  0.6900000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.epot.?,      38.0075992401519700, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.ekin.?,       1.6941659415406718, TEST_TOLERANCE);
    // zig fmt: on

    arena.deinit();
}

test "Diabatic Landau--Zener on Time-Linear Potential" {
    const opt = zinq.ClassicalDynamicsOptions{
        .initial_conditions = .{ .position = &.{0}, .momentum = &.{0}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .time_linear = .{} },
        .surface_hopping = .{ .landau_zener = .{} },
        .mass = 1,
        .trajectories = 100,
        .iterations = 2000,
        .time_step = 0.01,
        .adiabatic = false,
    };

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);

    const output = try zinq.classical_dynamics_run(f64, std.testing.io, opt, false, std.testing.allocator, arena.allocator());

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.pos.?.at(0),  0.1908475272076625, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.mom.?.at(0),  0.0091217734431337, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.pop.?.at(0),  0.7200000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.pop.?.at(1),  0.2800000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.epot.?,      44.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.ekin.?,       0.4526169090816229, TEST_TOLERANCE);
    // zig fmt: on

    arena.deinit();
}

test "Adiabatic Fewest Switches on Tully's First Potential" {
    const opt = zinq.ClassicalDynamicsOptions{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .tully_1 = .{} },
        .surface_hopping = .{ .fewest_switches = .{} },
        .mass = 2000,
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);

    const output = try zinq.classical_dynamics_run(f64, std.testing.io, opt, false, std.testing.allocator, arena.allocator());

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.pos.?.at(0), 13.4849803222532180, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.mom.?.at(0), 16.0917854719325900, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.pop.?.at(0),  0.4400000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.pop.?.at(1),  0.5600000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.epot.?,       0.0011999958406232, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.ekin.?,       0.0653446603416648, TEST_TOLERANCE);
    // zig fmt: on

    arena.deinit();
}

test "Diabatic Fewest Switches on Tully's First Potential" {
    const opt = zinq.ClassicalDynamicsOptions{
        .initial_conditions = .{ .position = &.{-10}, .momentum = &.{15}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .tully_1 = .{} },
        .surface_hopping = .{ .fewest_switches = .{} },
        .mass = 2000,
        .trajectories = 100,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);

    const output = try zinq.classical_dynamics_run(f64, std.testing.io, opt, false, std.testing.allocator, arena.allocator());

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.pos.?.at(0), 13.7651613977406570, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.mom.?.at(0), 16.4323403120651750, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.pop.?.at(0),  0.6700000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.pop.?.at(1),  0.3300000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.epot.?,       0.0034000002314440, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.ekin.?,       0.0680141209505341, TEST_TOLERANCE);
    // zig fmt: on

    arena.deinit();
}

test "Adiabatic Fewest Switches on Time-Linear Potential" {
    const opt = zinq.ClassicalDynamicsOptions{
        .initial_conditions = .{ .position = &.{0}, .momentum = &.{0}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .time_linear = .{} },
        .surface_hopping = .{ .fewest_switches = .{} },
        .mass = 1,
        .trajectories = 100,
        .iterations = 2000,
        .time_step = 0.01,
        .adiabatic = true,
    };

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);

    const output = try zinq.classical_dynamics_run(f64, std.testing.io, opt, false, std.testing.allocator, arena.allocator());

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.pos.?.at(0),  1.6763424299984810, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.mom.?.at(0),  0.1635017328919781, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.pop.?.at(0),  0.3600000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.pop.?.at(1),  0.6400000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.epot.?,      28.0055994401119800, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.ekin.?,       2.4083820506614500, TEST_TOLERANCE);
    // zig fmt: on

    arena.deinit();
}

test "Diabatic Fewest Switches on Time-Linear Potential" {
    const opt = zinq.ClassicalDynamicsOptions{
        .initial_conditions = .{ .position = &.{0}, .momentum = &.{0}, .gamma = &.{2}, .state = 1 },
        .potential = .{ .time_linear = .{} },
        .surface_hopping = .{ .fewest_switches = .{} },
        .mass = 1,
        .trajectories = 100,
        .iterations = 2000,
        .time_step = 0.01,
        .adiabatic = false,
    };

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);

    const output = try zinq.classical_dynamics_run(f64, std.testing.io, opt, false, std.testing.allocator, arena.allocator());

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.pos.?.at(0),  0.1908475272076625, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.mom.?.at(0),  0.0091217734431337, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.pop.?.at(0),  0.7800000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.pop.?.at(1),  0.2200000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.epot.?,      56.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.ekin.?,       0.4526169090816229, TEST_TOLERANCE);
    // zig fmt: on

    arena.deinit();
}
