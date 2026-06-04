// zig fmt: off
const std  = @import("std" );
const zinq = @import("zinq");
// zig fmt: on

const TEST_TOLERANCE = 1e-8;

test "RTP on Tully's First Potential" {
    const opt = zinq.QuantumDynamicsOptions{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 2048 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 1, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .tully_1 = .{} },
        .fft = .{ .plan = .estimate },
        .mass = 2000,
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);

    const output = try zinq.quantum_dynamics_run(f64, std.testing.io, opt, false, std.testing.allocator, arena.allocator());

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.items[0].pos.?.at(0), 13.4112520392165720, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].mom.?.at(0), 16.0102008489768330, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].pop.?.at(0),  0.4103832688920894, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].pop.?.at(1),  0.5896167311083441, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].epot.?,       0.0017923343653092, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].ekin.?,       0.0647076640849498, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].norm.?,       1.0000000000004314, TEST_TOLERANCE);
    // zig fmt: on

    arena.deinit();
}

test "RTP on Time-Linear Potential" {
    const opt = zinq.QuantumDynamicsOptions{
        .grid = .{ .bounds = &.{.{ -16, 16 }}, .npoint = 2048 },
        .initial_conditions = .{ .momentum = &.{0}, .position = &.{0}, .state = 1, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .time_linear = .{} },
        .fft = .{ .plan = .estimate },
        .mass = 1,
        .iterations = 20000,
        .time_step = 0.001,
        .adiabatic = true,
    };

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);

    const output = try zinq.quantum_dynamics_run(f64, std.testing.io, opt, false, std.testing.allocator, arena.allocator());

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.items[0].pos.?.at(0), -0.0200063015629076, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].mom.?.at(0),  0.0000000000002490, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].pop.?.at(0),  0.2846073338491507, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].pop.?.at(1),  0.7153926661529834, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].epot.?,      43.0871480756309200, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].ekin.?,       0.5000000000009120, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].norm.?,       1.0000000000021340, TEST_TOLERANCE);
    // zig fmt: on

    arena.deinit();
}

test "ITP on 1D HO Potential" {
    const opt = zinq.QuantumDynamicsOptions{
        .grid = .{ .bounds = &.{.{ -8, 8 }}, .npoint = 256 },
        .initial_conditions = .{ .momentum = &.{0}, .position = &.{1}, .gamma = &.{2} },
        .potential = .{ .harmonic = .{} },
        .imaginary = .{ .nstate = 2 },
        .fft = .{ .plan = .estimate },
        .mass = 1,
        .iterations = 1000,
        .time_step = 0.01,
    };

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);

    const output = try zinq.quantum_dynamics_run(f64, std.testing.io, opt, false, std.testing.allocator, arena.allocator());

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.items[0].pos.?.at(0), 0.0000605355096514, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].mom.?.at(0), 0.0000000000000011, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].pop.?.at(0), 1.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].epot.?,      0.2499968765473224, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].ekin.?,      0.2500031253240217, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].norm.?,      1.0000000000000000, TEST_TOLERANCE);
    // zig fmt: on

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.items[1].pos.?.at(0), -0.0000226986259582, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[1].mom.?.at(0), -0.0000000000000019, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[1].pop.?.at(0),  0.9999999999999990, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[1].epot.?,       0.7499906237514752, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[1].ekin.?,       0.7500093748913410, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[1].norm.?,       0.9999999999999990, TEST_TOLERANCE);
    // zig fmt: on

    arena.deinit();
}

test "ITP on 2D HO Potential" {
    const opt = zinq.QuantumDynamicsOptions{
        .grid = .{ .bounds = &.{ .{ -8, 8 }, .{ -8, 8 } }, .npoint = 128 },
        .initial_conditions = .{ .momentum = &.{ 0, 0 }, .position = &.{ 1, 1 }, .gamma = &.{ 2, 2 } },
        .potential = .{ .harmonic = .{ .k = &.{ 1, 1 } } },
        .imaginary = .{ .nstate = 2 },
        .fft = .{ .plan = .estimate },
        .mass = 1,
        .iterations = 1000,
        .time_step = 0.01,
    };

    var arena = std.heap.ArenaAllocator.init(std.testing.allocator);

    const output = try zinq.quantum_dynamics_run(f64, std.testing.io, opt, false, std.testing.allocator, arena.allocator());

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.items[0].pos.?.at(0),  0.0000605355096522, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].mom.?.at(0), -0.0000000000000028, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].pop.?.at(0),  1.0000000000000042, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].epot.?,       0.4999937530946449, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].ekin.?,       0.5000062506480457, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[0].norm.?,       1.0000000000000042, TEST_TOLERANCE);
    // zig fmt: on

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.items[1].pos.?.at(0), -0.0000113479399309, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[1].mom.?.at(0), -0.0000000000000007, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[1].pop.?.at(0),  1.0000000000000029, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[1].epot.?,       0.9999874983520464, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[1].ekin.?,       1.0000124994137747, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.items[1].norm.?,       1.0000000000000029, TEST_TOLERANCE);
    // zig fmt: on

    arena.deinit();
}
