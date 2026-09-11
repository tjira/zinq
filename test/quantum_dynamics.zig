// zig fmt: off
const std  = @import("std" );
const zinq = @import("zinq");
// zig fmt: on

const TEST_TOLERANCE = 1e-8;

test "Adiabatic RTP on Tully's First Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .tully_1 = .{} },
        .fft = .{ .plan = .estimate },
        .mass = &.{2000},
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.items[0].pos.?.at(0), 11.6269248968113810, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].mom.?.at(0), 14.0456023908310640, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(0),  0.6769375217267849, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(1),  0.3230624782736362, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].epot.?,      -0.0035387559533022, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].ekin.?,       0.0500387575031156, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].norm.?,       1.0000000000004217, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic RTP on Tully's First Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = false },
        .potential = .{ .tully_1 = .{} },
        .fft = .{ .plan = .estimate },
        .mass = &.{2000},
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.items[0].pos.?.at(0), 11.6269248968113480, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].mom.?.at(0), 14.0456023908310300, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(0),  0.3230626433492805, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(1),  0.6769373566511372, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].epot.?,      -0.0035387559533022, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].ekin.?,       0.0500387575031154, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].norm.?,       1.0000000000004183, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic RTP on Tully's Second Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .tully_2 = .{} },
        .fft = .{ .plan = .estimate },
        .mass = &.{2000},
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.items[0].pos.?.at(0), 12.4393051391604530, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].mom.?.at(0), 13.7764794654338230, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(0),  0.8733073000105006, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(1),  0.1266926999896673, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].epot.?,       0.0062932662885039, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].ekin.?,       0.0502066811056963, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].norm.?,       1.0000000000001665, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic RTP on Tully's Second Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = false },
        .potential = .{ .tully_2 = .{} },
        .fft = .{ .plan = .estimate },
        .mass = &.{2000},
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.items[0].pos.?.at(0), 12.4393078601794990, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].mom.?.at(0), 13.7764819035420650, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(0),  0.8739386183194002, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(1),  0.1260613816807633, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].epot.?,       0.0062933022922637, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].ekin.?,       0.0502066977305322, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].norm.?,       1.0000000000001630, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic RTP on Tully's Third Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .tully_3 = .{} },
        .fft = .{ .plan = .estimate },
        .mass = &.{2000},
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.items[0].pos.?.at(0), -12.0664752176863580, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].mom.?.at(0), -14.9298798323700700, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(0),   0.5701885511501665, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(1),   0.4298114488501240, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].epot.?,       -0.0000841350199430, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].ekin.?,        0.0559839768140773, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].norm.?,        1.0000000000002907, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic RTP on Tully's Third Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = false },
        .potential = .{ .tully_3 = .{} },
        .fft = .{ .plan = .estimate },
        .mass = &.{2000},
        .iterations = 3000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.items[0].pos.?.at(0), -11.6085920927679020, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].mom.?.at(0), -15.0674401601427230, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(0),   0.5749203324097494, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(1),   0.4250796675905397, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].epot.?,        0.0000953637213076, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].ekin.?,        0.0570046362590877, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].norm.?,        1.0000000000002878, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic RTP on Time-Linear Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -16, 16 }}, .npoint = 256 },
        .initial_conditions = .{ .momentum = &.{0}, .position = &.{0}, .state = 1, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .time_linear = .{} },
        .fft = .{ .plan = .estimate },
        .mass = &.{1},
        .iterations = 2000,
        .time_step = 0.01,
        .adiabatic = true,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.items[0].pos.?.at(0), -0.1600504125015748, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].mom.?.at(0), -0.0000000000000497, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(0),  0.2846081266515837, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(1),  0.7153918733486487, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].epot.?,      43.0869894832450560, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].ekin.?,       0.5000000000000960, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].norm.?,       1.0000000000002331, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic RTP on Time-Linear Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -16, 16 }}, .npoint = 256 },
        .initial_conditions = .{ .momentum = &.{0}, .position = &.{0}, .state = 1, .gamma = &.{2}, .adiabatic = false },
        .potential = .{ .time_linear = .{} },
        .fft = .{ .plan = .estimate },
        .mass = &.{1},
        .iterations = 2000,
        .time_step = 0.01,
        .adiabatic = false,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.items[0].pos.?.at(0), -0.1600504125015565, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].mom.?.at(0), -0.0000000000000480, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(0),  0.7327195219719144, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(1),  0.2672804780283119, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].epot.?,      44.8273168241485200, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].ekin.?,       0.5000000000000941, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].norm.?,       1.0000000000002260, TEST_TOLERANCE);
    // zig fmt: on
}

test "ITP on 1D HO Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -8, 8 }}, .npoint = 128 },
        .initial_conditions = .{ .momentum = &.{0}, .position = &.{1}, .gamma = &.{2} },
        .potential = .{ .harmonic = .{} },
        .imaginary = .{ .nstate = 2 },
        .fft = .{ .plan = .estimate },
        .mass = &.{1},
        .iterations = 1000,
        .time_step = 0.01,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.items[0].pos.?.at(0),  0.0000605355096521, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].mom.?.at(0), -0.0000000000000033, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(0),  1.0000000000000000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].epot.?,       0.2499968765473219, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].ekin.?,       0.2500031253240222, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].norm.?,       1.0000000000000000, TEST_TOLERANCE);
    // zig fmt: on

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.items[1].pos.?.at(0), -0.0000226986259616, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[1].mom.?.at(0),  0.0000000000000019, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[1].pop.?.at(0),  1.0000000000000004, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[1].epot.?,       0.7499906237514761, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[1].ekin.?,       0.7500093748913406, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[1].norm.?,       1.0000000000000004, TEST_TOLERANCE);
    // zig fmt: on
}

test "ITP on 2D HO Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{ .{ -8, 8 }, .{ -8, 8 } }, .npoint = 64 },
        .initial_conditions = .{ .momentum = &.{ 0, 0 }, .position = &.{ 1, 1 }, .gamma = &.{ 2, 2 } },
        .potential = .{ .harmonic = .{ .k = &.{ 1, 1 } } },
        .imaginary = .{ .nstate = 2 },
        .fft = .{ .plan = .estimate },
        .mass = &.{ 1, 1 },
        .iterations = 1000,
        .time_step = 0.01,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.items[0].pos.?.at(0),  0.0000605355096546, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].mom.?.at(0),  0.0000000000000002, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(0),  1.0000000000000016, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].epot.?,       0.4999937530946458, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].ekin.?,       0.5000062506480444, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].norm.?,       1.0000000000000016, TEST_TOLERANCE);
    // zig fmt: on

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.items[1].pos.?.at(0), -0.0000113479399280, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[1].mom.?.at(0),  0.0000000000000003, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[1].pop.?.at(0),  1.0000000000000016, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[1].epot.?,       0.9999874983520405, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[1].ekin.?,       1.0000124994137694, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[1].norm.?,       1.0000000000000016, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic RTP Cross Section on Tully's First Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .tully_1 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.001,
            .stop_norm = 1e-2,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.00025,
            .e_max = 0.225,
            .e_step = 0.001,
        },
        .fft = .{ .plan = .estimate },
        .mass = &.{2000},
        .iterations = 1000000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    var sum: f64 = 0;

    for (output.cross_section.?.data) |val| {
        sum += val;
    }

    // zig fmt: off
    try std.testing.expectApproxEqAbs(sum, 77.7645328896552400, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic RTP Cross Section on Tully's First Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = false },
        .potential = .{ .tully_1 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.001,
            .stop_norm = 1e-2,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.00025,
            .e_max = 0.225,
            .e_step = 0.001,
        },
        .fft = .{ .plan = .estimate },
        .mass = &.{2000},
        .iterations = 1000000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    var sum: f64 = 0;

    for (output.cross_section.?.data) |val| {
        sum += val;
    }

    // zig fmt: off
    try std.testing.expectApproxEqAbs(sum, 77.7646030591721000, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic RTP Cross Section on Tully's Second Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .tully_2 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.02,
            .stop_norm = 1e-2,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.04,
            .e_max = 2.5,
            .e_step = 0.01,
        },
        .fft = .{ .plan = .estimate },
        .mass = &.{2000},
        .iterations = 1000000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    var sum: f64 = 0;

    for (output.cross_section.?.data) |val| {
        sum += val;
    }

    // zig fmt: off
    try std.testing.expectApproxEqAbs(sum, 5.9770907740295485, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic RTP Cross Section on Tully's Second Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = false },
        .potential = .{ .tully_2 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.02,
            .stop_norm = 1e-2,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.04,
            .e_max = 2.5,
            .e_step = 0.01,
        },
        .fft = .{ .plan = .estimate },
        .mass = &.{2000},
        .iterations = 1000000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    var sum: f64 = 0;

    for (output.cross_section.?.data) |val| {
        sum += val;
    }

    // zig fmt: off
    try std.testing.expectApproxEqAbs(sum, 5.9773512785414140, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic RTP Cross Section on Tully's Third Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .tully_3 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.002,
            .stop_norm = 1e-2,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.0055,
            .e_max = 0.225,
            .e_step = 0.001,
        },
        .fft = .{ .plan = .estimate },
        .mass = &.{2000},
        .iterations = 1000000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    var sum: f64 = 0;

    for (output.cross_section.?.data) |val| {
        sum += val;
    }

    // zig fmt: off
    try std.testing.expectApproxEqAbs(sum, 75.5296966499151500, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic RTP Cross Section on Tully's Third Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = false },
        .potential = .{ .tully_3 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.002,
            .stop_norm = 1e-2,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.0055,
            .e_max = 0.225,
            .e_step = 0.001,
        },
        .fft = .{ .plan = .estimate },
        .mass = &.{2000},
        .iterations = 1000000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    var sum: f64 = 0;

    for (output.cross_section.?.data) |val| {
        sum += val;
    }

    // zig fmt: off
    try std.testing.expectApproxEqAbs(sum, 75.7432002737949300, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic RTP Partial Waves on Tully's First Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .tully_1 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.001,
            .stop_norm = 1e-2,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.00025,
            .e_max = 0.225,
            .e_step = 0.001,
        },
        .partial_waves = .{ .j_min = 0, .j_max = 10 },
        .fft = .{ .plan = .estimate },
        .mass = &.{2000},
        .iterations = 1000000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    var sum: f64 = 0;

    for (output.cross_section.?.data) |val| {
        sum += val;
    }

    // zig fmt: off
    try std.testing.expectApproxEqAbs(sum, 137.5260385613105100, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic RTP Partial Waves on Tully's First Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = false },
        .potential = .{ .tully_1 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.001,
            .stop_norm = 1e-2,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.00025,
            .e_max = 0.225,
            .e_step = 0.001,
        },
        .partial_waves = .{ .j_min = 0, .j_max = 10 },
        .fft = .{ .plan = .estimate },
        .mass = &.{2000},
        .iterations = 1000000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    var sum: f64 = 0;

    for (output.cross_section.?.data) |val| {
        sum += val;
    }

    // zig fmt: off
    try std.testing.expectApproxEqAbs(sum, 137.5116217529468700, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic RTP Partial Waves on Tully's Second Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .tully_2 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.02,
            .stop_norm = 1e-2,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.04,
            .e_max = 2.5,
            .e_step = 0.01,
        },
        .partial_waves = .{ .j_min = 0, .j_max = 10 },
        .fft = .{ .plan = .estimate },
        .mass = &.{2000},
        .iterations = 1000000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    var sum: f64 = 0;

    for (output.cross_section.?.data) |val| {
        sum += val;
    }

    // zig fmt: off
    try std.testing.expectApproxEqAbs(sum, 9.4108097828783630, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic RTP Partial Waves on Tully's Second Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = false },
        .potential = .{ .tully_2 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.02,
            .stop_norm = 1e-2,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.04,
            .e_max = 2.5,
            .e_step = 0.01,
        },
        .partial_waves = .{ .j_min = 0, .j_max = 10 },
        .fft = .{ .plan = .estimate },
        .mass = &.{2000},
        .iterations = 1000000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    var sum: f64 = 0;

    for (output.cross_section.?.data) |val| {
        sum += val;
    }

    // zig fmt: off
    try std.testing.expectApproxEqAbs(sum, 9.4103936491983500, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic RTP Partial Waves on Tully's Third Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .tully_3 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.002,
            .stop_norm = 1e-2,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.0055,
            .e_max = 0.225,
            .e_step = 0.001,
        },
        .partial_waves = .{ .j_min = 0, .j_max = 10 },
        .fft = .{ .plan = .estimate },
        .mass = &.{2000},
        .iterations = 1000000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    var sum: f64 = 0;

    for (output.cross_section.?.data) |val| {
        sum += val;
    }

    // zig fmt: off
    try std.testing.expectApproxEqAbs(sum, 135.4783472560165800, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic RTP Partial Waves on Tully's Third Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = false },
        .potential = .{ .tully_3 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.002,
            .stop_norm = 1e-2,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.0055,
            .e_max = 0.225,
            .e_step = 0.001,
        },
        .partial_waves = .{ .j_min = 0, .j_max = 10 },
        .fft = .{ .plan = .estimate },
        .mass = &.{2000},
        .iterations = 1000000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    var sum: f64 = 0;

    for (output.cross_section.?.data) |val| {
        sum += val;
    }

    // zig fmt: off
    try std.testing.expectApproxEqAbs(sum, 135.2753744156125300, TEST_TOLERANCE);
    // zig fmt: on
}
