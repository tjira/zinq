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
        .mass = &.{2000},
        .iterations = 2000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.items[0].pos.?.at(0),  4.6040561114198800, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].mom.?.at(0), 14.0462170873964460, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(0),  0.6769701557199891, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(1),  0.3230298442801738, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].epot.?,      -0.0035386177925500, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].ekin.?,       0.0500386193018480, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].norm.?,       1.0000000000001632, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic RTP on Tully's First Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = false },
        .potential = .{ .tully_1 = .{} },
        .mass = &.{2000},
        .iterations = 2000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.items[0].pos.?.at(0),  4.6040561114198780, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].mom.?.at(0), 14.0462170873964420, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(0),  0.3230054350415000, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(1),  0.6769945649586628, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].epot.?,      -0.0035386177925500, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].ekin.?,       0.0500386193018479, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].norm.?,       1.0000000000001634, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic RTP on Tully's Second Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .tully_2 = .{} },
        .mass = &.{2000},
        .iterations = 2000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.items[0].pos.?.at(0),  5.5353968492776180, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].mom.?.at(0), 13.9002720801034700, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(0),  0.8732796039592149, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(1),  0.1267203960409700, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].epot.?,       0.0059983168687150, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].ekin.?,       0.0505016317022329, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].norm.?,       1.0000000000001852, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic RTP on Tully's Second Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = false },
        .potential = .{ .tully_2 = .{} },
        .mass = &.{2000},
        .iterations = 2000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.items[0].pos.?.at(0),  5.5353983999932140, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].mom.?.at(0), 13.9002482628999040, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(0),  0.8721425547931240, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(1),  0.1278574452070559, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].epot.?,       0.0059985439593710, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].ekin.?,       0.0505014572587268, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].norm.?,       1.0000000000001796, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic RTP on Tully's Third Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .tully_3 = .{} },
        .mass = &.{2000},
        .iterations = 2000,
        .time_step = 1,
        .adiabatic = true,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.items[0].pos.?.at(0), -4.6179769489588830, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].mom.?.at(0), -14.7003583297921630, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(0),   0.5758258949699433, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(1),   0.4241741050303902, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].epot.?,        0.0014196688389071, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].ekin.?,        0.0544801693843232, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].norm.?,        1.0000000000003328, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic RTP on Tully's Third Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 512 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = false },
        .potential = .{ .tully_3 = .{} },
        .mass = &.{2000},
        .iterations = 2000,
        .time_step = 1,
        .adiabatic = false,
    };

    var output = try zinq.quantum_dynamics.run(f64, std.testing.io, opt, false, std.testing.allocator);
    defer output.deinit(std.testing.allocator);

    // zig fmt: off
    try std.testing.expectApproxEqAbs(output.observables.items[0].pos.?.at(0), -4.1180414692391670, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].mom.?.at(0), -14.5327049251101780, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(0),   0.6218105055945835, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].pop.?.at(1),   0.3781894944057727, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].epot.?,        0.0036813741032688, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].ekin.?,        0.0534186165764990, TEST_TOLERANCE);
    try std.testing.expectApproxEqAbs(output.observables.items[0].norm.?,        1.0000000000003577, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic RTP on Time-Linear Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -16, 16 }}, .npoint = 256 },
        .initial_conditions = .{ .momentum = &.{0}, .position = &.{0}, .state = 1, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .time_linear = .{} },
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
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 256 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .tully_1 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.001,
            .stop_norm = 0.1,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.005,
            .e_max = 0.225,
            .e_step = 0.01,
        },
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
    try std.testing.expectApproxEqAbs(sum, 32.1358941260022700, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic RTP Cross Section on Tully's First Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 256 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = false },
        .potential = .{ .tully_1 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.001,
            .stop_norm = 0.1,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.005,
            .e_max = 0.225,
            .e_step = 0.01,
        },
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
    try std.testing.expectApproxEqAbs(sum, 31.7477856126643400, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic RTP Cross Section on Tully's Second Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 256 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .tully_2 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.02,
            .stop_norm = 0.1,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.04,
            .e_max = 0.24,
            .e_step = 0.01,
        },
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
    try std.testing.expectApproxEqAbs(sum, 22.4851554312020350, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic RTP Cross Section on Tully's Second Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 256 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = false },
        .potential = .{ .tully_2 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.02,
            .stop_norm = 0.1,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.04,
            .e_max = 0.24,
            .e_step = 0.01,
        },
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
    try std.testing.expectApproxEqAbs(sum, 22.3986753668110620, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic RTP Cross Section on Tully's Third Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 256 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .tully_3 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.002,
            .stop_norm = 0.1,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.005,
            .e_max = 0.225,
            .e_step = 0.01,
        },
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
    try std.testing.expectApproxEqAbs(sum, 27.7701571325001770, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic RTP Cross Section on Tully's Third Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 256 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = false },
        .potential = .{ .tully_3 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.002,
            .stop_norm = 0.1,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.005,
            .e_max = 0.225,
            .e_step = 0.01,
        },
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
    try std.testing.expectApproxEqAbs(sum, 28.3459476347574080, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic RTP Partial Waves on Tully's First Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 256 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .tully_1 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.001,
            .stop_norm = 0.1,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.005,
            .e_max = 0.225,
            .e_step = 0.01,
        },
        .partial_waves = .{ .j_min = 0, .j_max = 1 },
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
    try std.testing.expectApproxEqAbs(sum, 1.9569985319299963, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic RTP Partial Waves on Tully's First Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 256 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = false },
        .potential = .{ .tully_1 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.001,
            .stop_norm = 0.1,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.005,
            .e_max = 0.225,
            .e_step = 0.01,
        },
        .partial_waves = .{ .j_min = 0, .j_max = 1 },
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
    try std.testing.expectApproxEqAbs(sum, 1.9526946359649326, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic RTP Partial Waves on Tully's Second Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 256 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .tully_2 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.02,
            .stop_norm = 0.1,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.04,
            .e_max = 0.24,
            .e_step = 0.01,
        },
        .partial_waves = .{ .j_min = 0, .j_max = 1 },
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
    try std.testing.expectApproxEqAbs(sum, 1.6933143443872576, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic RTP Partial Waves on Tully's Second Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 256 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = false },
        .potential = .{ .tully_2 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.02,
            .stop_norm = 0.1,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.04,
            .e_max = 0.24,
            .e_step = 0.01,
        },
        .partial_waves = .{ .j_min = 0, .j_max = 1 },
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
    try std.testing.expectApproxEqAbs(sum, 1.6893078938289168, TEST_TOLERANCE);
    // zig fmt: on
}

test "Adiabatic RTP Partial Waves on Tully's Third Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 256 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = true },
        .potential = .{ .tully_3 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.002,
            .stop_norm = 0.1,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.005,
            .e_max = 0.225,
            .e_step = 0.01,
        },
        .partial_waves = .{ .j_min = 0, .j_max = 1 },
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
    try std.testing.expectApproxEqAbs(sum, 2.1027485069149610, TEST_TOLERANCE);
    // zig fmt: on
}

test "Diabatic RTP Partial Waves on Tully's Third Potential" {
    const opt = zinq.quantum_dynamics.Options{
        .grid = .{ .bounds = &.{.{ -24, 32 }}, .npoint = 256 },
        .initial_conditions = .{ .momentum = &.{15}, .position = &.{-10}, .state = 0, .gamma = &.{2}, .adiabatic = false },
        .potential = .{ .tully_3 = .{} },
        .absorbing_potential = .{
            .bounds = &.{.{ -18, 18 }},
            .exponent = 0.002,
            .stop_norm = 0.1,
        },
        .flux_analysis = .{
            .flux_bounds = &.{.{ -15, 15 }},
            .e_min = 0.005,
            .e_max = 0.225,
            .e_step = 0.01,
        },
        .partial_waves = .{ .j_min = 0, .j_max = 1 },
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
    try std.testing.expectApproxEqAbs(sum, 2.3511800557613420, TEST_TOLERANCE);
    // zig fmt: on
}
