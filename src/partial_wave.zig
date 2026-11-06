//! Computes total and partial quantum reactive scattering cross sections via angular momentum expansion.

const std = @import("std");

const Allocator = std.mem.Allocator;

const FluxAnalysis = @import("flux_analysis.zig").FluxAnalysis;
const Matrix = @import("tensor.zig").Matrix;
const Potential = @import("potential.zig").Potential;

const printf = @import("read_write.zig").printf;

/// Configuration options for quantum partial wave scattering calculations.
pub const PartialWaveOptions = struct {
    log_interval: u32 = 1,
    max: ?u32 = null,
    step: u32 = 1,
    threads: ?u32 = null,
    write: PartialWaveWrite = .{},
};

/// Output file destinations for persisting computed partial wave scattering observables.
pub const PartialWaveWrite = struct {
    cross_section: ?[]const u8 = null,
};

/// Thread-safe counters and lock for parallel partial wave distribution.
pub fn PartialWaveContext(comptime T: type) type {
    return struct {
        next_j: std.atomic.Value(u32),
        in_mutx: std.Io.Mutex = .init,

        completed: std.atomic.Value(u32) = .init(0),
        mutx: std.atomic.Value(bool) = .init(false),

        j_max: u32,

        /// Propagates partial wavepackets concurrently and accumulates state-resolved reaction cross sections.
        pub fn worker(self: *@This(), io: std.Io, opt: anytype, log: bool, gpa: Allocator, thread_sigma: ?*Matrix(T), comptime initFn: anytype, comptime solveFn: anytype) void {
            const pw = opt.partial_waves.?;

            var arena = std.heap.ArenaAllocator.init(gpa);
            defer arena.deinit();

            const tasks = self.j_max / pw.step + 1;

            while (true) {
                const j = self.next_j.fetchAdd(pw.step, .monotonic);

                if (j > self.j_max) break;

                const timer = std.Io.Timestamp.now(io, .real);

                _ = arena.reset(.retain_capacity);

                var opt_j = opt;

                opt_j.j_quantum_number = j;
                opt_j.partial_waves = null;

                const alloc = arena.allocator();

                inline for (@typeInfo(@TypeOf(opt_j.write)).@"struct".field_names) |name| {
                    if (@field(opt_j.write, name)) |p| {
                        @field(opt_j.write, name) = injectAngularFname(p, j, alloc) catch null;
                    }
                }

                if (opt_j.flux_analysis) |*fa| if (fa.write.cross_section) |p| {
                    fa.write.cross_section = injectAngularFname(p, j, alloc) catch null;
                };

                self.in_mutx.lockUncancelable(io);

                var sim = initFn(T, io, opt_j, alloc) catch |err| {
                    self.in_mutx.unlock(io);

                    std.log.err("FAILED TO INITIALIZE SIMULATION FOR J={d} ({s})", .{ j, @errorName(err) });

                    continue;
                };

                self.in_mutx.unlock(io);

                defer sim.deinit(alloc);

                var res = solveFn(T, io, .{ .opt = opt_j, .sim = &sim, .eigs = 0, .log = false }, alloc) catch |err| {
                    std.log.err("SIMULATION FAILED FOR J={d} ({s})", .{ j, @errorName(err) });

                    continue;
                };

                defer res.deinit(alloc);

                var p_peak: T, var p_integ: T = .{ 0, 0 };

                if (res.cross_section) |*s_j| {
                    const fa = opt.flux_analysis.?;

                    const flux_res = integratePartialFlux(T, s_j, fa, opt.initial_conditions.state);

                    p_peak, p_integ = .{ flux_res.peak, flux_res.integ };

                    if (thread_sigma) |sigma| {
                        accumulatePartialCrossSection(T, sigma, s_j, fa, j, pw.step, opt.mass[0], gpa);
                    }
                }

                const r_0, const m_0 = .{ opt.initial_conditions.position[0], opt.mass[0] };

                const j_flt = @as(T, @floatFromInt(j));

                const v_cent = if (r_0 != 0) j_flt * (j_flt + 1) / (2 * m_0 * r_0 * r_0) else 0;

                const done = self.completed.fetchAdd(1, .monotonic) + 1;

                if (log and (done == 1 or (done - 1) % pw.log_interval == 0 or done == tasks)) {
                    self.logWorkerProgress(io, done - 1, j, v_cent, p_peak, p_integ, timer.untilNow(io, .real));
                }
            }
        }

        /// Accumulates energy-resolved reaction probability into the thread cross section matrix.
        fn accumulatePartialCrossSection(comptime U: type, sigma: *Matrix(U), s_j: *const Matrix(U), fa: anytype, j: u32, step: u32, mass_0: f64, gpa: Allocator) void {
            if (sigma.data.len == 0) {
                sigma.* = Matrix(U).initZero(s_j.nrow(), s_j.ncol(), gpa) catch return;
            }

            const dj = @as(U, @floatFromInt(step));

            const factor = dj * std.math.pi * (2 * @as(U, @floatFromInt(j)) + 1) / (2 * @as(U, @floatCast(mass_0)));

            const asymp: U = @floatCast(fa.asymptotic_offset);

            for (0..s_j.nrow()) |ei| {
                const E = FluxAnalysis(U).getCollisionEnergy(fa, ei) - asymp;

                if (E <= 0) continue;

                for (0..s_j.ncol()) |col| {
                    sigma.ptr(ei, col).* += (factor / E) * s_j.at(ei, col);
                }
            }
        }

        /// Integrates reaction probability across energies and finds peak probability.
        fn integratePartialFlux(comptime U: type, s_j: *const Matrix(U), fa: anytype, init_state: u32) struct { peak: U, integ: U } {
            var p_peak: U, var p_integ: U, var p_prev: U = .{ 0, 0, 0 };

            var e_prev: U = FluxAnalysis(U).getCollisionEnergy(fa, 0);

            for (0..s_j.nrow()) |ei| {
                const e_curr = FluxAnalysis(U).getCollisionEnergy(fa, ei);

                var p_row: U = 0;

                for (0..s_j.ncol()) |col| {
                    if (s_j.ncol() > 1 and col == init_state) continue;

                    p_row += s_j.at(ei, col);
                }

                if (p_row > p_peak) p_peak = p_row;

                if (ei > 0) {
                    p_integ += 0.5 * (p_prev + p_row) * (e_curr - e_prev);
                }

                p_prev, e_prev = .{ p_row, e_curr };
            }

            return .{ .peak = p_peak, .integ = p_integ };
        }

        /// Formats and logs partial wave worker progress to stdout under mutual exclusion.
        fn logWorkerProgress(self: *@This(), io: std.Io, idx: u32, j: u32, v_cent: T, p_peak: T, p_integ: T, elapsed: anytype) void {
            while (self.mutx.swap(true, .acquire)) {
                std.Thread.yield() catch {};
            }

            const fmt = "{d:06} {d:7} {d:12.6} {e:12.3} {e:12.3} {f}\n";

            printf(io, fmt, .{ idx, j, v_cent, p_peak, p_integ, elapsed }) catch {};

            self.mutx.store(false, .release);
        }
    };
}

/// Evaluates the maximum partial wave from the peak of the bare interaction potential barrier.
pub fn estimateMaxJ(comptime T: type, io: std.Io, opt: anytype, gpa: Allocator) !u32 {
    const fa, const npoint = .{ opt.flux_analysis.?, opt.grid.npoint };

    const dr = (opt.grid.bounds[0][1] - opt.grid.bounds[0][0]) / @as(f64, @floatFromInt(npoint - 1));

    const r_0, const r_flux = .{ opt.initial_conditions.position[0], fa.flux_bounds[0][1] };

    var pot = try Potential(T).init(io, opt.potential, gpa);
    defer pot.deinit(gpa);

    const v_buf = try gpa.alloc(T, pot.nstate() * pot.nstate());
    defer gpa.free(v_buf);

    const r_min, const r_max = .{ @min(r_0, r_flux), @max(r_0, r_flux) };

    const s_idx = opt.initial_conditions.state * (pot.nstate() + 1);

    var v_peak: f64, var r_peak: f64 = .{ -std.math.inf(f64), @abs(r_flux) };

    for (0..npoint) |i| {
        const r = opt.grid.bounds[0][0] + @as(f64, @floatFromInt(i)) * dr;

        if (r < r_min or r > r_max) continue;

        const r_t: T = @floatCast(r);

        try pot.eval(T, v_buf, &.{r_t}, 0);

        const v: f64 = @floatCast(v_buf[s_idx]);

        if (v > v_peak) v_peak, r_peak = .{ v, @max(@abs(r), dr) };
    }

    const k_eff = std.math.sqrt(2 * opt.mass[0] * @max(1.1 * fa.e_max - v_peak, 0));

    return @max(@as(u32, @intFromFloat(@ceil(k_eff * r_peak))) + 4, 1);
}

/// Injects the J partial wave index into output filenames.
pub fn injectAngularFname(path: []const u8, j: u32, allocator: Allocator) ![]const u8 {
    if (std.mem.indexOf(u8, path, "{J}")) |idx| {
        return try std.fmt.allocPrint(allocator, "{s}{d}{s}", .{ path[0..idx], j, path[idx + 3 ..] });
    }

    const last_dot = std.mem.lastIndexOfScalar(u8, path, '.');

    if (last_dot) |dot_idx| {
        return try std.fmt.allocPrint(allocator, "{s}_J{d}{s}", .{ path[0..dot_idx], j, path[dot_idx..] });
    }

    return try std.fmt.allocPrint(allocator, "{s}_J{d}", .{ path, j });
}
