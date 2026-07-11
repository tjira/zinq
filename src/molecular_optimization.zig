const std = @import("std");

const libint = @import("cimport.zig").libint;

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const Vector = @import("tensor.zig").Vector;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;

const printf = @import("read_write.zig").printf;

pub fn bfgs(comptime T: type, io: std.Io, runFn: anytype, opt: anytype, sys: *MolecularSystem(T), Pg: ?Matrix(T), log: bool, gpa: Allocator) anyerror!Matrix(T) {
    var run_opt, const is_hf = .{ opt, @hasField(@TypeOf(opt), "hartree_fock") };

    const bfgs_opt = opt.optimize.?.bfgs;

    run_opt.hessian, run_opt.optimize, run_opt.gradient = .{ null, null, bfgs_opt.gradient };

    if (comptime @hasField(@TypeOf(run_opt), "mulliken")) {
        run_opt.mulliken = false;
    }

    if (comptime @hasField(@TypeOf(run_opt), "response")) {
        run_opt.response = null;
    }

    if (comptime is_hf) {
        run_opt.hartree_fock.response = null;

        run_opt.hartree_fock.hessian, run_opt.hartree_fock.mulliken = .{ null, false };
    }

    if (bfgs_opt.iterations > 0 and log) {
        const fmt = "\nBFGS OPTIMIZATION\n{s:4} {s:20} {s:9} {s:9} {s:9}\n";

        try printf(io, fmt, .{ "ITER", "TOTAL ENERGY", "RMS(G)", "MAX(G)", "TIME" });
    }

    var guess_density: ?Matrix(T) = if (Pg) |guess| try guess.clone(gpa) else null;
    errdefer if (guess_density) |*p| p.deinit(gpa);

    var H = try Matrix(T).initZero(3 * sys.atoms.len, 3 * sys.atoms.len, gpa);
    defer H.deinit(gpa);

    for (0..3 * sys.atoms.len) |j| {
        H.ptr(j, j).* = 1.0;
    }

    var x_old = try Matrix(T).init(sys.atoms.len, 3, gpa);
    defer x_old.deinit(gpa);

    var g_old = try Matrix(T).init(sys.atoms.len, 3, gpa);
    defer g_old.deinit(gpa);

    for (0..bfgs_opt.iterations) |i| {
        var timer = std.Io.Timestamp.now(io, .real);

        var res = try runFn(T, io, run_opt, sys, guess_density, false, gpa);
        defer res.deinit(gpa);

        if (res.grad.len == 0) return error.GradientNotAvailable;

        const grad = res.grad[0];

        const max_grad = grad.max();
        const rms_grad = grad.rms();

        const elapsed = timer.untilNow(io, .real);

        if (log) {
            const fmt = "{d:4} {d:20.14} {e:9.3} {e:9.3} {s:9}\n";

            const time_str = try std.fmt.allocPrint(gpa, "{f}", .{elapsed});
            defer gpa.free(time_str);

            try printf(io, fmt, .{ i + 1, res.energy[0], rms_grad, max_grad, time_str });
        }

        const P = try (if (comptime is_hf) res.hartree_fock.P else res.P).clone(gpa);

        if (guess_density) |*p| p.deinit(gpa);

        guess_density = P;

        if (rms_grad < bfgs_opt.threshold) {
            break;
        }

        var x_new = try Matrix(T).init(sys.atoms.len, 3, gpa);
        defer x_new.deinit(gpa);

        for (0..x_new.nrow()) |j| for (0..x_new.ncol()) |k| {
            x_new.ptr(j, k).* = sys.coors[j * x_new.ncol() + k];
        };

        var g_new = try grad.clone(gpa);
        defer g_new.deinit(gpa);

        if (i > 0) {
            var s = try Vector(T).init(3 * sys.atoms.len, gpa);
            defer s.deinit(gpa);

            var y = try Vector(T).init(3 * sys.atoms.len, gpa);
            defer y.deinit(gpa);

            for (0..x_new.nrow()) |j| for (0..x_new.ncol()) |k| {
                const idx = j * x_new.ncol() + k;

                s.ptr(idx).* = x_new.at(j, k) - x_old.at(j, k);
                y.ptr(idx).* = g_new.at(j, k) - g_old.at(j, k);
            };

            var y_dot_s: T = 0;

            for (0..s.length()) |j| {
                y_dot_s += y.at(j) * s.at(j);
            }

            if (y_dot_s > 1e-12) {
                var u = try Vector(T).init(s.length(), gpa);
                defer u.deinit(gpa);

                for (0..u.length()) |j| {
                    var sum: T = 0;

                    for (0..y.length()) |k| {
                        sum += H.at(j, k) * y.at(k);
                    }

                    u.ptr(j).* = sum;
                }

                var y_dot_u: T = 0;

                for (0..y.length()) |j| {
                    y_dot_u += y.at(j) * u.at(j);
                }

                const rho, const c = .{ 1 / y_dot_s, (y_dot_u / y_dot_s + 1) / y_dot_s };

                for (0..H.nrow()) |j| for (0..H.ncol()) |k| {
                    const term1 = -rho * (s.at(j) * u.at(k) + u.at(j) * s.at(k));

                    H.ptr(j, k).* += term1 + c * s.at(j) * s.at(k);
                };
            }
        }

        for (0..x_new.nrow()) |j| for (0..x_new.ncol()) |k| {
            x_old.ptr(j, k).* = x_new.at(j, k);
            g_old.ptr(j, k).* = g_new.at(j, k);
        };

        var p = try Vector(T).init(3 * sys.atoms.len, gpa);
        defer p.deinit(gpa);

        for (0..p.length()) |j| {
            var sum: T = 0;

            for (0..g_new.nrow()) |row| for (0..g_new.ncol()) |col| {
                sum += H.at(j, row * g_new.ncol() + col) * g_new.at(row, col);
            };

            p.ptr(j).* = -sum;
        }

        for (0..p.length()) |j| {
            sys.coors[j] += bfgs_opt.step * p.at(j);
        }

        libint.libint_update_coords(sys.ptr, sys.coors.ptr);
    } else return error.OptimizationDidNotConverge;

    return guess_density.?;
}

pub fn steepestDescent(comptime T: type, io: std.Io, runFn: anytype, opt: anytype, sys: *MolecularSystem(T), Pg: ?Matrix(T), log: bool, gpa: Allocator) anyerror!Matrix(T) {
    var run_opt, const is_hf = .{ opt, @hasField(@TypeOf(opt), "hartree_fock") };

    const sd_opt = opt.optimize.?.steepest_descent;

    run_opt.hessian, run_opt.optimize, run_opt.gradient = .{ null, null, sd_opt.gradient };

    if (comptime @hasField(@TypeOf(run_opt), "mulliken")) {
        run_opt.mulliken = false;
    }

    if (comptime @hasField(@TypeOf(run_opt), "response")) {
        run_opt.response = null;
    }

    if (comptime is_hf) {
        run_opt.hartree_fock.response = null;

        run_opt.hartree_fock.hessian, run_opt.hartree_fock.mulliken = .{ null, false };
    }

    if (sd_opt.iterations > 0 and log) {
        const fmt = "\nSTEEPEST DESCENT OPTIMIZATION\n{s:4} {s:20} {s:9} {s:9} {s:9}\n";

        try printf(io, fmt, .{ "ITER", "TOTAL ENERGY", "RMS(G)", "MAX(G)", "TIME" });
    }

    var guess_density: ?Matrix(T) = if (Pg) |guess| try guess.clone(gpa) else null;
    errdefer if (guess_density) |*p| p.deinit(gpa);

    for (0..sd_opt.iterations) |i| {
        var timer = std.Io.Timestamp.now(io, .real);

        var res = try runFn(T, io, run_opt, sys, guess_density, false, gpa);
        defer res.deinit(gpa);

        if (res.grad.len == 0) return error.GradientNotAvailable;

        const grad = res.grad[0];

        const max_grad = grad.max();
        const rms_grad = grad.rms();

        const elapsed = timer.untilNow(io, .real);

        if (log) {
            const fmt = "{d:4} {d:20.14} {e:9.3} {e:9.3} {s:9}\n";

            const time_str = try std.fmt.allocPrint(gpa, "{f}", .{elapsed});
            defer gpa.free(time_str);

            try printf(io, fmt, .{ i + 1, res.energy[0], rms_grad, max_grad, time_str });
        }

        const P = try (if (comptime is_hf) res.hartree_fock.P else res.P).clone(gpa);

        if (guess_density) |*p| p.deinit(gpa);

        guess_density = P;

        if (rms_grad < sd_opt.threshold) {
            break;
        }

        for (0..grad.nrow()) |j| for (0..grad.ncol()) |k| {
            sys.coors[3 * j + k] -= sd_opt.step * grad.at(j, k);
        };

        libint.libint_update_coords(sys.ptr, sys.coors.ptr);
    } else return error.OptimizationDidNotConverge;

    return guess_density.?;
}
