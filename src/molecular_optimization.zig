const std = @import("std");

const libint = @import("cimport.zig").libint;

const Allocator = std.mem.Allocator;

const Matrix = @import("tensor.zig").Matrix;
const MolecularSystem = @import("molecular_system.zig").MolecularSystem;

const printf = @import("read_write.zig").printf;

pub fn steepestDescent(comptime T: type, io: std.Io, runFn: anytype, opt: anytype, sys: *MolecularSystem(T), Pg: ?Matrix(T), log: bool, gpa: Allocator) anyerror!Matrix(T) {
    var run_opt, const is_hf = .{ opt, @hasField(@TypeOf(opt), "hartree_fock") };

    run_opt.hessian, run_opt.optimize, run_opt.gradient = .{ null, null, opt.optimize.?.gradient };

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

    if (opt.optimize.?.iterations > 0 and log) {
        const fmt = "\nSTEEPEST DESCENT OPTIMIZATION\n{s:4} {s:20} {s:9} {s:9} {s:9}\n";

        try printf(io, fmt, .{ "ITER", "TOTAL ENERGY", "RMS(G)", "MAX(G)", "TIME" });
    }

    var guess_density: ?Matrix(T) = if (Pg) |guess| try guess.clone(gpa) else null;
    errdefer if (guess_density) |*p| p.deinit(gpa);

    for (0..opt.optimize.?.iterations) |i| {
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

            try printf(io, fmt, .{ i, res.energy[0], rms_grad, max_grad, time_str });
        }

        const P = try (if (comptime is_hf) res.hartree_fock.P else res.P).clone(gpa);

        if (guess_density) |*p| p.deinit(gpa);

        guess_density = P;

        if (rms_grad < opt.optimize.?.threshold) {
            break;
        }

        for (0..grad.nrow()) |j| for (0..grad.ncol()) |k| {
            sys.coors[3 * j + k] -= opt.optimize.?.step * grad.at(j, k);
        };

        libint.libint_update_coords(sys.ptr, sys.coors.ptr);
    } else return error.OptimizationDidNotConverge;

    return guess_density.?;
}
