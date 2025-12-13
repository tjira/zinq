//! File with all the necessary functions to compute the derivative of the energy with respect to nuclear coordinates.

const std = @import("std");

const classical_particle = @import("classical_particle.zig");
const device_write = @import("device_write.zig");
const error_context = @import("error_context.zig");
const real_matrix = @import("real_matrix.zig");

const ClassicalParticle = classical_particle.ClassicalParticle;
const ErrorContext = error_context.ErrorContext;
const RealMatrix = real_matrix.RealMatrix;

const print = device_write.print;

/// Calculate the partial derivative of the energy with respect to nuclear coordinates using finite differences and assign it to the provided result.
pub fn firstPartialDerivative(comptime T: type, result: *T, opt: anytype, system: ClassicalParticle(T), efunc: anytype, i: usize, enable_printing: bool, allocator: std.mem.Allocator) !void {
    var timer = try std.time.Timer.start();

    var sysp1 = try system.clone(allocator); defer sysp1.deinit(allocator); sysp1.position.ptr(i).* += opt.gradient.?.numeric.step;
    var sysm1 = try system.clone(allocator); defer sysm1.deinit(allocator); sysm1.position.ptr(i).* -= opt.gradient.?.numeric.step;

    var outp1 = try efunc(T, opt, sysp1, false, allocator); const Ep1 = outp1.energy; outp1.deinit(allocator);
    var outm1 = try efunc(T, opt, sysm1, false, allocator); const Em1 = outm1.energy; outm1.deinit(allocator);

    result.* = (Ep1 - Em1) / (2 * opt.gradient.?.numeric.step);

    if (enable_printing) try print("{d:3}/{d:3} {d:20.14} {D}\n", .{i + 1, system.position.len, result.*, timer.read()});
}

/// Function that returns a gradient given the system and energy function.
pub fn nuclearGradient(comptime T: type, opt: anytype, system: ClassicalParticle(T), efunc: anytype, method: []const u8, enable_printing: bool, allocator: std.mem.Allocator) !RealMatrix(T) {
    var grad = try RealMatrix(T).init(system.position.len / 3, 3, allocator);

    if (enable_printing) try print("\n{s} NUMERICAL GRADIENT:\n{s:7} {s:20} {s:4}\n", .{method, "INDEX", "GRADIENT ELEMENT", "TIME"});

    var pool: std.Thread.Pool = undefined; try pool.init(.{.n_jobs = opt.gradient.?.numeric.nthread, .allocator = allocator}); var error_ctx = ErrorContext{};

    const parallel_function = struct {
        pub fn call(result: *T, params: anytype, err_ctx: *ErrorContext) void {
            if (err_ctx.err == null) firstPartialDerivative(T, result, params[0], params[1], params[2], params[3], params[4], params[5]) catch |err| {
                err_ctx.capture(err); return;
            };
        }
    }.call;

    for (0..system.position.len) |i| {

        const params = .{opt, system, efunc, i, enable_printing, allocator}; const result = grad.ptr(i / 3, i % 3);

        if (opt.gradient.?.numeric.nthread == 1) parallel_function(result, params, &error_ctx) else try pool.spawn(parallel_function, .{result, params, &error_ctx});
    }

    pool.deinit(); if (error_ctx.err) |err| return err;

    return grad;
}

/// Returns a Hessian given the system and energy function.
pub fn nuclearHessian(comptime T: type, opt: anytype, system: ClassicalParticle(T), efunc: anytype, method: []const u8, enable_printing: bool, allocator: std.mem.Allocator) !RealMatrix(T) {
    var hess = try RealMatrix(T).init(system.position.len, system.position.len, allocator);

    if (enable_printing) try print("\n{s} NUMERICAL HESSIAN:\n{s:7} {s:20} {s:4}\n", .{method, "INDEX", "HESSIAN ELEMENT", "TIME"});

    var pool: std.Thread.Pool = undefined; try pool.init(.{.n_jobs = opt.hessian.?.numeric.nthread, .allocator = allocator}); var error_ctx = ErrorContext{};

    const parallel_function = struct {
        pub fn call(result: *T, params: anytype, err_ctx: *ErrorContext) void {
            if (err_ctx.err == null) secondPartialDerivative(T, result, params[0], params[1], params[2], params[3], params[4], params[5], params[6]) catch |err| {
                err_ctx.capture(err); return;
            };
        }
    }.call;

    for (0..hess.rows) |i| for (i..hess.cols) |j| {

        const params = .{opt, system, efunc, i, j, enable_printing, allocator}; const result = hess.ptr(i, j);

        if (opt.hessian.?.numeric.nthread == 1) parallel_function(result, params, &error_ctx) else try pool.spawn(parallel_function, .{result, params, &error_ctx});
    };

    pool.deinit(); if (error_ctx.err) |err| return err;

    for (0..hess.rows) |i| for (i + 1..hess.cols) |j| {
        hess.ptr(j, i).* = hess.at(i, j);
    };

    return hess;
}

/// Calculate the mixed second derivative of the energy with respect to nuclear coordinates using finite differences and assign it to the provided result.
pub fn secondPartialDerivative(comptime T: type, result: *T, opt: anytype, system: ClassicalParticle(T), efunc: anytype, i: usize, j: usize, enable_printing: bool, allocator: std.mem.Allocator) !void {
    var timer = try std.time.Timer.start();

    const index_of_derivative = i * system.position.len - (i * (i - 1)) / 2 + (j - i) + 1;
    const total_derivatives = system.position.len * (system.position.len - 1) / 2 + system.position.len;

    var sysp2   = try system.clone(allocator); defer   sysp2.deinit(allocator); sysp2.position.ptr(i).*   += opt.hessian.?.numeric.step; sysp2.position.ptr(j).*   += opt.hessian.?.numeric.step;
    var sysp1m1 = try system.clone(allocator); defer sysp1m1.deinit(allocator); sysp1m1.position.ptr(i).* += opt.hessian.?.numeric.step; sysp1m1.position.ptr(j).* -= opt.hessian.?.numeric.step;
    var sysm1p1 = try system.clone(allocator); defer sysm1p1.deinit(allocator); sysm1p1.position.ptr(i).* -= opt.hessian.?.numeric.step; sysm1p1.position.ptr(j).* += opt.hessian.?.numeric.step;
    var sysm2   = try system.clone(allocator); defer   sysm2.deinit(allocator); sysm2.position.ptr(i).*   -= opt.hessian.?.numeric.step; sysm2.position.ptr(j).*   -= opt.hessian.?.numeric.step;

    var outp2   = try efunc(T, opt, sysp2,   false, allocator); const Ep2   = outp2.energy;     outp2.deinit(allocator);
    var outp1m1 = try efunc(T, opt, sysp1m1, false, allocator); const Ep1m1 = outp1m1.energy; outp1m1.deinit(allocator);
    var outm1p1 = try efunc(T, opt, sysm1p1, false, allocator); const Em1p1 = outm1p1.energy; outm1p1.deinit(allocator);
    var outm2   = try efunc(T, opt, sysm2,   false, allocator); const Em2   = outm2.energy;     outm2.deinit(allocator);

    result.* = (Ep2 - Ep1m1 - Em1p1 + Em2) / (4 * opt.hessian.?.numeric.step * opt.hessian.?.numeric.step);

    if (enable_printing) try print("{d:3}/{d:3} {d:20.14} {D}\n", .{index_of_derivative, total_derivatives, result.*, timer.read()});
}
