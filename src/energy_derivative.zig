//! File with all the necessary functions to compute the derivative of the energy with respect to nuclear coordinates.

const std = @import("std");

const classical_particle = @import("classical_particle.zig");
const device_write = @import("device_write.zig");
const real_matrix = @import("real_matrix.zig");

const ClassicalParticle = classical_particle.ClassicalParticle;
const RealMatrix = real_matrix.RealMatrix;

const print = device_write.print;

/// Calculate the partial derivative of the energy with respect to nuclear coordinates using finite differences and assign it to the provided result.
pub fn firstPartialDerivative(comptime T: type, io: std.Io, result: *T, opt: anytype, system: ClassicalParticle(T), efunc: anytype, i: usize, enable_printing: bool, allocator: std.mem.Allocator) !void {
    var timer = std.Io.Timestamp.now(io, .real);

    var sysp1 = try system.clone(allocator);
    defer sysp1.deinit(allocator);
    sysp1.position.ptr(i).* += opt.gradient.?.numeric.step;

    var sysm1 = try system.clone(allocator);
    defer sysm1.deinit(allocator);
    sysm1.position.ptr(i).* -= opt.gradient.?.numeric.step;

    var outp1 = try efunc(T, io, opt, sysp1, false, allocator);
    const Ep1 = outp1.energy;
    outp1.deinit(allocator);

    var outm1 = try efunc(T, io, opt, sysm1, false, allocator);
    const Em1 = outm1.energy;
    outm1.deinit(allocator);

    result.* = (Ep1 - Em1) / (2 * opt.gradient.?.numeric.step);

    if (enable_printing) try print(io, "{d:3}/{d:3} {d:20.14} {f}\n", .{ i + 1, system.position.len, result.*, timer.untilNow(io, .real) });
}

/// Function that returns a gradient given the system and energy function.
pub fn nuclearGradient(comptime T: type, io: std.Io, opt: anytype, system: ClassicalParticle(T), efunc: anytype, method: []const u8, enable_printing: bool, allocator: std.mem.Allocator) !RealMatrix(T) {
    var grad = try RealMatrix(T).init(system.position.len / 3, 3, allocator);

    if (enable_printing) try print(io, "\n{s} NUMERICAL GRADIENT:\n{s:7} {s:20} {s:4}\n", .{ method, "INDEX", "GRADIENT ELEMENT", "TIME" });

    for (0..system.position.len) |i| {
        const result = grad.ptr(i / 3, i % 3);

        try firstPartialDerivative(T, io, result, opt, system, efunc, i, enable_printing, allocator);
    }

    return grad;
}

/// Returns a Hessian given the system and energy function.
pub fn nuclearHessian(comptime T: type, io: std.Io, opt: anytype, system: ClassicalParticle(T), efunc: anytype, method: []const u8, enable_printing: bool, allocator: std.mem.Allocator) !RealMatrix(T) {
    var hess = try RealMatrix(T).init(system.position.len, system.position.len, allocator);

    if (enable_printing) try print(io, "\n{s} NUMERICAL HESSIAN:\n{s:7} {s:20} {s:4}\n", .{ method, "INDEX", "HESSIAN ELEMENT", "TIME" });

    for (0..hess.rows) |i| for (i..hess.cols) |j| {
        const result = hess.ptr(i, j);

        try secondPartialDerivative(T, io, result, opt, system, efunc, i, j, enable_printing, allocator);
    };

    for (0..hess.rows) |i| for (i + 1..hess.cols) |j| {
        hess.ptr(j, i).* = hess.at(i, j);
    };

    return hess;
}

/// Calculate the mixed second derivative of the energy with respect to nuclear coordinates using finite differences and assign it to the provided result.
pub fn secondPartialDerivative(comptime T: type, io: std.Io, result: *T, opt: anytype, system: ClassicalParticle(T), efunc: anytype, i: usize, j: usize, enable_printing: bool, allocator: std.mem.Allocator) !void {
    var timer = std.Io.Timestamp.now(io, .real);

    const index_of_derivative = i * system.position.len - (i * (i - 1)) / 2 + (j - i) + 1;
    const total_derivatives = system.position.len * (system.position.len - 1) / 2 + system.position.len;

    var sysp2 = try system.clone(allocator);
    defer sysp2.deinit(allocator);
    sysp2.position.ptr(i).* += opt.hessian.?.numeric.step;
    sysp2.position.ptr(j).* += opt.hessian.?.numeric.step;

    var sysp1m1 = try system.clone(allocator);
    defer sysp1m1.deinit(allocator);
    sysp1m1.position.ptr(i).* += opt.hessian.?.numeric.step;
    sysp1m1.position.ptr(j).* -= opt.hessian.?.numeric.step;

    var sysm1p1 = try system.clone(allocator);
    defer sysm1p1.deinit(allocator);
    sysm1p1.position.ptr(i).* -= opt.hessian.?.numeric.step;
    sysm1p1.position.ptr(j).* += opt.hessian.?.numeric.step;

    var sysm2 = try system.clone(allocator);
    defer sysm2.deinit(allocator);
    sysm2.position.ptr(i).* -= opt.hessian.?.numeric.step;
    sysm2.position.ptr(j).* -= opt.hessian.?.numeric.step;

    var outp2 = try efunc(T, io, opt, sysp2, false, allocator);
    const Ep2 = outp2.energy;
    outp2.deinit(allocator);

    var outp1m1 = try efunc(T, io, opt, sysp1m1, false, allocator);
    const Ep1m1 = outp1m1.energy;
    outp1m1.deinit(allocator);

    var outm1p1 = try efunc(T, io, opt, sysm1p1, false, allocator);
    const Em1p1 = outm1p1.energy;
    outm1p1.deinit(allocator);

    var outm2 = try efunc(T, io, opt, sysm2, false, allocator);
    const Em2 = outm2.energy;
    outm2.deinit(allocator);

    result.* = (Ep2 - Ep1m1 - Em1p1 + Em2) / (4 * opt.hessian.?.numeric.step * opt.hessian.?.numeric.step);

    if (enable_printing) try print(io, "{d:3}/{d:3} {d:20.14} {f}\n", .{ index_of_derivative, total_derivatives, result.*, timer.untilNow(io, .real) });
}
