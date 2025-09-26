//! File with all the necessary functions to compute the derivative of the energy with respect to nuclear coordinates.

const std = @import("std");

const classical_particle = @import("classical_particle.zig");
const device_write = @import("device_write.zig");
const real_matrix = @import("real_matrix.zig");

const ClassicalParticle = classical_particle.ClassicalParticle;
const RealMatrix = real_matrix.RealMatrix;

const print = device_write.print;

/// Function that returns a gradient given the system and energy function.
pub fn nuclearGradient(comptime T: type, opt: anytype, system: ClassicalParticle(T), efunc: anytype, method: []const u8, enable_printing: bool, allocator: std.mem.Allocator) !RealMatrix(T) {
    var grad = try RealMatrix(T).init(system.position.len / 3, 3, allocator);

    if (enable_printing) try print("\n{s} NUMERICAL GRADIENT:\n{s:7} {s:20} {s:4}\n", .{method, "INDEX", "GRADIENT ELEMENT", "TIME"});

    for (0..system.position.len) |i| {

        var timer = try std.time.Timer.start();

        if (enable_printing) try print("{d:3}/{d:3}", .{i, system.position.len});

        var sysp1 = try system.clone(); defer sysp1.deinit(); sysp1.position.ptr(i).* += opt.gradient.?.numeric.step;
        var sysm1 = try system.clone(); defer sysm1.deinit(); sysm1.position.ptr(i).* -= opt.gradient.?.numeric.step;

        var outp1 = try efunc(T, opt, sysp1, false, allocator); const Ep1 = outp1.energy; outp1.deinit();
        var outm1 = try efunc(T, opt, sysm1, false, allocator); const Em1 = outm1.energy; outm1.deinit();

        grad.ptr(i / 3, i % 3).* = (Ep1 - Em1) / (2 * opt.gradient.?.numeric.step);

        if (enable_printing) try print(" {d:20.14} {D}\n", .{grad.at(i / 3, i % 3), timer.read()});
    }

    return grad;
}

/// Returns a Hessian given the system and energy function.
pub fn nuclearHessian(comptime T: type, opt: anytype, system: ClassicalParticle(T), efunc: anytype, method: []const u8, enable_printing: bool, allocator: std.mem.Allocator) !RealMatrix(T) {
    var hess = try RealMatrix(T).init(system.position.len, system.position.len, allocator); var k: usize = 0;

    if (enable_printing) try print("\n{s} NUMERICAL HESSIAN:\n{s:7} {s:20} {s:4}\n", .{method, "INDEX", "HESSIAN ELEMENT", "TIME"});

    for (0..hess.rows) |i| for (i..hess.cols) |j| {

        var timer = try std.time.Timer.start();

        if (enable_printing) try print("{d:3}/{d:3}", .{k + 1, hess.rows * (hess.rows - 1) / 2 + hess.rows});

        var sysp2   = try system.clone(); defer   sysp2.deinit(); sysp2.position.ptr(i).*   += opt.hessian.?.numeric.step; sysp2.position.ptr(j).*   += opt.hessian.?.numeric.step;
        var sysp1m1 = try system.clone(); defer sysp1m1.deinit(); sysp1m1.position.ptr(i).* += opt.hessian.?.numeric.step; sysp1m1.position.ptr(j).* -= opt.hessian.?.numeric.step;
        var sysm1p1 = try system.clone(); defer sysm1p1.deinit(); sysm1p1.position.ptr(i).* -= opt.hessian.?.numeric.step; sysm1p1.position.ptr(j).* += opt.hessian.?.numeric.step;
        var sysm2   = try system.clone(); defer   sysm2.deinit(); sysm2.position.ptr(i).*   -= opt.hessian.?.numeric.step; sysm2.position.ptr(j).*   -= opt.hessian.?.numeric.step;

        var outp2   = try efunc(T, opt, sysp2,   false, allocator); const Ep2   = outp2.energy;     outp2.deinit();
        var outp1m1 = try efunc(T, opt, sysp1m1, false, allocator); const Ep1m1 = outp1m1.energy; outp1m1.deinit();
        var outm1p1 = try efunc(T, opt, sysm1p1, false, allocator); const Em1p1 = outm1p1.energy; outm1p1.deinit();
        var outm2   = try efunc(T, opt, sysm2,   false, allocator); const Em2   = outm2.energy;     outm2.deinit();

        hess.ptr(i, j).* = (Ep2 - Ep1m1 - Em1p1 + Em2) / (4 * opt.hessian.?.numeric.step * opt.hessian.?.numeric.step); hess.ptr(j, i).* = hess.at(i, j); k += 1;

        if (enable_printing) try print(" {d:20.14} {D}\n", .{hess.at(i, j), timer.read()});
    };

    return hess;
}
