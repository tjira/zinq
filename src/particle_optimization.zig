//! File with all the necessary functions to compute the derivative of the energy with respect to nuclear coordinates.

const std = @import("std");

const energy_derivative = @import("energy_derivative.zig");
const error_handling = @import("error_handling.zig");
const classical_particle = @import("classical_particle.zig");
const device_write = @import("device_write.zig");
const real_matrix = @import("real_matrix.zig");

const ClassicalParticle = classical_particle.ClassicalParticle;
const RealMatrix = real_matrix.RealMatrix;

const nuclearGradient = energy_derivative.nuclearGradient;
const throw = error_handling.throw;
const print = device_write.print;

/// Optimizes the positions of a set of classical particles using the steepest descent method and a gradient.
pub fn particleSteepestDescent(comptime T: type, opt: anytype, system: ClassicalParticle(T), efunc: anytype, method: []const u8, enable_printing: bool, allocator: std.mem.Allocator) !ClassicalParticle(T) {
    if (opt.gradient == null) return throw(ClassicalParticle(T), "NO GRADIENT METHOD SPECIFIED FOR OPTIMIZATION", .{});

    var optimized_system = try system.clone();

    if (enable_printing) try print("\n{s} GEOMETRY OPTIMIZATION:\n{s:4} {s:20} {s:4}\n", .{method, "ITER", "GRADIENT NORM", "TIME"});

    for (0..opt.optimize.?.maxiter + 1) |i| {

        if (i == opt.optimize.?.maxiter) return throw(ClassicalParticle(T), "MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED", .{});

        var timer = try std.time.Timer.start();

        if (enable_printing) try print("{d:4}", .{i + 1});

        var G = try nuclearGradient(T, opt, optimized_system, efunc, method, false, allocator);

        if (enable_printing) try print(" {d:20.14} {D}\n", .{G.asVector().norm(), timer.read()});

        if (G.asVector().norm() < opt.optimize.?.threshold) break;

        G.muls(opt.optimize.?.step);

        try optimized_system.position.sub(G.asVector());
    }

    return optimized_system;
}
