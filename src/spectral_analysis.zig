//! Computes absorption or emission spectra from wavefunction autocorrelation functions via Fourier transform.

const std = @import("std");

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;

const Vector = @import("tensor.zig").Vector;

const fft = @import("fourier_transform.zig").fft;

/// Calculates the spectrum by Fourier transforming a windowed and zero-padded autocorrelation function.
pub fn calcSpectrum(comptime T: type, acf: Vector(Complex(T)), dt: T, padding: usize, thresh: T, gpa: Allocator) !struct { Vector(Complex(T)), Vector(T) } {
    const decay = -@log(thresh);

    const n = try std.math.ceilPowerOfTwo(usize, @max((acf.length() + padding - 1) * 2, 2));

    var acfp = try Vector(Complex(T)).initZero(n / 2 + 1, gpa);
    errdefer acfp.deinit(gpa);

    var spec = try Vector(T).init(n, gpa);
    errdefer spec.deinit(gpa);

    for (0..acf.length()) |i| {
        const x = @as(T, @floatFromInt(i)) / @as(T, @floatFromInt(acf.length()));

        const window = std.math.exp(-decay * x * x);

        acfp.ptr(i).* = Complex(T).init(acf.at(i).re * window, acf.at(i).im * window);
    }

    var full = try gpa.alloc(Complex(T), n);
    defer gpa.free(full);

    for (0..n / 2 + 1) |k| {
        full[k] = acfp.at(k);
    }

    for (n / 2 + 1..n) |k| {
        const orig = full[n - k];

        full[k] = Complex(T).init(orig.re, -orig.im);
    }

    fft(T, full, full, 1);

    const n_flt = @as(T, @floatFromInt(n));

    for (0..n) |i| {
        spec.ptr(i).* = full[i].re * n_flt;
    }

    for (0..n / 2) |i| {
        const j, const temp = .{ i + n / 2, spec.at(i) };

        spec.ptr(i).*, spec.ptr(j).* = .{ spec.at(j) * dt, temp * dt };
    }

    return .{ acfp, spec };
}
