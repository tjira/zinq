//! Computes absorption or emission spectra from wavefunction autocorrelation functions via Fourier transform.

const std = @import("std");

const fftw = @import("cimport.zig").fftw;

const Allocator = std.mem.Allocator;
const Complex = std.math.Complex;

const Vector = @import("tensor.zig").Vector;

/// Calculates the spectrum by Fourier transforming a windowed and zero-padded autocorrelation function.
pub fn calcSpectrum(comptime T: type, acf: Vector(Complex(T)), dt: T, padding: usize, thresh: T, gpa: Allocator) !struct { Vector(Complex(T)), Vector(T) } {
    const flags, const decay = .{ fftw.FFTW_ESTIMATE | fftw.FFTW_PRESERVE_INPUT, -@log(thresh) };

    var acfp = try Vector(Complex(T)).initZero(acf.length() + padding, gpa);
    errdefer acfp.deinit(gpa);

    var spec = try Vector(T).init((acfp.length() - 1) * 2, gpa);
    errdefer spec.deinit(gpa);

    const inptr = @as([*c]fftw.fftw_complex, @ptrCast(acfp.data.ptr));

    const plan = fftw.fftw_plan_dft_c2r_1d(@intCast(spec.length()), inptr, spec.data.ptr, flags);
    defer fftw.fftw_destroy_plan(plan);

    for (0..acf.length()) |i| {
        const x = @as(T, @floatFromInt(i)) / @as(T, @floatFromInt(acf.length()));

        const window = std.math.exp(-decay * x * x);

        acfp.ptr(i).* = Complex(T).init(acf.at(i).re * window, acf.at(i).im * window);
    }

    fftw.fftw_execute(plan);

    for (0..spec.length() / 2) |i| {
        const j, const temp = .{ i + spec.length() / 2, spec.at(i) };

        spec.ptr(i).*, spec.ptr(j).* = .{ spec.at(j) * dt, temp * dt };
    }

    return .{ acfp, spec };
}
