//! Fourier transform module.

const std = @import("std");

const array_functions = @import("array_functions.zig");
const complex_vector = @import("complex_vector.zig");
const error_handling = @import("error_handling.zig");
const global_variables = @import("global_variables.zig");
const math_functions = @import("math_functions.zig");
const real_vector = @import("real_vector.zig");
const strided_complex_vector = @import("strided_complex_vector.zig");

const Complex = std.math.Complex;
const ComplexVector = complex_vector.ComplexVector;
const StridedComplexVector = strided_complex_vector.StridedComplexVector;
const RealVector = real_vector.RealVector;

const prod = array_functions.prod;
const revk = math_functions.revk;
const throw = error_handling.throw;

const TEST_TOLERANCE = global_variables.TEST_TOLERANCE;

/// Fast Fourier transform for a one-dimensional array. The factor argument is the value in the exponent of the Fourier transform. Factor -1 corresponds to the forward Fourier transform, while factor 1 corresponds to the inverse Fourier transform.
pub fn cfft1(comptime T: type, vector: *StridedComplexVector(T), factor: i32) !void {
    const n = vector.len; const bit_count: u6 = @intCast(std.math.log2(n));

    if (std.math.pow(usize, 2, @intCast(bit_count)) != n) {
        return throw(void, "FOURIER TRANSFORM ONLY IMPLEMENTED FOR LENGTHS THAT ARE POWERS OF 2", .{});
    }

    for (0..n) |i| {

        const j = revk(i, bit_count);

        if (i < j) {
            std.mem.swap(std.math.Complex(T), vector.ptr(j), vector.ptr(i));
        }
    }

    for (0..bit_count) |stage| {

        const block_size = std.math.pow(usize, 2, stage + 1);
        const angle = 2 * std.math.pi * @as(T, @floatFromInt(factor)) / @as(T, @floatFromInt(block_size));
        const twiddle_base = std.math.complex.exp(Complex(T).init(0, angle));

        for (0..n / block_size) |j| {

            var omega = Complex(T).init(1, 0);

            for (0..block_size / 2) |k| {

                const t = omega.mul(vector.at(j * block_size + k + block_size / 2));

                const u = vector.at(j * block_size + k);

                vector.ptr(j * block_size + k                 ).* = u.add(t);
                vector.ptr(j * block_size + k + block_size / 2).* = u.sub(t);

                omega = omega.mul(twiddle_base);
            }
        }
    }

    if (factor > 0) for (0..n) |i| {
        vector.ptr(i).* = vector.at(i).div(std.math.Complex(T).init(@as(T, @floatFromInt(n)), 0));
    };
}

/// Fast Fourier transform for an n-dimensional array. The factor argument is the value in the exponent of the Fourier transform. Factor -1 corresponds to the forward Fourier transform, while factor 1 corresponds to the inverse Fourier transform.
pub fn cfftn(comptime T: type, vector: *ComplexVector(T), shape: []const usize, factor: i32) !void {
    const N = prod(usize, shape);

    if (std.math.pow(usize, shape[0], shape.len) != N) {
        return throw(void, "FOURIER TRANSFORM ONLY IMPLEMENTED FOR CUBIC ARRAYS", .{});
    }

    for (0..shape.len) |i| {

        const stride = prod(usize, shape[0..i]);

        for (0..N / shape[i]) |j| {

            var offset: usize = 0; var index: usize = 0;

            for (0..shape.len) |k| if (k != i) {

                const block = prod(usize, shape[k + 1..]) / (if (k < i) shape[i] else 1);

                offset += (j / block % shape[k]) * prod(usize, shape[0..k]); index += 1;
            };
            
            var axis_slice = StridedComplexVector(T){.data = vector.data, .len = shape[i], .stride = stride, .zero = offset};

            try cfft1(T, &axis_slice, factor);
        }
    }
}

test "cfft1" {
    var v = try ComplexVector(f64).init(8, std.testing.allocator); defer v.deinit();

    v.ptr(0).* = Complex(f64).init(0, 0);
    v.ptr(1).* = Complex(f64).init(1, 0);
    v.ptr(2).* = Complex(f64).init(2, 0);
    v.ptr(3).* = Complex(f64).init(3, 0);
    v.ptr(4).* = Complex(f64).init(4, 0);
    v.ptr(5).* = Complex(f64).init(5, 0);
    v.ptr(6).* = Complex(f64).init(6, 0);
    v.ptr(7).* = Complex(f64).init(7, 0);

    const v_fft = try ComplexVector(f64).init(8, std.testing.allocator); defer v_fft.deinit();

    v_fft.ptr(0).* = Complex(f64).init(28.00000000000000,  0.00000000000000);
    v_fft.ptr(1).* = Complex(f64).init(-4.00000000000000,  9.65685424949238);
    v_fft.ptr(2).* = Complex(f64).init(-4.00000000000000,  4.00000000000000);
    v_fft.ptr(3).* = Complex(f64).init(-4.00000000000000,  1.65685424949238);
    v_fft.ptr(4).* = Complex(f64).init(-4.00000000000000,  0.00000000000000);
    v_fft.ptr(5).* = Complex(f64).init(-4.00000000000000, -1.65685424949238);
    v_fft.ptr(6).* = Complex(f64).init(-4.00000000000000, -4.00000000000000);
    v_fft.ptr(7).* = Complex(f64).init(-4.00000000000000, -9.65685424949238);

    var v_strided = v.asStrided();

    try cfft1(f64, &v_strided, -1);

    try std.testing.expect(v.eq(v_fft, TEST_TOLERANCE));
}
