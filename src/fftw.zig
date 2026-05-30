const std = @import("std");

const fftw = @cImport(@cInclude("fftw3.h"));

pub fn fftn(comptime T: type, arr: []T, shape: []i32, sign: i32) !void {
    const ptr = @as([*c]fftw.fftw_complex, @ptrCast(arr.ptr));

    const plan = fftw.fftw_plan_dft(@intCast(shape.len), shape.ptr, ptr, ptr, sign, fftw.FFTW_ESTIMATE);
    defer fftw.fftw_destroy_plan(plan);

    fftw.fftw_execute(plan);

    if (sign == 1) for (arr) |*e| {
        e.*.re /= @as(f64, @floatFromInt(arr.len));
        e.*.im /= @as(f64, @floatFromInt(arr.len));
    };
}
