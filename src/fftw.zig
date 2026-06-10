const std = @import("std");

const fftw = @cImport(@cInclude("fftw3.h"));

const primType = @import("value.zig").primType;

pub fn FftPlan(comptime T: type) type {
    return struct {
        plan: fftw.fftw_plan,

        sign: i32,

        pub fn init(arr: []T, shape: []i32, sign: i32, mode: u32) !@This() {
            if (comptime primType(T) != f64) @compileError("FFTW NOW ONLY SUPPORTS F64 NUMBERS");

            const ptr = @as([*c]fftw.fftw_complex, @ptrCast(arr.ptr));

            const plan = fftw.fftw_plan_dft(@intCast(shape.len), shape.ptr, ptr, ptr, sign, mode);

            return .{ .plan = plan orelse return error.PlanCreationFailed, .sign = sign };
        }

        pub fn deinit(self: @This()) void {
            fftw.fftw_destroy_plan(self.plan);
        }

        pub fn clone(self: @This()) !@This() {
            const plan = fftw.fftw_copy_plan(self.plan);

            return .{ .plan = plan orelse return error.PlanDuplicationFailed, .sign = self.sign };
        }

        pub fn execute(self: @This(), arr: []T) void {
            const ptr = @as([*c]fftw.fftw_complex, @ptrCast(arr.ptr));

            fftw.fftw_execute_dft(self.plan, ptr, ptr);

            if (self.sign == 1) for (arr) |*e| {
                e.*.re /= @as(f64, @floatFromInt(arr.len));
                e.*.im /= @as(f64, @floatFromInt(arr.len));
            };
        }
    };
}
