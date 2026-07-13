//! Wrapper around FFTW for performing fast Fourier transforms on multi-dimensional complex arrays.

const fftw = @cImport(@cInclude("fftw3.h"));

const primType = @import("value.zig").primType;

/// Returns an FFT plan type that encapsulates execution details for discrete Fourier transforms.
pub fn FftPlan(comptime T: type) type {
    return struct {
        plan: fftw.fftw_plan,

        sign: i32,

        /// Initializes a discrete Fourier transform plan for the specified shape, direction, and planning mode.
        pub fn init(arr: []T, shape: []i32, sign: i32, mode: u32) !@This() {
            if (comptime primType(T) != f64) @compileError("FFTW NOW ONLY SUPPORTS F64 NUMBERS");

            const ptr = @as([*c]fftw.fftw_complex, @ptrCast(arr.ptr));

            const plan = fftw.fftw_plan_dft(@intCast(shape.len), shape.ptr, ptr, ptr, sign, mode);

            return .{ .plan = plan orelse return error.PlanCreationFailed, .sign = sign };
        }

        /// Destroys the FFTW plan, releasing all internal resources and plans.
        pub fn deinit(self: @This()) void {
            fftw.fftw_destroy_plan(self.plan);
        }

        /// Creates a duplicate of the existing FFTW plan with identical transform properties.
        pub fn clone(self: @This()) !@This() {
            const plan = fftw.fftw_copy_plan(self.plan);

            return .{ .plan = plan orelse return error.PlanDuplicationFailed, .sign = self.sign };
        }

        /// Executes the DFT in-place and applies normalization (1/N) if computing an inverse transform.
        pub fn execute(self: @This(), arr: []T) void {
            const ptr = @as([*c]fftw.fftw_complex, @ptrCast(arr.ptr));

            fftw.fftw_execute_dft(self.plan, ptr, ptr);

            if (self.sign == 1) for (0..arr.len) |i| {
                arr[i].re /= @as(f64, @floatFromInt(arr.len));
                arr[i].im /= @as(f64, @floatFromInt(arr.len));
            };
        }
    };
}
