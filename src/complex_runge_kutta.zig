//! Runge-Kutta integrator with complex numbers.

const std = @import("std");

const complex_vector = @import("complex_vector.zig");

const Complex = std.math.complex.Complex;
const ComplexVector = complex_vector.ComplexVector;

pub fn ComplexRungeKutta(comptime T: type) type {
    return struct {
        k1: ComplexVector(T),
        k2: ComplexVector(T),
        k3: ComplexVector(T),
        k4: ComplexVector(T),
        y1: ComplexVector(T),
        y2: ComplexVector(T),
        y3: ComplexVector(T),

        /// Initialize the Runge-Kutta struct.
        pub fn init(size: usize, allocator: std.mem.Allocator) !@This() {
            return @This(){
                .k1 = try ComplexVector(T).initZero(size, allocator),
                .k2 = try ComplexVector(T).initZero(size, allocator),
                .k3 = try ComplexVector(T).initZero(size, allocator),
                .k4 = try ComplexVector(T).initZero(size, allocator),
                .y1 = try ComplexVector(T).initZero(size, allocator),
                .y2 = try ComplexVector(T).initZero(size, allocator),
                .y3 = try ComplexVector(T).initZero(size, allocator)
            };
        }

        /// Free the memory allocated for the Runge-Kutta struct.
        pub fn deinit(self: @This()) void {
            self.k1.deinit();
            self.k2.deinit();
            self.k3.deinit();
            self.k4.deinit();
            self.y1.deinit();
            self.y2.deinit();
            self.y3.deinit();
        }

        /// Propagate the wavefunction coefficients using the Runge-Kutta method.
        pub fn propagate(self: *@This(), amplitudes: *ComplexVector(T), function: anytype, parameters: anytype, time_step: T) void {
            self.k1.zero(); self.k2.zero(); self.k3.zero(); self.k4.zero();

            function(&self.k1, amplitudes.*, parameters);

            for (0..amplitudes.len) |j| {
                self.y1.ptr(j).* = amplitudes.at(j).add(self.k1.at(j).mul(Complex(T).init(time_step / 2, 0)));
            }

            function(&self.k2, self.y1, parameters);

            for (0..amplitudes.len) |j| {
                self.y2.ptr(j).* = amplitudes.at(j).add(self.k2.at(j).mul(Complex(T).init(time_step / 2, 0)));
            }

            function(&self.k3, self.y2, parameters);

            for (0..amplitudes.len) |j| {
                self.y3.ptr(j).* = amplitudes.at(j).add(self.k3.at(j).mul(Complex(T).init(time_step, 0)));
            }

            function(&self.k4, self.y3, parameters);

            for (0..amplitudes.len) |j| {

                const ksum = self.k1.at(j).add(self.k2.at(j).mul(Complex(T).init(2, 0))).add(self.k3.at(j).mul(Complex(T).init(2, 0))).add(self.k4.at(j));

                amplitudes.ptr(j).* = amplitudes.at(j).add(ksum.mul(Complex(T).init(time_step / 6, 0)));
            }
        }
    };
}
