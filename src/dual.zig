//! Implements dual numbers to enable automatic differentiation for scientific and numerical computations.

const std = @import("std");

/// Returns a dual number type for computing value and first derivative of a scalar function.
pub fn ScalarDual(comptime T: type) type {
    return struct {
        val: T,
        der: T,

        /// Initializes a dual number with a real value and its corresponding first-order derivative.
        pub fn init(val: T, der: T) @This() {
            return .{ .val = val, .der = der };
        }

        /// Computes the absolute value and its derivative using the sign function derivative step.
        pub fn abs(self: @This()) @This() {
            return .{ .val = @abs(self.val), .der = if (self.val >= 0) self.der else -self.der };
        }

        /// Computes the sum of two dual numbers using the linearity of differentiation: d(u+v) = du + dv.
        pub fn add(self: @This(), other: @This()) @This() {
            return .{ .val = self.val + other.val, .der = self.der + other.der };
        }

        /// Adds a constant scalar to a dual number, leaving the derivative term unchanged.
        pub fn adds(self: @This(), scalar: T) @This() {
            return .{ .val = self.val + scalar, .der = self.der };
        }

        /// Divides two dual numbers using the quotient rule: d(u/v) = (du*v - u*dv) / v^2.
        pub fn div(self: @This(), other: @This()) @This() {
            const der = (self.der * other.val - self.val * other.der) / (other.val * other.val);

            return .{ .val = self.val / other.val, .der = der };
        }

        /// Divides a dual number by a constant scalar, scaling both value and derivative.
        pub fn divs(self: @This(), scalar: T) @This() {
            return .{ .val = self.val / scalar, .der = self.der / scalar };
        }

        /// Computes the exponential of a dual number using the chain rule: d(e^u) = e^u * du.
        pub fn exp(self: @This()) @This() {
            const val = std.math.exp(self.val);

            return .{ .val = val, .der = val * self.der };
        }

        /// Multiplies two dual numbers using the product rule: d(uv) = u*dv + v*du.
        pub fn mul(self: @This(), other: @This()) @This() {
            return .{ .val = self.val * other.val, .der = self.val * other.der + other.val * self.der };
        }

        /// Multiplies a dual number by a constant scalar, scaling both value and derivative.
        pub fn muls(self: @This(), scalar: T) @This() {
            return .{ .val = self.val * scalar, .der = self.der * scalar };
        }

        /// Computes the difference of two dual numbers: d(u-v) = du - dv.
        pub fn sub(self: @This(), other: @This()) @This() {
            return .{ .val = self.val - other.val, .der = self.der - other.der };
        }

        /// Subtracts a constant scalar from a dual number, leaving the derivative term unchanged.
        pub fn subs(self: @This(), scalar: T) @This() {
            return .{ .val = self.val - scalar, .der = self.der };
        }
    };
}
