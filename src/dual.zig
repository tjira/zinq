const std = @import("std");

pub fn ScalarDual(comptime T: type) type {
    return struct {
        val: T,
        der: T,

        pub fn init(val: T, der: T) @This() {
            return .{ .val = val, .der = der };
        }

        pub fn add(self: @This(), other: @This()) @This() {
            return .{ .val = self.val + other.val, .der = self.der + other.der };
        }

        pub fn adds(self: @This(), scalar: T) @This() {
            return .{ .val = self.val + scalar, .der = self.der };
        }

        pub fn sub(self: @This(), other: @This()) @This() {
            return .{ .val = self.val - other.val, .der = self.der - other.der };
        }

        pub fn subs(self: @This(), scalar: T) @This() {
            return .{ .val = self.val - scalar, .der = self.der };
        }

        pub fn mul(self: @This(), other: @This()) @This() {
            return .{ .val = self.val * other.val, .der = self.val * other.der + other.val * self.der };
        }

        pub fn muls(self: @This(), scalar: T) @This() {
            return .{ .val = self.val * scalar, .der = self.der * scalar };
        }

        pub fn div(self: @This(), other: @This()) @This() {
            const der = (self.der * other.val - self.val * other.der) / (other.val * other.val);

            return .{ .val = self.val / other.val, .der = der };
        }

        pub fn divs(self: @This(), scalar: T) @This() {
            return .{ .val = self.val / scalar, .der = self.der / scalar };
        }

        pub fn exp(self: @This()) @This() {
            const val = std.math.exp(self.val);

            return .{ .val = val, .der = val * self.der };
        }

        pub fn abs(self: @This()) @This() {
            return .{ .val = @abs(self.val), .der = if (self.val >= 0) self.der else -self.der };
        }

        pub fn sqrt(self: @This()) @This() {
            const val = @sqrt(self.val);

            return .{ .val = val, .der = if (val == 0) 0 else self.der / (2 * val) };
        }

        pub fn cos(self: @This()) @This() {
            return .{ .val = std.math.cos(self.val), .der = -self.der * std.math.sin(self.val) };
        }

        pub fn sin(self: @This()) @This() {
            return .{ .val = std.math.sin(self.val), .der = self.der * std.math.cos(self.val) };
        }

        pub fn clamp(self: @This(), min_val: T, max_val: T) @This() {
            if (self.val < min_val) {
                return .{ .val = min_val, .der = 0 };
            }

            if (self.val > max_val) {
                return .{ .val = max_val, .der = 0 };
            }

            return self;
        }
    };
}
