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
    };
}
