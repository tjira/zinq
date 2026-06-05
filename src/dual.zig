pub fn Dual(comptime T: type, N: usize) type {
    return struct {
        // zig fmt: off
        val: T, grad: [N]T,
        // zig fmt: on

        pub fn init(val: T, i: usize) @This() {
            var grad: [N]T = undefined;

            inline for (0..N) |i| {
                grad[i] = 0;
            }

            grad[i] = 1;

            return .{ .val = val, .grad = grad };
        }

        pub fn add(self: @This(), other: @This()) @This() {
            var grad: [N]T = undefined;

            inline for (0..N) |i| {
                grad[i] = self.grad[i] + other.grad[i];
            }

            return .{ .val = self.val + other.val, .grad = grad };
        }

        pub fn mul(self: @This(), other: @This()) @This() {
            var grad: [N]T = undefined;

            inline for (0..N) |i| {
                grad[i] = self.val * other.grad[i] + other.val * self.grad[i];
            }

            return .{ .val = self.val * other.val, .grad = grad };
        }

        pub fn muls(self: @This(), scalar: T) @This() {
            var grad: [N]T = undefined;

            inline for (0..N) |i| {
                grad[i] = self.grad[i] * scalar;
            }

            return .{ .val = self.val * scalar, .grad = grad };
        }
    };
}
