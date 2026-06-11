const std = @import("std");

const Complex = std.math.Complex;

const ScalarDual = @import("dual.zig").ScalarDual;

pub fn Value(comptime T: type) type {
    const U = primType(T);

    return struct {
        val: T,

        pub fn init(val: T) @This() {
            return .{ .val = val };
        }

        pub fn fromFloat(val: U) @This() {
            if (comptime isFloat(T)) {
                return init(val);
            }

            if (comptime isDual(T) or isComplex(T)) {
                return init(T.init(val, 0));
            }

            @compileError("TYPE '" ++ @typeName(T) ++ "' UNSUPPORTED");
        }

        pub fn add(self: @This(), other: @This()) @This() {
            if (comptime isFloat(T)) {
                return init(self.val + other.val);
            }

            if (comptime isDual(T) or isComplex(T)) {
                return init(self.val.add(other.val));
            }

            @compileError("TYPE '" ++ @typeName(T) ++ "' UNSUPPORTED");
        }

        pub fn adds(self: @This(), scalar: U) @This() {
            if (comptime isFloat(T)) {
                return init(self.val + scalar);
            }

            if (comptime isDual(T)) {
                return init(self.val.adds(scalar));
            }

            if (comptime isComplex(T)) {
                return init(self.val.add(T.init(scalar, 0)));
            }

            @compileError("TYPE '" ++ @typeName(T) ++ "' UNSUPPORTED");
        }

        pub fn sub(self: @This(), other: @This()) @This() {
            if (comptime isFloat(T)) {
                return init(self.val - other.val);
            }

            if (comptime isDual(T) or isComplex(T)) {
                return init(self.val.sub(other.val));
            }

            @compileError("TYPE '" ++ @typeName(T) ++ "' UNSUPPORTED");
        }

        pub fn subs(self: @This(), scalar: U) @This() {
            if (comptime isFloat(T)) {
                return init(self.val - scalar);
            }

            if (comptime isDual(T)) {
                return init(self.val.subs(scalar));
            }

            if (comptime isComplex(T)) {
                return init(self.val.sub(T.init(scalar, 0)));
            }

            @compileError("TYPE '" ++ @typeName(T) ++ "' UNSUPPORTED");
        }

        pub fn mul(self: @This(), other: @This()) @This() {
            if (comptime isFloat(T)) {
                return init(self.val * other.val);
            }

            if (comptime isDual(T) or isComplex(T)) {
                return init(self.val.mul(other.val));
            }

            @compileError("TYPE '" ++ @typeName(T) ++ "' UNSUPPORTED");
        }

        pub fn muls(self: @This(), scalar: U) @This() {
            if (comptime isFloat(T)) {
                return init(self.val * scalar);
            }

            if (comptime isDual(T)) {
                return init(self.val.muls(scalar));
            }

            if (comptime isComplex(T)) {
                return init(self.val.mul(T.init(scalar, 0)));
            }

            @compileError("TYPE '" ++ @typeName(T) ++ "' UNSUPPORTED");
        }

        pub fn div(self: @This(), other: @This()) @This() {
            if (comptime isFloat(T)) {
                return init(self.val / other.val);
            }

            if (comptime isDual(T) or isComplex(T)) {
                return init(self.val.div(other.val));
            }

            @compileError("TYPE '" ++ @typeName(T) ++ "' UNSUPPORTED");
        }

        pub fn divs(self: @This(), scalar: U) @This() {
            if (comptime isFloat(T)) {
                return init(self.val / scalar);
            }

            if (comptime isDual(T)) {
                return init(self.val.divs(scalar));
            }

            if (comptime isComplex(T)) {
                return init(self.val.div(T.init(scalar, 0)));
            }

            @compileError("TYPE '" ++ @typeName(T) ++ "' UNSUPPORTED");
        }

        pub fn neg(self: @This()) @This() {
            if (comptime isFloat(T)) {
                return init(-self.val);
            }

            if (comptime isDual(T)) {
                return init(T.init(-self.val.val, -self.val.der));
            }

            if (comptime isComplex(T)) {
                return init(T.init(-self.val.re, -self.val.im));
            }

            @compileError("TYPE '" ++ @typeName(T) ++ "' UNSUPPORTED");
        }

        pub fn abs(self: @This()) @This() {
            if (comptime isFloat(T)) {
                return init(@abs(self.val));
            }

            if (comptime isDual(T)) {
                return init(self.val.abs());
            }

            if (comptime isComplex(T)) {
                return init(T.init(std.math.complex.abs(self.val), 0));
            }

            @compileError("TYPE '" ++ @typeName(T) ++ "' UNSUPPORTED");
        }

        pub fn exp(self: @This()) @This() {
            if (comptime isFloat(T)) {
                return init(std.math.exp(self.val));
            }

            if (comptime isDual(T)) {
                return init(self.val.exp());
            }

            if (comptime isComplex(T)) {
                return init(std.math.complex.exp(self.val));
            }

            @compileError("TYPE '" ++ @typeName(T) ++ "' UNSUPPORTED");
        }

        pub fn sign(self: @This()) @This() {
            if (comptime isFloat(T)) {
                return init(std.math.sign(self.val));
            }

            if (comptime isDual(T)) {
                return init(T.init(std.math.sign(self.val.val), 0));
            }

            if (comptime isComplex(T)) {
                const a = std.math.complex.abs(self.val);

                if (a == 0) return init(T.init(0, 0));

                return init(self.val.div(T.init(a, 0)));
            }

            @compileError("TYPE '" ++ @typeName(T) ++ "' UNSUPPORTED");
        }
    };
}

pub fn isComplex(comptime T: type) bool {
    if (@typeInfo(T) == .@"struct") {
        return @hasField(T, "re") and @hasField(T, "im");
    }

    return false;
}

pub fn isDual(comptime T: type) bool {
    if (@typeInfo(T) == .@"struct") {
        return @hasField(T, "val") and @hasField(T, "der");
    }

    return false;
}

pub fn isFloat(comptime T: type) bool {
    return @typeInfo(T) == .float;
}

pub fn primType(comptime T: type) type {
    return if (isFloat(T)) T else @typeInfo(T).@"struct".fields[0].type;
}
