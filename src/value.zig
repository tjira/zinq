//! Polymorphic value type wrapping floating-point, complex, and dual numbers for algebraic operations.

const std = @import("std");

const Complex = std.math.Complex;

/// Returns a polymorphic value wrapper providing unified algebraic interfaces for multiple numeric types.
pub fn Value(comptime T: type) type {
    const U = primType(T);

    return struct {
        val: T,

        /// Initializes a Value wrapper with the given numeric data.
        pub fn init(val: T) @This() {
            return .{ .val = val };
        }

        /// Computes the absolute value (modulus) of the wrapped value.
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

        /// Computes the sum of two values.
        pub fn add(self: @This(), other: @This()) @This() {
            if (comptime isFloat(T)) {
                return init(self.val + other.val);
            }

            if (comptime isDual(T) or isComplex(T)) {
                return init(self.val.add(other.val));
            }

            @compileError("TYPE '" ++ @typeName(T) ++ "' UNSUPPORTED");
        }

        /// Adds a primitive scalar to the value.
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

        /// Divides the value by another value.
        pub fn div(self: @This(), other: @This()) @This() {
            if (comptime isFloat(T)) {
                return init(self.val / other.val);
            }

            if (comptime isDual(T) or isComplex(T)) {
                return init(self.val.div(other.val));
            }

            @compileError("TYPE '" ++ @typeName(T) ++ "' UNSUPPORTED");
        }

        /// Divides the value by a primitive scalar.
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

        /// Computes the exponential of the value.
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

        /// Promotes a primitive float to a Value wrapper of type T.
        pub fn fromFloat(val: U) @This() {
            if (comptime isFloat(T)) {
                return init(val);
            }

            if (comptime isDual(T) or isComplex(T)) {
                return init(T.init(val, 0));
            }

            @compileError("TYPE '" ++ @typeName(T) ++ "' UNSUPPORTED");
        }

        /// Computes the product of two values.
        pub fn mul(self: @This(), other: @This()) @This() {
            if (comptime isFloat(T)) {
                return init(self.val * other.val);
            }

            if (comptime isDual(T) or isComplex(T)) {
                return init(self.val.mul(other.val));
            }

            @compileError("TYPE '" ++ @typeName(T) ++ "' UNSUPPORTED");
        }

        /// Multiplies the value by a primitive scalar.
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

        /// Computes the additive inverse (negation) of the value.
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

        /// Computes the sign or phase (normalized value) of the wrapped number.
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

        /// Computes the difference between two values.
        pub fn sub(self: @This(), other: @This()) @This() {
            if (comptime isFloat(T)) {
                return init(self.val - other.val);
            }

            if (comptime isDual(T) or isComplex(T)) {
                return init(self.val.sub(other.val));
            }

            @compileError("TYPE '" ++ @typeName(T) ++ "' UNSUPPORTED");
        }

        /// Subtracts a primitive scalar from the value.
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
    };
}

/// Returns the underlying primitive float type (e.g. f64) of a scalar, complex, or dual number.
pub fn primType(comptime T: type) type {
    return if (isFloat(T)) T else @typeInfo(T).@"struct".fields[0].type;
}

/// Returns true if the type T is a complex number representation.
pub fn isComplex(comptime T: type) bool {
    if (@typeInfo(T) == .@"struct") {
        return @hasField(T, "re") and @hasField(T, "im");
    }

    return false;
}

/// Returns true if the type T is a dual number representation.
pub fn isDual(comptime T: type) bool {
    if (@typeInfo(T) == .@"struct") {
        return @hasField(T, "val") and @hasField(T, "der");
    }

    return false;
}

/// Returns true if the type T is a native floating-point type.
fn isFloat(comptime T: type) bool {
    return @typeInfo(T) == .float;
}
