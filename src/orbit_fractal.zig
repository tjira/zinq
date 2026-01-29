//! Orbit fractal generator module.

const std = @import("std");

const fractal_generator = @import("fractal_generator.zig");

const Complex = std.math.Complex;

/// Orbit fractal types.
pub fn OrbitFractal(comptime T: type) type {
    return union(enum) {
        buffalo: Buffalo(T),
        burningship: BurningShip(T),
        julia: Julia(T),
        mandelbrot: Mandelbrot(T),
        manowar: Manowar(T),
        phoenix: Phoenix(T),
    };
}

/// Buffalo fractal structure.
pub fn Buffalo(comptime T: type) type {
    return struct {
        z: ?Complex(T) = Complex(T).init(0, 0), zp: ?Complex(T) = Complex(T).init(0, 0),

        /// Initialize function that returns the initial z and zp values.
        pub fn init(_: @This(), _: Complex(T)) struct{Complex(T), Complex(T)} {
            return .{Complex(T).init(0, 0), Complex(T).init(0, 0)};
        }
        
        /// Burningship iteration function.
        pub fn iterate(_: @This(), p: Complex(T), z: Complex(T), _: Complex(T)) Complex(T) {
            const zabs = Complex(T).init(@abs(z.re), @abs(z.im));

            return zabs.mul(zabs).sub(zabs).add(p);
        }
    };
}

/// Burning ship fractal structure.
pub fn BurningShip(comptime T: type) type {
    return struct {
        z: ?Complex(T) = Complex(T).init(0, 0), zp: ?Complex(T) = Complex(T).init(0, 0),

        /// Initialize function that returns the initial z and zp values.
        pub fn init(_: @This(), _: Complex(T)) struct{Complex(T), Complex(T)} {
            return .{Complex(T).init(0, 0), Complex(T).init(0, 0)};
        }
        
        /// Burningship iteration function.
        pub fn iterate(_: @This(), p: Complex(T), z: Complex(T), _: Complex(T)) Complex(T) {
            const zabs = Complex(T).init(@abs(z.re), @abs(z.im));

            return zabs.mul(zabs).add(p);
        }
    };
}

/// Julia fractal structure.
pub fn Julia(comptime T: type) type {
    return struct {
        z: ?Complex(T) = null, zp: ?Complex(T) = Complex(T).init(0, 0),

        c: Complex(T),

        /// Init function that returns the initial z and zp values.
        pub fn init(_: @This(), p: Complex(T)) struct{Complex(T), Complex(T)} {
            return .{p, Complex(T).init(0, 0)};
        }

        /// Julia iteration function.
        pub fn iterate(self: @This(), _: Complex(T), z: Complex(T), _: Complex(T)) Complex(T) {
            return z.mul(z).add(self.c);
        }
    };
}

/// Mandelbrot fractal structure.
pub fn Mandelbrot(comptime T: type) type {
    return struct {
        z: ?Complex(T) = Complex(T).init(0, 0), zp: ?Complex(T) = Complex(T).init(0, 0),

        /// Init function that returns the initial z and zp values.
        pub fn init(_: @This(), _: Complex(T)) struct{Complex(T), Complex(T)} {
            return .{Complex(T).init(0, 0), Complex(T).init(0, 0)};
        }
        
        /// Mandelbrot iteration function.
        pub fn iterate(_: @This(), p: Complex(T), z: Complex(T), _: Complex(T)) Complex(T) {
            return z.mul(z).add(p);
        }
    };
}

/// Manowar fractal structure.
pub fn Manowar(comptime T: type) type {
    return struct {
        z: ?Complex(T) = null, zp: ?Complex(T) = null,

        /// Initialize function that returns the initial z and zp values.
        pub fn init(_: @This(), p: Complex(T)) struct{Complex(T), Complex(T)} {
            return .{p, p};
        }
        
        /// Manowar iteration function.
        pub fn iterate(_: @This(), p: Complex(T), z: Complex(T), zp: Complex(T)) Complex(T) {
            return z.mul(z).add(zp).add(p);
        }
    };
}

/// Phoenix fractal structure.
pub fn Phoenix(comptime T: type) type {
    return struct {
        z: ?Complex(T) = null, zp: ?Complex(T) = Complex(T).init(0, 0),

        c: Complex(T),

        /// Init function that returns the initial z and zp values.
        pub fn init(self: @This(), p: Complex(T)) struct{Complex(T), Complex(T)} {
            return .{p, self.c};
        }
        
        /// Phoenix iteration function.
        pub fn iterate(_: @This(), _: Complex(T), z: Complex(T), zp: Complex(T)) Complex(T) {
            return z.mul(z).sub(zp.mul(Complex(T).init(0.5, 0))).add(Complex(T).init(0.5667, 0));
        }
    };
}
