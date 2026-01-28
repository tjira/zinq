//! Orbit fractal generator module.

const std = @import("std");

const fractal_generator = @import("fractal_generator.zig");
const image = @import("image.zig");
const rgb = @import("rgb.zig");

const Complex = std.math.Complex;
const Image = image.Image;
const RGB = rgb.RGB;

/// Fractal types.
pub fn Fractal(comptime T: type) type {
    return union(enum) {
        mandelbrot: Mandelbrot(T),
    };
}

/// Mandelbrot fractal structure.
pub fn Mandelbrot(comptime T: type) type {
    return struct {
        z:  ?Complex(T) = Complex(T).init(0, 0),
        zp: ?Complex(T) = Complex(T).init(0, 0),
        
        /// Mandelbrot iteration function.
        pub fn iterate(_: @This(), p: Complex(T), z: Complex(T), _: Complex(T)) Complex(T) {
            return z.mul(z).add(p);
        }
    };
}

/// Algorithm types.
pub fn Algorithm(comptime T: type) type {
    return union(enum) {
        escape: Escape(T),
    };
}

/// Escape time algorithm structure.
pub fn Escape(comptime T: type) type {
    return struct {
        bailout: T = 4, maxiter: u32 = 100,
    };
}

/// Coloring types.
pub fn Coloring(comptime T: type) type {
    return union(enum) {
        gradient: Gradient(T),
        solid: Solid(T),
    };
}

/// Gradient coloring structure.
pub fn Gradient(comptime T: type) type {
    return struct {
        start: RGB = .{.r =   0, .g =   0, .b =   0},
        end:   RGB = .{.r = 255, .g = 255, .b = 255},

        /// Color getter from a fractal normalized iteration value.
        pub fn get(self: @This(), value: T) RGB {
            const r = @as(u8, @intFromFloat(@as(T, @floatFromInt(self.start.r)) + value * (@as(T, @floatFromInt(self.end.r)) - @as(T, @floatFromInt(self.start.r)))));
            const g = @as(u8, @intFromFloat(@as(T, @floatFromInt(self.start.g)) + value * (@as(T, @floatFromInt(self.end.g)) - @as(T, @floatFromInt(self.start.g)))));
            const b = @as(u8, @intFromFloat(@as(T, @floatFromInt(self.start.b)) + value * (@as(T, @floatFromInt(self.end.b)) - @as(T, @floatFromInt(self.start.b)))));

            return RGB{.r = r, .g = g, .b = b};
        }
    };
}

/// Solid coloring structure.
pub fn Solid(comptime T: type) type {
    return struct {
        rgb: RGB = .{.r = 255, .g = 255, .b = 255},

        /// Color getter from a fractal normalized iteration value.
        pub fn get(self: @This(), _: T) RGB {
            return self.rgb;
        }
    };
}

/// Orbit fractal generator.
pub fn OrbitFractalGenerator(comptime T: type) type {
    return struct {
        fractal: Fractal(T), algorithm: Algorithm(T), coloring: Coloring(T), center: Complex(T), zoom: T, smooth: bool,

        /// Initialize the orbit fractal generator.
        pub fn init(opt: std.meta.TagPayload(fractal_generator.Options(T).Category, .orbit)) @This() {
            const fractal = switch (opt.fractal) {
                .mandelbrot => Fractal(T){.mandelbrot = Mandelbrot(T){}}
            };

            const algorithm = switch (opt.algorithm) {
                .escape => |alg| Algorithm(T){.escape = alg}
            };

            const coloring = switch (opt.coloring) {
                .solid => |col| Coloring(T){.solid = col},
                .gradient => |col| Coloring(T){.gradient = col}
            };

            const center = Complex(T).init(opt.center[0], opt.center[1]);

            return @This(){
                .fractal = fractal, .algorithm = algorithm, .coloring = coloring, .center = center, .zoom = opt.zoom, .smooth = opt.smooth
            };
        }

        /// Orbit fractal painting function.
        pub fn paint(self: @This(), canvas: *Image) void {
            switch (self.algorithm) {
                .escape => |opt| switch (self.fractal) {inline else => |fractal| {switch (self.coloring) {inline else => |coloring| {
                    self.escape(canvas, fractal, coloring, opt.bailout, opt.maxiter);
                }}}}
            }
        }

        /// Escape time algorithm.
        pub fn escape(self: @This(), canvas: *Image, fractal: anytype, coloring: anytype, bailout: T, maxiter: u32) void {
            const hf: T = @floatFromInt(canvas.height); const wf: T = @floatFromInt(canvas.width);

            for (0..canvas.height) |i| for (0..canvas.width) |j| {

                const im = -self.center.im + (3.0 * (@as(T, @floatFromInt(i)) + 0.5) - 1.5 * hf) / self.zoom / hf;
                const re =  self.center.re + (3.0 * (@as(T, @floatFromInt(j)) + 0.5) - 1.5 * wf) / self.zoom / hf;

                const p = Complex(T).init(re, im); var iter: u32 = 0;

                var z  = if (fractal.z ) |z | z  else p;
                var zp = if (fractal.zp) |zp| zp else p;

                while (z.squaredMagnitude() <= bailout * bailout and iter < maxiter) : (iter += 1) {
                    const zt = z; z = fractal.iterate(p, z, zp); zp = zt;
                }

                if (iter < maxiter) {

                    var value: T = @floatFromInt(iter);

                    if (self.smooth) value -= std.math.log2(std.math.log2(z.magnitude()));

                    canvas.ptr(i, j).* = coloring.get(value / @as(T, @floatFromInt(maxiter)));
                }
            };
        }
    };
}
