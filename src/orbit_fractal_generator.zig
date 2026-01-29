//! Orbit fractal generator module.

const std = @import("std");

const rgb = @import("rgb.zig");
const image = @import("image.zig");
const orbit_fractal = @import("orbit_fractal.zig");
const fractal_generator = @import("fractal_generator.zig");

const Buffalo = orbit_fractal.Buffalo;
const BurningShip = orbit_fractal.BurningShip;
const Complex = std.math.Complex;
const Image = image.Image;
const Julia = orbit_fractal.Julia;
const Mandelbrot = orbit_fractal.Mandelbrot;
const Manowar = orbit_fractal.Manowar;
const OrbitFractal = orbit_fractal.OrbitFractal;
const Phoenix = orbit_fractal.Phoenix;
const RGB = rgb.RGB;

/// Algorithm type from the fractal generator options.
pub fn Algorithm(comptime T: type) type {
    return std.meta.TagPayload(fractal_generator.Options(T).Category, .orbit).Algorithm;
}

/// Coloring type from the fractal generator options.
pub fn Coloring(comptime T: type) type {
    return std.meta.TagPayload(fractal_generator.Options(T).Category, .orbit).Coloring;
}

/// Orbit fractal generator.
pub fn OrbitFractalGenerator(comptime T: type) type {
    return struct {
        fractal: OrbitFractal(T), algorithm: Algorithm(T), coloring: Coloring(T), center: Complex(T), zoom: T,

        /// Initialize the orbit fractal generator.
        pub fn init(opt: anytype) @This() {
            const fractal = switch (opt.fractal) {
                .buffalo => OrbitFractal(T){.buffalo = Buffalo(T){}},
                .burningship => OrbitFractal(T){.burningship = BurningShip(T){}},
                .julia => OrbitFractal(T){.julia = Julia(T){.c = opt.fractal.julia.c}},
                .mandelbrot => OrbitFractal(T){.mandelbrot = Mandelbrot(T){}},
                .manowar => OrbitFractal(T){.manowar = Manowar(T){}},
                .phoenix => OrbitFractal(T){.phoenix = Phoenix(T){.c = opt.fractal.phoenix.c}},
            };

            var generator = @This(){
                .fractal = fractal, .algorithm = opt.algorithm, .coloring = opt.coloring, .center = opt.center, .zoom = opt.zoom
            };

            if (opt.coloring == .gradient and opt.coloring.gradient.seed != null) {

                var seed: usize = @as(usize, opt.coloring.gradient.seed.?);

                if (seed == 0) seed = @intCast(std.time.milliTimestamp());

                var split_mix = std.Random.SplitMix64.init(seed); var rng = std.Random.DefaultPrng.init(split_mix.next());

                var random = rng.random();

                generator.coloring.gradient.start = RGB{.r = random.int(u8), .g = random.int(u8), .b = random.int(u8)};
                generator.coloring.gradient.end   = RGB{.r = random.int(u8), .g = random.int(u8), .b = random.int(u8)};
            }

            if (opt.coloring == .periodic and opt.coloring.periodic.seed != null) {

                var seed: usize = @as(usize, opt.coloring.periodic.seed.?);

                if (seed == 0) seed = @intCast(std.time.milliTimestamp());

                var split_mix = std.Random.SplitMix64.init(seed); var rng = std.Random.DefaultPrng.init(split_mix.next());

                var random = rng.random();

                for (0..3) |i| {
                    generator.coloring.periodic.amplitude[i] += 10 * std.math.pi * random.float(T);
                    generator.coloring.periodic.frequency[i] +=  2 * std.math.pi * random.float(T);
                }
            }

            if (opt.coloring == .solid and opt.coloring.solid.seed != null) {

                var seed: usize = @as(usize, opt.coloring.solid.seed.?);

                if (seed == 0) seed = @intCast(std.time.milliTimestamp());

                var split_mix = std.Random.SplitMix64.init(seed); var rng = std.Random.DefaultPrng.init(split_mix.next());

                var random = rng.random();

                generator.coloring.solid.rgb = RGB{.r = random.int(u8), .g = random.int(u8), .b = random.int(u8)};
            }

            return generator;
        }

        /// Orbit fractal painting function.
        pub fn paint(self: @This(), canvas: *Image) void {
            const coloring = switch (self.coloring) {
                .gradient => &@This().gradient, .periodic => &@This().periodic, .solid => &@This().solid
            };

            switch (self.algorithm) {
                .escape => |opt| switch (self.fractal) {inline else => |fractal| {
                    self.escape(canvas, fractal, coloring, opt.bailout, opt.maxiter, opt.smooth);
                }},
                .orbitrap => |opt| switch (self.fractal) {inline else => |fractal| {
                    self.orbitrap(canvas, fractal, coloring, opt.bailout, opt.maxiter, opt.fill);
                }},
            }
        }

        /// Escape time algorithm.
        pub fn escape(self: @This(), canvas: *Image, fractal: anytype, coloring: anytype, bailout: T, maxiter: u32, smooth: bool) void {
            const hf: T = @floatFromInt(canvas.height); const wf: T = @floatFromInt(canvas.width);

            for (0..canvas.height) |i| for (0..canvas.width) |j| {

                const im = -self.center.im + (3.0 * (@as(T, @floatFromInt(i)) + 0.5) - 1.5 * hf) / self.zoom / hf;
                const re =  self.center.re + (3.0 * (@as(T, @floatFromInt(j)) + 0.5) - 1.5 * wf) / self.zoom / hf;

                const p = Complex(T).init(re, im); var iter: u32 = 0;

                var z, var zp = fractal.init(p);

                while (z.squaredMagnitude() <= bailout * bailout and iter < maxiter) : (iter += 1) {
                    const zt = z; z = fractal.iterate(p, z, zp); zp = zt;
                }

                if (iter < maxiter) {

                    var value: T = @floatFromInt(iter);

                    if (smooth) value -= std.math.log2(std.math.log2(z.magnitude()));

                    canvas.ptr(i, j).* = coloring(self, value / @as(T, @floatFromInt(maxiter)));
                }
            };
        }

        /// Orbitrap algorithm.
        pub fn orbitrap(self: @This(), canvas: *Image, fractal: anytype, coloring: anytype, bailout: T, maxiter: u32, fill: bool) void {
            const hf: T = @floatFromInt(canvas.height); const wf: T = @floatFromInt(canvas.width);

            const trap = switch (self.algorithm.orbitrap.trap) {
                .circle => &@This().circle,
                .disk => &@This().disk,
                .linear => &@This().linear,
                .minabs => &@This().minabs,
                .mindc => &@This().mindc
            };

            for (0..canvas.height) |i| for (0..canvas.width) |j| {

                const im = -self.center.im + (3.0 * (@as(T, @floatFromInt(i)) + 0.5) - 1.5 * hf) / self.zoom / hf;
                const re =  self.center.re + (3.0 * (@as(T, @floatFromInt(j)) + 0.5) - 1.5 * wf) / self.zoom / hf;

                const p = Complex(T).init(re, im); var iter: u32 = 0;

                var z, var zp = fractal.init(p); var dist = std.math.floatMax(T);

                while (z.squaredMagnitude() <= bailout * bailout and iter < maxiter) : (iter += 1) {

                    const zt = z; z = fractal.iterate(p, z, zp); zp = zt; 

                    if (trap(z) < dist) dist = z.magnitude();
                }

                if (iter < maxiter or !fill) {

                    const value: T = switch (self.coloring) {
                        .gradient => 1.0 / (1 + 5 * dist),
                        .periodic => 0.1 * std.math.log10(dist),
                        .solid => dist
                    };

                    canvas.ptr(i, j).* = coloring(self, value);
                }
            };
        }

        /// Gradient coloring function.
        pub fn gradient(self: @This(), value: T) RGB {
            const coloring = self.coloring.gradient;

            const r = @as(u8, @intFromFloat(@as(T, @floatFromInt(coloring.start.r)) + value * (@as(T, @floatFromInt(coloring.end.r)) - @as(T, @floatFromInt(coloring.start.r)))));
            const g = @as(u8, @intFromFloat(@as(T, @floatFromInt(coloring.start.g)) + value * (@as(T, @floatFromInt(coloring.end.g)) - @as(T, @floatFromInt(coloring.start.g)))));
            const b = @as(u8, @intFromFloat(@as(T, @floatFromInt(coloring.start.b)) + value * (@as(T, @floatFromInt(coloring.end.b)) - @as(T, @floatFromInt(coloring.start.b)))));

            return RGB{.r = r, .g = g, .b = b};
        }

        /// Priodic coloring function.
        pub fn periodic(self: @This(), value: T) RGB {
            const coloring = self.coloring.periodic;

            const r = @as(u8, @intFromFloat((std.math.sin(coloring.amplitude[0] * value + coloring.frequency[0]) + 1) * 127.5));
            const g = @as(u8, @intFromFloat((std.math.sin(coloring.amplitude[1] * value + coloring.frequency[1]) + 1) * 127.5));
            const b = @as(u8, @intFromFloat((std.math.sin(coloring.amplitude[2] * value + coloring.frequency[2]) + 1) * 127.5));

            return RGB{.r = r, .g = g, .b = b};
        }

        /// Solid coloring function.
        pub fn solid(self: @This(), _: T) RGB {
            const coloring = self.coloring.solid;

            return coloring.rgb;
        }

        /// Circular trap for orbitrap algorithm.
        pub fn circle(z: Complex(T)) T {
            return @abs(z.magnitude() - 1);
        }

        /// Disk trap for orbitrap algorithm.
        pub fn disk(z: Complex(T)) T {
            return z.magnitude();
        }

        /// Linear trap for orbitrap algorithm.
        pub fn linear(z: Complex(T)) T {
            return 0.70710678 * @min(@abs(z.re - z.im), @abs(z.re + z.im));
        }

        /// Minimum of absolute values trap for orbitrap algorithm.
        pub fn minabs(z: Complex(T)) T {
            return @min(@abs(z.re), @abs(z.im));
        }

        /// Minimum between disk and circle trap for orbitrap algorithm.
        pub fn mindc(z: Complex(T)) T {
            return @min(disk(z), circle(z));
        }
    };
}
