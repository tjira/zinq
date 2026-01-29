//! Main file for the fractal generator application.

const std = @import("std");

const device_write = @import("device_write.zig");
const image = @import("image.zig");
const orbit_fractal = @import("orbit_fractal.zig");
const orbit_fractal_generator = @import("orbit_fractal_generator.zig");
const rgb = @import("rgb.zig");

const exportImageAsPPM = device_write.exportImageAsPPM;
const print = device_write.print;
const printJson = device_write.printJson;

const Complex = std.math.Complex;
const Image = image.Image;
const OrbitFractalGenerator = orbit_fractal_generator.OrbitFractalGenerator;
const RGB = rgb.RGB;

/// The options for the fractal generator.
pub fn Options(comptime T: type) type {
    return struct {
        pub const Category = union(enum) {
            orbit: struct {
                pub const Algorithm = union(enum) {
                    escape: struct{
                        bailout: T = 8,
                        maxiter: u32 = 64,
                    }
                };
                pub const Coloring = union(enum) {
                    solid: struct {
                        rgb: RGB = .{.r = 255, .g = 255, .b = 255},
                        seed: ?u32 = null
                    },
                    gradient: struct {
                        start: RGB = .{.r =   0, .g =   0, .b =   0},
                        end:   RGB = .{.r = 255, .g = 255, .b = 255},
                        seed: ?u32 = null
                    },
                    periodic: struct {
                        amplitude: [3]T = .{31.93, 30.38, 11.08},
                        frequency: [3]T = .{ 6.26,  5.86,  0.80},
                        seed: ?u32 = null
                    },
                };

                fractal: union(enum) {
                    buffalo: struct{},
                    burningship: struct{},
                    julia: struct {c: Complex(T) = Complex(T).init(0, 1)},
                    mandelbrot: struct{},
                    manowar: struct{},
                    phoenix: struct{c: Complex(T) = Complex(T).init(0, 0)}
                },

                algorithm: Algorithm = .{.escape = .{}},
                coloring: Coloring = .{.periodic = .{}},
                center: Complex(T) = Complex(T).init(0, 0),
                zoom: T = 1,
                smooth: bool = true,
            }
        };

        background: RGB = .{.r = 0, .g = 0, .b = 0},
        resolution: [2]u32 = .{1920, 1080},
        output: []const u8 = "fractal.ppm",
        category: Category,
    };
}

/// Output of the fractal generator.
pub fn Output(comptime _: type) type {
    return struct {
        canvas: Image,

        /// Initialize the output structure.
        pub fn init(height: usize, width: usize, allocator: std.mem.Allocator) !@This() {
            return @This(){
                .canvas = try Image.init(height, width, allocator),
            };
        }

        /// Deinitialize the output.
        pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
            self.canvas.deinit(allocator);
        }
    };
}

/// Run the fractal generator with the given options.
pub fn run(comptime T: type, opt: Options(T), enable_printing: bool, allocator: std.mem.Allocator) !Output(T) {
    if (enable_printing) try printJson(opt);

    var timer = try std.time.Timer.start();

    if (enable_printing) try print("\nINITIALIZING CANVAS:", .{});

    var output = try Output(T).init(opt.resolution[1], opt.resolution[0], allocator);

    output.canvas.fill(opt.background);

    if (enable_printing) try print(" {D}\nGENERATING FRACTALS:", .{timer.read()}); timer.reset();

    switch (opt.category) {
        .orbit => |orbit| {
            OrbitFractalGenerator(T).init(orbit).paint(&output.canvas);
        }
    }

    if (enable_printing) try print(" {D}\nWRITING IMG TO DISK:", .{timer.read()}); timer.reset();

    try exportImageAsPPM(opt.output, output.canvas);

    if (enable_printing) try print(" {D}\n", .{timer.read()});

    return output;
}
