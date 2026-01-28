//! Main file for the fractal generator application.

const std = @import("std");

const device_write = @import("device_write.zig");
const image = @import("image.zig");
const orbit_fractal = @import("orbit_fractal.zig");
const rgb = @import("rgb.zig");

const exportImageAsPPM = device_write.exportImageAsPPM;
const print = device_write.print;
const printJson = device_write.printJson;

const Complex = std.math.Complex;
const Image = image.Image;
const OrbitFractalGenerator = orbit_fractal.OrbitFractalGenerator;
const RGB = rgb.RGB;

/// The options for the fractal generator.
pub fn Options(comptime T: type) type {
    return struct {
        pub const Category = union(enum) {
            orbit: struct {
                algorithm: union(enum) {
                    escape: orbit_fractal.Escape(T),
                } = .{.escape = .{}},
                coloring: union(enum) {
                    solid: orbit_fractal.Solid(T),
                    gradient: orbit_fractal.Gradient(T),
                } = .{.gradient = .{}},
                fractal: union(enum) {
                    mandelbrot: struct{}
                },
                center: [2]T = .{0, 0},
                zoom: T = 1,
                smooth: bool = true,
            }
        };

        background: [3]u8 = .{0, 0, 0},
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

    output.canvas.fill(RGB{.r = opt.background[0], .g = opt.background[1], .b = opt.background[2]});

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
