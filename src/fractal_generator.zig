//! Main file for the fractal generator application.

const std = @import("std");

const device_write = @import("device_write.zig");
const image = @import("image.zig");
const rgb = @import("rgb.zig");

const exportImageAsPPM = device_write.exportImageAsPPM;
const print = device_write.print;
const printJson = device_write.printJson;

const Complex = std.math.Complex;
const Image = image.Image;
const RGB = rgb.RGB;

/// The options for the fractal generator.
pub fn Options(comptime T: type) type {
    return struct {
        pub const Category = union(enum) {
            orbit: struct {
                fractal: union(enum) {
                    mandelbrot: struct{}
                },
                center: [2]T = .{0, 0},
                zoom: T = 1,
                maxiter: u32 = 100,
                escape: T = 2,
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

    const center = Complex(T).init(opt.category.orbit.center[0], opt.category.orbit.center[1]); const zoom: T = 1;

    if (enable_printing) try print(" {D}\nGENERATING FRACTALS:", .{timer.read()}); timer.reset();

    for (0..output.canvas.height) |i| for (0..output.canvas.width) |j| {

        const hf: T = @floatFromInt(output.canvas.height); const wf: T = @floatFromInt(output.canvas.width);

        const im = -center.im + (3.0 * (@as(T, @floatFromInt(i)) + 0.5) - 1.5 * hf) / zoom / hf;
        const re =  center.re + (3.0 * (@as(T, @floatFromInt(j)) + 0.5) - 1.5 * wf) / zoom / hf;

        const p = Complex(T).init(re, im); var z = Complex(T).init(0, 0);

        var iter: u32 = 0;

        while (z.squaredMagnitude() <= opt.category.orbit.escape * opt.category.orbit.escape and iter < opt.category.orbit.maxiter) : (iter += 1) {
            z = z.mul(z); z = z.add(p);
        }

        if (iter < opt.category.orbit.maxiter) {
            output.canvas.ptr(i, j).* = RGB{.r = 255, .g = 255, .b = 255};
        }
    };

    if (enable_printing) try print(" {D}\nWRITING IMG TO DISK:", .{timer.read()}); timer.reset();

    try exportImageAsPPM(opt.output, output.canvas);

    if (enable_printing) try print(" {D}\n", .{timer.read()});

    return output;
}
