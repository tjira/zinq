//! Image class.

const std = @import("std");

const rgb = @import("rgb.zig");

const RGB = rgb.RGB;

/// Image class
pub const Image = struct {
    data: []u8, height: usize, width: usize,

    /// Initialize the image.
    pub fn init(height: usize, width: usize, allocator: std.mem.Allocator) !@This() {
        return @This(){.data = try allocator.alloc(u8, height * width * 3), .height = height, .width = width};
    }

    /// Deinitialize the image.
    pub fn deinit(self: *@This(), allocator: std.mem.Allocator) void {
        allocator.free(self.data);
    }

    /// Get the RGB value at the given row and column.
    pub fn at(self: @This(), row: usize, col: usize) RGB {
        const index = (row * self.width + col) * 3;

        return RGB{
            .r = self.data[index + 0],
            .g = self.data[index + 1],
            .b = self.data[index + 2],
        };
    }

    /// Returns a pointer to the RGB value at the given row and column.
    pub fn ptr(self: @This(), row: usize, col: usize) *RGB {
        const index = (row * self.width + col) * 3;

        return @ptrCast(&self.data[index]);
    }

    /// Fill the image with the given RGB color.
    pub fn fill(self: *@This(), color: RGB) void {
        for (0..self.height) |i| for (0..self.width) |j| {
            self.ptr(i, j).* = color;
        };
    }
};
