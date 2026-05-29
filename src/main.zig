const std = @import("std");
const writer = @import("writer.zig");

const printf = writer.printf;

pub fn main(init: std.process.Init) !void {
    // const gpa = init.gpa;
    const io = init.io;

    try printf(io, "Hello World!\n", .{});
}
